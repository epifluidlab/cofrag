# Calculate statistical distance between two binned segments.
#
# Each bin contains a number of reads aligned to the reference genome. The length distribution
# may be used to infer information about chromosome 3D conformation. To achieve this, we start
# from calculating the statistical distance between any two binned segments

library(dplyr)
library(doParallel)
library(rtracklayer)
library(Repitools)
library(Rsamtools)
library(stringr)
library(foreach)
library(doParallel)

library(logging)
basicConfig()

calc_bfp <-
  function(bam_file, gr, bin_size) {
    # Load the aligned reads from bam_file, and calculate binned fragmentation profile
    # Parameters:
    #   * bam_file: must be sorted and indexed
    #   * gr: a GRanges object indicating the region of interest
    #   * bin_size
    #
    # Return: list(c(162, 163), c(122), ...)
    
    # Argument check
    stopifnot(width(gr) %% bin_size == 0)
    
    # Load aligned reads
    aligned_reads <-
      scanBam(bam_file,
              param = ScanBamParam(which = gr,
                                   what = c("pos", "isize")))[[1]]
    
    # Distribute the aligned reads into bins
    endpoints <- seq(start(gr), end(gr), by = bin_size)
    
    # For example
    # [,1]     [,2]     [,3]     [,4]     [,5]
    # [1,] 16050001 16050021 16050041 16050061 16050081
    # [2,] 16050020 16050040 16050060 16050080 16050100
    bin_layout <- sapply(1:length(endpoints), function(idx) {
      start_pos <- endpoints[idx]
      if (idx == length(endpoints))
        end_pos <- end(gr)
      else
        end_pos <- endpoints[idx + 1] - 1
      c(start_pos, end_pos)
    })
    
    # =====
    bin_idx <- 1
    
    # The stop points determine how to assign reads to bins
    # For example:
    #      [,1] [,2] [,3]
    # [1,]    1  NA  281
    # [2,]  280  NA  587
    # Means: aligned_reads[1:280] => bin#1, aligned_reads[281:432] => bin#2
    binned_reads_stop_points <-
      matrix(rep(NA, 2 * ncol(bin_layout)),
             nrow = 2,
             ncol = ncol(bin_layout))
    
    for (idx in 1:length(aligned_reads$pos)) {
      # Scan the aligned reads (sorted bam) to calcuate the binning
      pos <- aligned_reads$pos[idx]
      if (between(pos, bin_layout[1, bin_idx], bin_layout[2, bin_idx])) {
        if (idx == 1)
          binned_reads_stop_points[1, bin_idx] <- 1
        else if (idx == length(aligned_reads$pos))
          binned_reads_stop_points[2, bin_idx] <- idx
        next
      }
      
      # pos falls beyond current bin. Finalize the bin.
      if (idx != 1)
        binned_reads_stop_points[2, bin_idx] <- idx - 1
      
      # Try to find where it falls, by incrementally
      # scan the bin_layout matrix
      while (TRUE) {
        bin_idx <- bin_idx + 1
        stopifnot(bin_idx <= ncol(bin_layout))
        if (between(pos, bin_layout[1, bin_idx], bin_layout[2, bin_idx])) {
          binned_reads_stop_points[1, bin_idx] <- idx
          break
        }
      }
      if (idx == length(aligned_reads$pos))
        binned_reads_stop_points[2, bin_idx] <- idx
    }
    
    apply(binned_reads_stop_points, 2, function(x) {
      if (is.na(x[1]) || is.na(x[2]))
        integer(0)
      else
        abs(aligned_reads$isize[x[1]:x[2]])
    })
    
    # =====
    
    # Another way of doing the task.
    # apply(bin_layout, 2, function(x) {
    #   # x: the genomic range: (16050401, 16050420)
    #   abs(aligned_reads$isize[aligned_reads$pos >= x[1] &
    #                             aligned_reads$pos <= x[2]])
    # })
  }

# Perform a Two-sample K-S test for statistical distance calculation
ks_distance <- function(s1, s2) {
  # We reuqire at least 100 samples for both collections
  min_vol <- 100
  if (min(length(s1), length(s2)) < min_vol)
    NA
  else {
    pv <- ks.test(s1, s2)$p.value
    # pv shouldn't be zero
    min_pv <- 5e-320
    if (pv < min_pv)
      pv <- min_pv
    - log10(pv)
  }
  
}

# Cut the genomic range a:b from bfp
clip_bfp <- function(bfp, a, b, gr_start, bin_size) {
  stopifnot((a - gr_start) %% bin_size == 0)
  stopifnot((b + 1 - gr_start) %% bin_size == 0)
  bfp[((a - gr_start) / bin_size + 1):((b + 1 - gr_start) / bin_size)]
}

calc_distance_by_block <-
  function(bfp_pair, stat_distance_func) {
    # Calculate the distance matrix within one block
    # Parameters:
    #   * bfp_pair: a pair of subregion binned fragmentation profile, corresponding
    #       to region_i and region_j respectively. Each element is a list of integer vectors, such as:
    #       list(c(167, 125), c(134, 159, 147)). Each integer represents the insert size of a fragment.
    #   * stat_distance_func: a function to calculate the statistical distance between two samples
    
    # Argument check
    stopifnot(length(bfp_pair) == 2)
    
    # For each pair of bins, perform the KS-test and use the p-value is the
    # statistical distance
    logdebug("Start calculating the distance matrix for the block")
    
    bin_nums <- c(length(bfp_pair[[1]]), length(bfp_pair[[2]]))
    # Create a permutation of all possible bin pairs
    bin_pairs <- expand.grid(1:bin_nums[1], 1:bin_nums[2])
    dist_vec <-
      sapply(1:nrow(bin_pairs), function(l) {
        idx_i <- bin_pairs[l, ]$Var1
        idx_j <- bin_pairs[l, ]$Var2
        
        v1 <- bfp_pair[[1]][[idx_i]]
        v2 <- bfp_pair[[2]][[idx_j]]
        
        # Perform a Two-sample statistical distance calculation
        stat_distance_func(v1, v2)
      })
    
    logdebug("Completed calculating the distance matrix for the block")
    
    dist_vec
  }

calc_distance <-
  function(bam_file,
           gr,
           bin_size,
           block_size,
           nthreads = 1) {
    # Calculate the distance matrix
    # Parameters:
    # * bam_file: must be indexed. Alignd to a single chromosome.
    # * gr: GRanges object stating the genomic region of interest
    # * bin_zize
    # * block_size: the calculation is performed based on "blocks". Notice that
    #     block_size should be multiples of bin_size
    
    # Argument check
    # Both the gr range and block_size should be multiples of bin_size
    gr <- GRanges(gr)
    stopifnot(width(gr) %% bin_size == 0)
    stopifnot(block_size %% bin_size == 0)
    stopifnot(block_size <= width(gr))
    
    # Divide the genome of interest into blocks
    gr_start <- start(gr)
    gr_end <- end(gr)
    block_num <- round((gr_end - gr_start) / block_size)
    
    loginfo(
      sprintf(
        "Calculation summary: BAM: %s Range: %s:%d-%d, bin_size=%d, block_size=%d, nthreads=%d",
        bam_file,
        seqnames(gr),
        start(gr),
        end(gr),
        bin_size,
        block_size,
        nthreads
      )
    )
    
    loginfo(sprintf(
      "Region-of-interest has been divided into %d sub-regions",
      block_num
    ))
    
    # Example: c(IRanges(1, 100), IRanges(101, 200))
    block_layout <- lapply(1:block_num, function(i) {
      if (i < block_num)
        rg <- IRanges(start = gr_start + (i - 1) * block_size,
                      end = gr_start + i * block_size - 1)
      else
        rg <-
          IRanges(start = gr_start + (i - 1) * block_size, end = gr_end)
      loginfo(sprintf(
        "Sub-region #%d: %s:%d-%d, width: %d",
        i,
        seqnames(gr),
        start(rg),
        end(rg),
        width(rg)
      ))
      rg
    })
    
    # Calculate the binned fragmentation profile for the entire genomic range
    loginfo("Calculating the binned fragmentation profile")
    bfp <- calc_bfp(bam_file, gr, bin_size)
    
    # Calculate the distance matrixes block by block. Row first.
    # dist_matrix_list <- list()
    
    block_pairs <- list()
    block_idx <- 0
    for (row_idx in 1:block_num) {
      for (col_idx in row_idx:block_num) {
        block_idx <- block_idx + 1
        block_pairs[[block_idx]] <- c(row_idx, col_idx)
      }
    }
    
    loginfo("Calculating the distance matrix")

    gr_name = as.character(seqnames(gr))
    cl <- makeCluster(nthreads)
    registerDoParallel(cl)
    dist_matrix_list <-
      foreach (
        pair = block_pairs,
        .combine = "list",
        .multicombine = TRUE,
        .export = c(
          "loginfo",
          "logdebug",
          "calc_distance_by_block",
          "clip_bfp",
          "ks_distance",
          "start",
          "end"
        )
      ) %dopar% {
        row_idx <- pair[1]
        col_idx <- pair[2]
        
        bfp_pair <-
          list(
            clip_bfp(
              bfp,
              start(block_layout[[row_idx]]),
              end(block_layout[[row_idx]]),
              gr_start,
              bin_size
            ),
            clip_bfp(
              bfp,
              start(block_layout[[col_idx]]),
              end(block_layout[[col_idx]]),
              gr_start,
              bin_size
            )
          )
        
        calc_distance_by_block(bfp_pair, ks_distance)
      }
    
    stopCluster(cl)
    
    loginfo("Completed calculating the distance matrix")
    
    dist_matrix_list
  }