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
library(logging)
library(foreach)
library(doParallel)

basicConfig()

split_interval <- function(a, b, n) {
  # Split the interval [a, b] into n equal sub-intervals
  #
  # It returns a matrix:
  #
  # > split_interval(1, 87, 7)
  #
  # [1]  1 13 26 38 50 62 75 87
  # [,1] [,2] [,3] [,4] [,5] [,6] [,7]
  # [1,]    1   13   26   38   50   62   75
  # [2,]   12   25   37   49   61   75   87
  stop_points <- round(seq(a, b, by = (b - a) / n))
  sapply(1:n, function(i) {
    if (i < n - 1) {
      c(stop_points[i], stop_points[i + 1] - 1)
    } else {
      # The last interval
      c(stop_points[i], stop_points[i + 1])
    }
  })
}

get_fragment_len <- function(bam_file, grange) {
  which <- GRanges(grange)
  what <- c("isize")
  param <- ScanBamParam(which = which, what = what)
  scanBam(bam_file, param = param)
}

calc_binned_insert_lengths <-
  function(bam_file, gr, bin_size) {
    # Load the aligned reads from bam_file, and calculate binned insert lengths
    # Parameters:
    #   * gr: a GRanges object.
    #
    # Return: list(c(162, 163), c(122), ...)
    
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
    
    # apply(bin_layout, 2, function(x) {
    #   # x: the genomic range: (16050401, 16050420)
    #   abs(aligned_reads$isize[aligned_reads$pos >= x[1] &
    #                             aligned_reads$pos <= x[2]])
    # })
  }

calc_distance_by_block <-
  function(bfp_pair) {
    # Calculate the distance matrix within one block
    # Parameters:
    #   * bam_file: must be indexed.
    #   * gr: GRanges objects stating the genomic regions this block matches.
    #   * chr_name: e.g. 22
    #   * block_range: the starting and ending genomic coordinates. For example,
    #       the following is a 5kb block: list(c(200000, 205000), c(300000, 305000)).
    #       The interval is inclusive on both ends.
    #   * bfp_pair: a pair of subregion binned fragmentation profile, corresponding
    #       to region_i and region_j respectively.
    
    # Argument check
    # Two regions should be of the same chromosome
    # stopifnot(length(gr) == 2)
    # stopifnot(identical(seqnames(gr)[1], seqnames(gr)[2]))
    stopifnot(length(bfp_pair) == 2)
    
    # # Two blocks should be of the same size
    # stopifnot(width(gr)[1] == width(gr)[2])
    # block_size must be a multiple of bin_size
    # stopifnot(width(gr)[1] %% bin_size == 0)
    
    # loginfo("Calculate binned insert lengths")
    
    # # Aligned reads will be put into bins
    # # It has two objects, each one a numeric vector representing the insert lengths
    # # [[1]]: c(167, 157, 163, ...)
    # binned_insert_lengths <- list()
    # for (i in 1:2) {
    #   # Avoid duplicate loading for diagonal blocks
    #   if (i == 2 && identical(gr[1], gr[2])) {
    #     binned_insert_lengths[[2]] <- binned_insert_lengths[[1]]
    #     break
    #   }
    #
    #   binned_insert_lengths[[i]] <-
    #     calc_binned_insert_lengths(bam_file, gr[i], bin_size)
    # }
    
    # loginfo("Completed calculating binned insert lengths")
    
    # For each pair of bins, perform the KS-test and use the p-value is the
    # statistical distance
    
    # library(foreach)
    # library(doParallel)
    
    logdebug("Start calculating the distance matrix for the block")
    
    # cl <- makeCluster(nthreads)
    # registerDoParallel(cl)
    
    bin_nums <- c(length(bfp_pair[[1]]), length(bfp_pair[[2]]))
    # Create a permutation of all possible bin pairs
    bin_pairs <- expand.grid(1:bin_nums[1], 1:bin_nums[2])
    dist_vec <-
      sapply(1:nrow(bin_pairs), function(l) {
        # foreach(l = 1:nrow(bin_pairs), .combine = c) %dopar% {
        idx_i <- bin_pairs[l, ]$Var1
        idx_j <- bin_pairs[l, ]$Var2
        
        v1 <- bfp_pair[[1]][[idx_i]]
        v2 <- bfp_pair[[2]][[idx_j]]
        
        # We need at least samples to apply K-S test
        min_vol <- 10
        if (min(length(v1), length(v2)) < min_vol)
          NA
        else {
          pv <- ks.test(v1, v2)$p.value
          # pv shouldn't be zero
          min_pv <- 10 ^ (-10)
          if (pv < min_pv)
            pv <- min_pv
          - log10(pv)
        }
      })
    # stopCluster(cl)
    
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
    bfp <- calc_binned_insert_lengths(bam_file, gr, bin_size)
    
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
    cl <- makeCluster(nthreads)
    registerDoParallel(cl)
    
    # Cut the genomic range a:b from bfp
    clip_bfp <- function(a, b) {
      stopifnot((a - gr_start) %% bin_size == 0)
      stopifnot((b + 1 - gr_start) %% bin_size == 0)
      bfp[((a - gr_start) / bin_size + 1):((b + 1 - gr_start) / bin_size)]
    }
    
    gr_name = as.character(seqnames(gr))
    dist_matrix_list <-
      foreach (
        pair = block_pairs,
        .combine = "list",
        .multicombine = TRUE,
        .export = c("logdebug", "loginfo", "calc_distance_by_block")
      ) %dopar% {
        row_idx <- pair[1]
        col_idx <- pair[2]
        
        loginfo(
          sprintf(
            "Process Block (%d, %d), %s:%d-%d vs. %s:%d-%d",
            row_idx,
            col_idx,
            gr_name,
            start(block_layout[[row_idx]]),
            end(block_layout[[row_idx]]),
            gr_name,
            start(block_layout[[col_idx]]),
            end(block_layout[[col_idx]])
          )
        )
        
        bfp_pair <-
          list(clip_bfp(start(block_layout[[row_idx]]), end(block_layout[[row_idx]])),
               clip_bfp(start(block_layout[[col_idx]]), end(block_layout[[col_idx]])))
        
        calc_distance_by_block(bfp_pair)
        
        #
        # gr_block <-
        #   GRanges(seqnames(gr), IRanges(
        #     start = c(start(block_layout[[row_idx]]), start(block_layout[[col_idx]])),
        #     end = c(end(block_layout[[row_idx]]), end(block_layout[[col_idx]]))
        #   ))
        #
        #
        # loginfo(sprintf("Process Block (%d, %d)", row_idx, col_idx))
        # dist_matrix_list[[block_id]] <-
        #   calc_distance_by_block(bam_file, gr_block, bin_size, binned_insert_lengths, nthreads = nthreads)
      }
    
    stopCluster(cl)
    
    
    
    # block_id <- 0
    # for (row_idx in 1:block_num) {
    #   for (col_idx in row_idx:block_num) {
    #     block_id <- block_id + 1
    #     gr_block <-
    #       GRanges(seqnames(gr), IRanges(
    #         start = c(start(block_layout[[row_idx]]), start(block_layout[[col_idx]])),
    #         end = c(end(block_layout[[row_idx]]), end(block_layout[[col_idx]]))
    #       ))
    #     loginfo(sprintf("Process Block (%d, %d)", row_idx, col_idx))
    #     dist_matrix_list[[block_id]] <-
    #       calc_distance_by_block(bam_file, gr_block, bin_size, binned_insert_lengths, nthreads = nthreads)
    #   }
    # }
    
    dist_matrix_list
    loginfo("Completed calculating the distance matrix")
  }


process_arguments <- function() {
  # Default value
  arguments <- list(nthreads = 1)
  args_list <- commandArgs(trailingOnly = TRUE)
  idx <- 1
  while (idx <= length(args_list)) {
    action <- args_list[idx]
    if (action == "--grange") {
      idx <- idx + 1
      arguments$grange <- args_list[idx]
    } else if (action == "--bin-size") {
      idx <- idx + 1
      arguments$bin_size <- as.numeric(args_list[idx])
    } else if (action == "--block-size") {
      idx <- idx + 1
      arguments$block_size <- as.numeric(args_list[idx])
    } else if (action == "-o") {
      idx <- idx + 1
      arguments$output_file <- args_list[idx]
    } else if (action == "-n") {
      idx <- idx + 1
      arguments$nthreads <- as.integer(args_list[idx])
    } else if (idx == length(args_list)) {
      arguments$bam_file <- args_list[idx]
      break
    }
    idx <- idx + 1
  }
  arguments
}

main <- function(arguments) {
  dist_matrix_list <-
    calc_distance(
      arguments$bam_file,
      GRanges(arguments$grange),
      arguments$bin_size,
      arguments$block_size,
      arguments$nthreads
    )
  
  save(dist_matrix_list,
       file = arguments$output_file,
       ascii = TRUE)
}

main(process_arguments())