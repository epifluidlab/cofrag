# Calculate statistical distance between two binned segments.
#
# Each bin contains a number of reads aligned to the reference genome. The length distribution
# may be used to infer information about chromosome 3D conformation. To achieve this, we start
# from calculating the statistical distance between any two binned segments

library(foreach)
library(doParallel)
library(logging)
library(tidyverse)

requireNamespace("GenomicRanges")

source(here::here("genomic_matrix.R"))


# Read fragments data from connection
# Return a 
.import_fragments <- function(conn) {
  frag_data <-
    read_delim(
      conn,
      delim = "\t",
      col_names = c("chr", "start", "end", "mapq", "strand"),
      col_types = cols(
        col_factor(),
        col_integer(),
        col_integer(),
        col_integer(),
        col_character()
      )
    )
  
  # Only work for intra-chromosome cases
  chr_names <- unique(frag_data$chr)
  stopifnot(length(chr_names) == 1)
  chr <- as.character(unique(frag_data$chr))[1]
  
  loginfo("Finished importing fragments")
  loginfo(str_interp("Total fragments #: ${nrow(frag_data)}"))
  
  range_start <- head(frag_data, n = 1)$start
  range_end <- tail(frag_data, n = 1)$end
  loginfo(str_interp("Range: ${chr}:${range_start + 1}-${range_end}"))
  
  frag_data %>% mutate(length = end - start)
}


# Calculate the bin index for a certain coordinate (pos)
.calc_bin_index <- function(pos, gr_start, bin_size) {
  as.integer((pos - gr_start) %/% bin_size + 1)
}


# Calculate binned fragmentation profile
# Parameters:
#   * frag_data: BED-format data frame representing fragments
#   * gr: a GRanges object indicating the region of interest
#   * bin_size
#
# Return a data frame containing fragment lengths for each bin:
# # A tibble: 704 x 3
#     bin_idx bin_start frag          
#       <int>     <int> <list>        
#   1       1  16100000 <int [18,163]>
#   2       2  16150000 <int [39,340]>
.calc_bfp <-
  function(frag_data, gr, bin_size) {

    
    stopifnot(GenomicRanges::width(gr) %% bin_size == 0)
    
    # Build the bin layout
    gr_start <- GenomicRanges::start(gr) - 1
    bin_layout <- sapply(1:(GenomicRanges::width(gr) %/% bin_size), function(v) {
      as.integer(gr_start + (v - 1) * bin_size)
    })
    
    # bin1: index of the bin which contains the fragment start
    # bin2: index of the bin which contains the fragment end
    # cross_bin: logical value indicating whether the fragment straddles across adjacent bins
    frag_data <- frag_data %>% mutate(
      length = end - start,
      bin1 = as.integer(.calc_bin_index(start, gr_start, bin_size)),
      bin2 = as.integer(.calc_bin_index(end - 1, gr_start, bin_size)),
      cross_bin = (bin1 != bin2)
      )
    
    frag_coll_1 <- frag_data %>% group_by(bin1) %>% summarize(frag_coll = list(length))
    # For cross-bin fragments, group them by bins. When retrieving fragments in
    # a certain bin, we should not only consider the fragments starting in the
    # bin, but also those ending in the bin. Thus, we need to combine
    # frag_coll_1 and frag_coll_2
    frag_coll_2 <- frag_data %>% filter(cross_bin == TRUE) %>% group_by(bin2) %>% summarize(frag_coll = list(length))
    
    # A helper function which retrieve all fragments for a specified bin
    retrieve_frag_coll <- function(data, bin_name, bin_idx) {
      results <- data %>% filter(.data[[bin_name]] == bin_idx)
      if (nrow(results) > 0)
        (results$frag_coll)[[1]]
      else
        integer(0)
    }
    
    as_tibble(list(
      bin_idx = seq_along(bin_layout), 
      bin_start = bin_layout,
      frag = lapply(seq_along(bin_layout), function(v) {
        coll1 <- retrieve_frag_coll(frag_coll_1, "bin1", v)
        coll2 <- retrieve_frag_coll(frag_coll_2, "bin2", v)
        c(coll1, coll2)
      })
    ))
  }

# Perform a Two-sample K-S test for statistical distance calculation
ks_distance <- function(s1, s2, min_samples = 100) {
  # We reuqire at least a certain minimal number of samples
  # The minimal meaningful p-value is (< 2e-16). The corresponding -log10 value is 15.7
  # Therefore, here we limit the score range to [0, 15.7)
  if (min(length(s1), length(s2)) < min_samples)
    return(NA)
  
  score_cap <- -log10(2e-16)
  min(-log10(ks.test(s1, s2)$p.value), score_cap)
}


# Perform a Two-sample K-S test for statistical distance calculation
cucconi_distance <- function(s1, s2, min_samples = 1000) {
  source(here::here("nonpar/cucconi.test.R"))
  source(here::here("nonpar/cucconi.teststat.R"))
  source(here::here("nonpar/cucconi.dist.perm.R"))
  source(here::here("nonpar/cucconi.dist.boot.R"))
  # We reuqire at least a certain minimal number of samples
  if (min(length(s1), length(s2)) < min_samples)
    NA
  else {
    if (length(s1) > 10000)
      s1 <- sample(s1, 10000)
    if (length(s2) > 10000)
      s2 <- sample(s2, 10000)
    
    pv <- cucconi.test(s1, s2)$p.value
    # pv shouldn't be zero
    min_pv <- 5e-320
    if (pv < min_pv)
      pv <- min_pv
    - log10(pv)
  }
}


# Calculate the distance matrix between paired intervals
.calc_distance_paired_intervals <- function(bfp, gr, bin_size, interval1, interval2, distance_func) {
  # Argument check
  # gr and intervals should be multiples of bin_size
  stopifnot(all(list(gr, interval1, interval2) %>% map_lgl(function(v) {
    GenomicRanges::width(v) %% bin_size == 0
  })))
  # Intervals should be sub-regions of gr
  stopifnot(all(list(interval1, interval2) %>% map_lgl(function(v) {
    gr_start <- GenomicRanges::start(gr) - 1
    gr_end <- GenomicRanges::end(gr)
    interv_start <- GenomicRanges::start(v) - 1
    interv_end <- GenomicRanges::end(v)
    all(between(c(interv_start, interv_end), gr_start, gr_end))
  })))
  # stopifnot(GenomicRanges::width(interval1) == GenomicRanges::width(interval2))
  
  # Bin index list for each interval
  bin_idx <- list(interval1, interval2) %>% map(function(interv) {
    num <- GenomicRanges::width(interv) %/% bin_size
    start_idx <- (GenomicRanges::start(interv) - GenomicRanges::start(gr)) %/% bin_size + 1
    seq(start_idx, start_idx + num - 1)
  })
  
  # Generate a list of paired bins, where bin_idx1 <= bin_idx2 (only process the upper-right half of the matrix)
  paired_bins <- expand.grid(bin_idx[[1]], bin_idx[[2]]) %>%
    transmute(bin1 = Var1, bin2 = Var2) %>%
    filter(bin1 <= bin2)
  
  chr_name <- as.character(GenomicRanges::seqnames(gr))
  
  # # Stats
  # # Total number of fragments in each interval
  # num_frag_interval1 <- bfp[paired_bins$bin1,]$frag %>% map_int(function(x) length(x)) %>% sum
  # num_frag_interval2 <- bfp[paired_bins$bin2,]$frag %>% map_int(function(x) length(x)) %>% sum
  # logdebug(str_interp("# of fragments in each interval: ${num_frag_interval1} vs. ${num_frag_interval2}"), logger = "worker")
  # logdebug(str_interp("# of paired bins: ${nrow(paired_bins)}"), logger = "worker")
  
  results <- 1:nrow(paired_bins) %>% map(function(row_idx) {
    bin_idx1 <- paired_bins[row_idx,]$bin1
    bin_idx2 <- paired_bins[row_idx,]$bin2
    
    frag1 <- bfp[bin_idx1,]$frag[[1]]
    bin_start1 <- bfp[bin_idx1,]$bin_start
    frag2 <- bfp[bin_idx2,]$frag[[1]]
    bin_start2 <- bfp[bin_idx2,]$bin_start
    
    score <- distance_func(frag1, frag2)
    list(
      start1 = bin_start1,
      start2 = bin_start2,
      score = score
    )
  })
  # Build the data frame
  as_tibble(matrix(unlist(results), ncol = 3, byrow = TRUE)) %>%
    transmute(start1 = as.integer(V1), start2 = as.integer(V2), score = V3)
}


# In order to calculate the distance matrix, we first need to infer the genomic
# range from the fragments, and make paddings as necessary so that: the entire 
# genomic range is multiple of bin_size
.make_paddings <- function(frag_data, bin_size) {
  # Infer the genomic range
  chr_name <- frag_data[1,]$chr
  range_start <- frag_data[1,]$start
  range_end <- tail(frag_data, n = 1)$end
  
  # Start is round to multiples of bin_size. Notice that range_start is 0-based and semi-inclusive
  # For example: padded range: 1000-2500 (BED coordinates, bin_size=500)
  logdebug(str_interp("range_start: $[d]{range_start}"))
  padded_start <- range_start %/% bin_size * bin_size
  padded_end <- (range_end - padded_start) %/% bin_size * bin_size + padded_start
  logdebug(str_interp("padded_start: $[d]{padded_start}, bin_size: $[d]{bin_size}, padded_end: $[d]{padded_end}"))
  if (padded_end < range_end)
    padded_end <- padded_end + bin_size
  
  list(start = padded_start, end = padded_end)
}


.build_block_layout <- function(gr, block_size) {
  gr_start <- GenomicRanges::start(gr) - 1
  gr_end <- GenomicRanges::end(gr)
  block_start <- as.integer(seq(gr_start, gr_end, by = block_size) %>% keep(function(x) x < gr_end))
  block_end <- (block_start + block_size) %>% map_int(function(x) {
    as.integer(ifelse(x < gr_end, x, gr_end))
  })
  as_tibble(list(start = block_start, end = block_end))
}


# Divide the genomic range into blocks, and return a list of paired blocks from the upper-right half of the matrix
# Example:
#       i     j start1  end1 start2  end2
#   <int> <int>  <int> <int>  <int> <int>
# 1     1     1      0   300      0   300
# 2     1     2      0   300    300   600
# 3     2     2    300   600    300   600
.build_block_pairs <- function(gr, block_size) {
  gr_start <- GenomicRanges::start(gr) - 1
  gr_end <- GenomicRanges::end(gr)
  num_blocks <- length(seq(gr_start, gr_end, by = block_size) %>% keep(function(x) x < gr_end))
  
  block_layout <- .build_block_layout(gr, block_size)
  
  as_tibble(expand.grid(1:num_blocks, 1:num_blocks) %>% filter(Var2 >= Var1)) %>% 
    transmute(i = Var1, j = Var2, 
              start1 = block_layout$start[i], end1 = block_layout$end[i],
              start2 = block_layout$start[j], end2 = block_layout$end[j])
}


calc_distance <-
  function(frag_data,
           bin_size = 50000,
           block_size = 500000,
           ncores = 1,
           metrics = "ks",
           opts = list()) {
    
    # Calculate the distance matrix
    # Parameters:
    # * frag_data: a BED-format data frame containing fragments
    # * bin_zize
    # * block_size: the calculation is performed based on "blocks". Notice that
    #     block_size should be multiples of bin_size
    # * ncores: multi-processing
    # * opts: other options for the calculation
    
    # Argument check
    # block_size should be multiple of bin_size
    if (block_size %% bin_size != 0)
      stop(str_interp("Block size ${block_size} should be multiple of bin size ${bin_size}"))
    
    padded_range <- .make_paddings(frag_data, bin_size)
    gr <- GenomicRanges::GRanges(str_interp("${frag_data[1,]$chr}:${padded_range$start + 1}-${padded_range$end}"))
    loginfo(str_interp("Genomic range padded to: ${gr}"))
    
    # Calculate the binned fragment profile (BFP)
    bfp <- .calc_bfp(frag_data, gr, bin_size)
    
    block_pairs <- .build_block_pairs(gr, block_size)
    
    interval_pairs <- 1:nrow(block_pairs) %>% map(function(idx) {
      pair <- block_pairs[idx,]
      interval1 <- GenomicRanges::GRanges(GenomicRanges::seqnames(gr), str_interp("${pair[['start1']] + 1}-${pair[['end1']]}"))
      interval2 <- GenomicRanges::GRanges(GenomicRanges::seqnames(gr), str_interp("${pair[['start2']] + 1}-${pair[['end2']]}"))
      list(interval1, interval2)
    })
    
    distance_func <- switch(metrics,
      "ks" = ks_distance,
      "cucconi" = cucconi_distance,
      stop(str_interp("Invalid metrics: ${metrics}"))
    )
    
    loginfo(str_interp("Registering cluster for ${ncores} cores"))
    logger_level <- getLogger()[['level']]
        
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    results <-
      foreach (
        pair_index = 1:length(interval_pairs),
        .combine = "rbind",
        # .multicombine = TRUE,
        .export = c(
          "loginfo",
          "str_interp",
          ".calc_distance_paired_intervals"
        )
      ) %dopar% {
        library(logging)
        library(tidyverse)
       
        # Create a logger for the worker
        logger <- getLogger('worker')
        addHandler(writeToFile, logger = "worker", file = "/dev/stderr")
        setLevel(logger_level, container = "worker")
        
        pair <- interval_pairs[[pair_index]]
        interval1 <- pair[[1]]
        interval2 <- pair[[2]]
        loginfo(str_interp("Processing block #${pair_index} / ${length(interval_pairs)}: ${interval1} vs. ${interval2}"), logger = "worker")
        .calc_distance_paired_intervals(bfp, gr, bin_size, interval1, interval2, distance_func)
      }
    stopCluster(cl)
    
    genomic_matrix(results, frag_data[1,]$chr, padded_range$start, padded_range$end, bin_size)
  }


# Entrypoint for distance calculation
calc_distance_cli <- function(options) {
  loginfo("Command line arguments:")
  loginfo(paste(
    names(options),
    options,
    sep = "=",
    collapse = "; "
  ))
  
  # Import fragmentation data from stdin
  frag_data <-
    .import_fragments(conn = fifo("/dev/stdin")) %>% filter(mapq >= options$`min-mapq` & length <= options$`max-frag-size`)
  
  distance_matrix <- calc_distance(
    frag_data,
    bin_size = options$`bin-size`,
    block_size = options$`block-size`,
    ncores = options$`num-cores`,
    metrics = options$metrics,
    opts = list(min_samples = options$`min-samples`)
  )
  
  dump_genomic_matrix(distance_matrix)
}