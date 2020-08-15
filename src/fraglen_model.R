# Model based on the statistical distance between fragment-length distribution
#
# Calculate statistical distance between two binned segments.
#
# Each bin contains a number of reads aligned to the reference genome. The length distribution
# may be used to infer information about chromosome 3D conformation. To achieve this, we start
# from calculating the statistical distance between any two binned segments


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
cucconi_distance <- function(s1, s2, min_samples = 100) {
  source(here::here("src/nonpar/cucconi.test.R"))
  source(here::here("src/nonpar/cucconi.teststat.R"))
  source(here::here("src/nonpar/cucconi.dist.perm.R"))
  source(here::here("src/nonpar/cucconi.dist.boot.R"))
  # We reuqire at least a certain minimal number of samples
  if (min(length(s1), length(s2)) < min_samples)
    NA
  else {
    if (length(s1) > 10000)
      s1 <- sample(s1, 10000)
    if (length(s2) > 10000)
      s2 <- sample(s2, 10000)
    
    score_cap <- -log10(2e-16)
    min(-log10(cucconi.test(s1, s2)$p.value), score_cap)
  }
}


# Calculate the distance matrix between paired intervals
fraglen_model <- function(bfp, gr, bin_size, interval1, interval2, metrics, ...) {
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
  
  # Determine the distance function
  distance_func <- switch(metrics,
                          "ks" = ks_distance,
                          "cucconi" = cucconi_distance,
                          stop(str_interp("Invalid metrics: ${metrics}")))
  
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
    
    frag1 <- bfp[bin_idx1,]$frag[[1]][[1]]
    bin_start1 <- bfp[bin_idx1,]$bin_start
    frag2 <- bfp[bin_idx2,]$frag[[1]][[1]]
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