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
  
  # score_cap <- -log10(2e-16)
  # min(-log10(ks.test(s1, s2)$p.value), score_cap)
  ks.test(s1, s2)$p.value
}


# Perform a Two-sample K-S test for statistical distance calculation
cucconi_distance <- function(s1, s2, min_samples = 100) {
  # We reuqire at least a certain minimal number of samples
  if (min(length(s1), length(s2)) < min_samples)
    return(NA)
  
  tpepler.nonpar::cucconi.test(s1, s2)$p.value
}


# Calculate the distance matrix between paired intervals
fraglen_model <- function(bfp, gr, bin_size, interval1, interval2, 
                          metrics = "ks", bootstrap = 1, 
                          subsample = NULL, ...) {
  args <- list(...)
  logger_name <- if (is.null(args$logger)) "worker" else args$logger
  logger <- getLogger(name = logger_name)
  
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
    
    pvalues<- 1:bootstrap %>% map_dbl(function(iter) {
      if (iter == 1 || iter %% 50 == 0)
        logdebug(str_interp(
          "Bootstrap #$[d]{iter}: ${chr_name}:$[d]{bin_start1 + 1}-$[d]{bin_start1 + bin_size} vs. ${chr_name}:$[d]{bin_start2 + 1}-$[d]{bin_start2 + bin_size}"), 
          logger = logger_name)
      
      min_sample_cnt <- 10
      if (list(frag1, frag2) %>% map_int(length) %>% min() < min_sample_cnt)
        NA
      else if (!is.null(subsample)) {
        # Only use fix_frag_cnt fragments
        frag1 <- sample(frag1, subsample, replace = TRUE)
        frag2 <- sample(frag2, subsample, replace = TRUE)
        distance_func(frag1, frag2)
      } else {
        distance_func(frag1, frag2)
      }
    })
    # scores <- pvalues %>% map_dbl(function(v) -log10(max(c(2e-16, v))))
    
    # do.call(rbind, 1:2 %>% map(function(idx) { list(start1 = idx, start2 = idx+1, score = runif(3)) %>% as_tibble() }))
    list(
      start1 = bin_start1,
      start2 = bin_start2,
      score = -log10(max(c(2e-16, median(pvalues, na.rm = TRUE)))),
      score2 = -log10(pmax(2e-16, pvalues)),
      pvalue = pvalues,
      bootstrap = 1:bootstrap,
      frag_cnt1 = ifelse(is.null(subsample), length(frag1), subsample),
      frag_cnt2 = ifelse(is.null(subsample), length(frag2), subsample)
    ) %>% as_tibble()
  })
  do.call(rbind, results)
  # # Build the data frame
  # as_tibble(matrix(unlist(results), ncol = 3, byrow = TRUE)) %>%
  #   transmute(start1 = as.integer(V1), start2 = as.integer(V2), score = V3)
}