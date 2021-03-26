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
    return(list(score = NA, pvalue = NA))
  
  # score_cap <- -log10(2e-16)
  # min(-log10(ks.test(s1, s2)$p.value), score_cap)
  pv <- ks.test(s1, s2)$p.value
  list(
    score = -log10(2e-16) + log10(max(2e-16, pv)),
    pvalue = pv
  )
}


# Perform a Two-sample K-S test for statistical distance calculation
cucconi_distance <- function(s1, s2, min_samples = 100) {
  # We reuqire at least a certain minimal number of samples
  if (min(length(s1), length(s2)) < min_samples)
    return(list(score = NA, pvalue = NA))
  
  pv <- tpepler.nonpar::cucconi.test(s1, s2)$p.value
  list(
    score = -log10(2e-16) + log10(max(2e-16, pv)),
    pvalue = pv
  )
}


# Calculate the total variance distance
tv_distance <- function(s1, s2, 
                        length_from = 100,
                        length_to = 350,
                        min_samples = 100) {
  if (min(length(s1), length(s2)) < min_samples)
    return(list(score = NA))
  
  pdfs <- list(s1, s2) %>% map(~ {
    s <- .x
    cdf <- ecdf(s)
    data <- length_from:length_to %>% map_dbl(cdf)
    tail(data, n = -1) - head(data, n = -1)
  })
  
  list(score = max(abs(pdfs[[1]] - pdfs[[2]])))
}


ks_stat_distance <- function(s1, s2, 
                        length_from = 100,
                        length_to = 350,
                        min_samples = 100) {
  if (min(length(s1), length(s2)) < min_samples)
    return(list(score = NA))
  
  cdfs <- list(s1, s2) %>% map(~ {
    s <- .x
    cdf <- ecdf(s)
    length_from:length_to %>% map_dbl(cdf)
  })
  
  list(score = max(abs(cdfs[[1]] - cdfs[[2]])))
}

cvm_stat_distance <- function(s1, s2, 
                           length_from = 100,
                           length_to = 350,
                           min_samples = 100) {
  if (min(length(s1), length(s2)) < min_samples)
    return(list(score = NA))
  
  cdfs <- list(s1, s2) %>% map(~ {
    s <- .x
    cdf <- ecdf(s)
    length_from:length_to %>% map_dbl(cdf)
  })
  
  list(score = (cdfs[[1]] - cdfs[[2]]) ** 2 %>% sum())
}

ad_stat_distance <- function(s1, s2, 
                              length_from = 100,
                              length_to = 350,
                              min_samples = 100) {
  if (min(length(s1), length(s2)) < min_samples)
    return(list(score = NA))
  
  cdfs <- list(s1, s2) %>% map(~ {
    s <- .x
    cdf <- ecdf(s)
    length_from:length_to %>% map_dbl(cdf)
  })
  
  a = (cdfs[[1]] - cdfs[[2]]) ** 2
  b = cdfs[[1]] * (1 - cdfs[[2]])
  non_zeros = (b != 0)
  list(score = sum(a[non_zeros] / b[non_zeros]))
  # sum(a)
}

cucconi_stat_distance <- function(s1, s2) {
  x = s1
  y = s2
  m = length(x)
  n = length(y)
  
  N <- m + n
  S <- rank(c(x, y))[(m + 1):N]
  denom <- sqrt(m * n * (N + 1) * (2 * N + 1) * (8 * N + 11) / 5)
  U <- (6 * sum(S^2) - n * (N + 1) * (2 * N + 1)) / denom
  V <- (6 * sum((N + 1 - S)^2) - n * (N + 1) * (2 * N + 1)) / denom
  rho <- (2 * (N^2 - 4)) / ((2 * N + 1) * (8 * N + 11)) - 1
  C <- (U^2 + V^2 - 2 * rho * U * V) / (2 * (1 - rho^2))
  
  list(score = C)
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
                          "cucconi.stat" = cucconi_stat_distance,
                          "tv" = tv_distance,
                          "ks.stat" = ks_stat_distance,
                          "cvm.stat" = cvm_stat_distance,
                          "ad.stat" = ad_stat_distance,
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
  
  1:nrow(paired_bins) %>% map_dfr(function(row_idx) {
    bin_idx1 <- paired_bins[row_idx,]$bin1
    bin_idx2 <- paired_bins[row_idx,]$bin2
    
    frag1 <- bfp[bin_idx1,]$frag[[1]][[1]]
    n_frag1 <- length(frag1)
    bin_start1 <- bfp[bin_idx1,]$bin_start
    frag2 <- bfp[bin_idx2,]$frag[[1]][[1]]
    n_frag2 <- length(frag2)
    bin_start2 <- bfp[bin_idx2,]$bin_start
    
    # Minimal # of fragments
    if (min(n_frag1, n_frag2) < 100)
      return(NULL)
    
    1:bootstrap %>% map_dfr(function(iter) {
      if (iter == 1 || iter %% 5 == 0)
        logdebug(str_interp(
          "Bootstrap #$[d]{iter}: ${chr_name}:$[d]{bin_start1 + 1}-$[d]{bin_start1 + bin_size} vs. ${chr_name}:$[d]{bin_start2 + 1}-$[d]{bin_start2 + bin_size}"), 
          logger = logger_name)
      
      if (!is.null(subsample)) {
        # Only use fix_frag_cnt fragments
        frag1 <- sample(frag1, subsample, replace = TRUE)
        frag2 <- sample(frag2, subsample, replace = TRUE)
      }
      distance_func(frag1, frag2) %>% as_tibble() %>% mutate(
        bootstrap = as.integer(iter),
        frag_cnt1 = n_frag1, # length(frag1),
        frag_cnt2 = n_frag2, # length(frag2)
        frag_ss_cnt1 = length(frag1),
        frag_ss_cnt2 = length(frag2)
        )
    }) %>%
      mutate(
        start1 = bin_start1,
        start2 = bin_start2,
        score2 = score,
        score = median(score2, na.rm = TRUE)
      ) %>%
      select(start1, start2, score, score2, bootstrap, frag_cnt1, frag_cnt2, everything())
  }) %>%
    na.omit()
}