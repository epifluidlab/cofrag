calc_disp_matrix <- function(nrow, ncol) {
  # Return the displacement matrix, i.e. abs(row_index - col_index)
  # For example:
  #      [,1] [,2] [,3]
  # [1,]    0    1    2
  # [2,]    1    0    1
  # [3,]    2    1    0
  row_matrix <- matrix(rep(1:nrow, ncol), nrow = nrow, ncol = ncol)
  col_matrix <-
    matrix(rep(1:ncol, 1, each = nrow), nrow = nrow, ncol = ncol)
  abs(row_matrix - col_matrix)
}

# Summarize matrix entries grouped by abs(i - j) distance.
# This function is the foundation of calculating the frequency curve
# and the proximity curve
calc_grouped_summary <- function(data) {
  # data: a square matrix
  # return a data frame which is the summary metrics for different distances
  
  # The matrix should be square
  stopifnot(nrow(data) == ncol(data))
  
  library(dplyr)
  
  disp_matrix <- calc_disp_matrix(nrow(data), ncol(data))
  as_tibble(data.frame(score = c(data), dist = c(disp_matrix))) %>% group_by(dist) %>% summarise(
    mean = mean(score, na.rm = TRUE),
    median = median(score, na.rm = TRUE),
    min = min(score, na.rm = TRUE),
    max = max(score, na.rm = TRUE),
    n = sum(!is.na(score)),
    var = var(score, na.rm = TRUE),
    sum = sum(score, na.rm = TRUE)
  )
}

calc_contact_prob <- function(data, nthread = 1) {
  sm <- calc_grouped_summary(data, nthread)
  sm %>% select()
}

calc_freq_curve <- function(data, bin_size, nthread = 1) {
  # Calculate the contact frequency curve
  # The matrix should be square
  stopifnot(nrow(data) == ncol(data))
  
  disp_matrix <- calc_disp_matrix(nrow(data), ncol(data))
  
  # sapply(0:(nrow(data) - 1), function(d)
  #   mean(data[disp_matrix == d], na.rm = TRUE))
  
  library(foreach)
  library(doParallel)
  
  cl <- makeCluster(nthread)
  registerDoParallel(cl)
  
  tasks_range <-
    round(seq(
      from = 0,
      to = nrow(data),
      length.out = nthread + 1
    ))
  results <-
    foreach(
      task_id = 1:nthread,
      .combine = "c",
      .multicombine = T,
      .export = "str_match"
    ) %dopar% {
      task_start <- tasks_range[task_id]
      task_end <- tasks_range[task_id + 1] - 1
      
      r <- lapply(task_start:task_end, function(d) {
        # All entries with d distance
        data_d_dist <- data[disp_matrix == d]
        data_d_dist <- data_d_dist[!is.na(data_d_dist)]
        data_summary <- summary(data_d_dist, na.rm = T)
        data_summary$n <- length(data_d_dist)
        unlist(data_summary)
      })
      r
    }
  
  stopCluster(cl)
  
  # results is a list of rows. Each row is a numeric vector representing "summary"
  # Here we convert it to a data.frame
  df <- data.frame(matrix(unlist(results), nrow = length(results), byrow = TRUE))
  colnames(df) <- names(results[[1]])
  df[["Dist."]] <- (1:nrow(df) - 1) * bin_size
  df
}

calc_prox_curve <- calc_freq_curve

# Load Hi-C data in _obs_data format and get grouped summary
calc_hic_prob <- function(hic_obs) {
  tot <- sum(hic_obs$contact, na.rm = T)
  hic_obs %>% mutate(dist = abs(i - j)) %>% group_by(dist) %>% summarise(s = sum(contact, na.rm = TRUE)) %>% mutate(dist = dist, prob = s / tot) %>% select(dist, prob)
}


correct_prox_matrix <- function(pm, prox_curve, contact_prob) {
  n1 <- nrow(prox_curve)
  n2 <- nrow(contact_prob)
  if (n1 > n2)
    ratio <- c(contact_prob$prob, rep(NA, n1 - n2)) / prox_curve$score
  else
    ratio <- contact_prob$prob / c(prox_curve$score, rep(NA, n2 - n1))
  ratio[is.infinite(ratio) | is.nan(ratio)] <- NA
  
  disp_matrix <- calc_disp_matrix(length(ratio), length(ratio))
  matrix(ratio[disp_matrix + 1], nrow = length(ratio)) * pm
}