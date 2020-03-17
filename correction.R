calc_disp_matrix <- function(nrow, ncol) {
  # Return the displacement matrix, i.e. abs(row_index - col_index)
  row_matrix <- matrix(rep(1:nrow, ncol), nrow = nrow, ncol = ncol)
  col_matrix <-
    matrix(rep(1:ncol, 1, each = nrow), nrow = nrow, ncol = ncol)
  abs(row_matrix - col_matrix)
}

calc_freq_curve <- function(data, nthread = 1) {
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
        tmp <- unlist(data_summary)
        tmp
      })
      r
    }
  
  stopCluster(cl)
  
  results
}