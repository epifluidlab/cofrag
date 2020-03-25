load_distance_matrix <- function(file_name, distance_cap = NULL) {
  dm <- load_matrix(file_name)
  # Cap entries with extremely large distance
  if (!is.null(distance_cap))
    dm[dm > distance_cap] <- distance_cap
  dm
}

# Load a matrix from file
load_matrix <- function(file_name) {
  dm <- as.matrix(read.table(file_name, as.is = TRUE))
  colnames(dm) <- rownames(dm)
  dm
}

# Convert a distance matrix to a proximity matrix.
# A proximity matrix is simply a linear transformed distance matrix on a 0-1 scale
dist2prox <- function(data) 1 - data / max(data, na.rm = TRUE)


# Convert the matrix to "short format" which is required by Pre for .hic generation
# https://github.com/aidenlab/juicer/wiki/Pre#short-format
matrix_to_short <- function(data, output_file, nthread = 1) {
  # Parameters:
  # * data: a square symmetric matrix. The row names and column names should be an identical
  #     vector, each element specify the corresponding genomic coordinate.
  # * output_file: name of the output file.
  
  # Argument check
  stopifnot(nrow(data) == ncol(data))
  rnames <- rownames(data)
  cnames <- colnames(data)
  stopifnot(identical(rnames, cnames))
  stopifnot(!is.null(rnames))
  stopifnot(!is.null(output_file))
  
  library(stringr)
  
  gcoord_split <- function(s) {
    r <- str_match(s, "([^:]+):(\\d+)")
    chr <- r[2]
    pos <- r[3]
    stopifnot(!is.na(chr))
    stopifnot(!is.na(pos))
    c(chr, pos)
  }
  
  # Construct the tasks assignment
  library(dplyr)
  # Only process upper triangle entries
  tasks <- as.matrix(expand.grid(1:nrow(data), 1:ncol(data)) %>% filter(Var1 <= Var2))
  
  library(foreach)
  library(doParallel)
  
  cl <- makeCluster(nthread)
  registerDoParallel(cl)
  
  tasks_range <- round(seq(from = 1, to = nrow(tasks) + 1, length.out = nthread + 1))
  output_lines <- foreach(task_id = 1:nthread, .combine = "c", .export = "str_match") %dopar% {
    task_start <- tasks_range[task_id]
    task_end <- tasks_range[task_id + 1] - 1
    lines <- sapply(task_start:task_end, function(idx) {
      row_idx <- tasks[idx, 1]
      col_idx <- tasks[idx, 2]
      
      if (is.na(data[row_idx, col_idx]))
        return(NA)
      
      row_gcoord <- gcoord_split(rnames[row_idx])
      col_gcoord <- gcoord_split(cnames[col_idx])
      sprintf("0 %s %s 0 0 %s %s 1 %f",
              row_gcoord[1],
              row_gcoord[2],
              col_gcoord[1],
              col_gcoord[2],
              data[row_idx, col_idx])
    })
    lines[!is.na(lines)]
  }
  stopCluster(cl)
  
  fileConn <- file(output_file)
  writeLines(output_lines, fileConn)
  close(fileConn)
}
