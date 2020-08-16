library(stringr)

# Convert a genomic matrix to BEDPE format
# A genomic matrix is a matrix representing data associated with paired genomic intervals, such as HiC data
# row_coords and col_coords are lists which indicate the genomic coordinates. 
# For example: (chr = "chr22", start = 1000, end = 2000)
# If row_coords and col_coords are NULL, the information will be inferred from rownames/colnames of data
genomic_matrix_to_bed <- function(data, con = stdout(), row_coords = NULL, col_coords = NULL) {
  infer_coords <- function(s) {
    s <- str_trim(s)
    # Example: chr22:2,258,501-2,258,750
    # Notice: 1-based inclusive interval
    match <- str_match(s, "^([^:]+):([[:digit:],]+)-([[:digit:],]+)$")
    
    chr <- match[2]
    start <- as.integer(str_replace_all(match[3], ",", ""))
    end <- as.integer(str_replace_all(match[4], ",", ""))
    if (is.na(match[1]) || !is.integer(start) || !is.integer(end))
      stop(str_interp("Invalid coords: ${s}"))
    
    return(list(chr = chr, start = start, end = end))
  }
  
  if (is.null(row_coords))
    row_coords <- lapply(rownames(data), infer_coords)
  if (is.null(col_coords))
    col_coords <- lapply(colnames(data), infer_coords)
  
  # Dump the matrix by column
  for (col_index in 1:ncol(data)) {
    cc <- col_coords[[col_index]]
    for (row_index in 1:nrow(data)) {
      rc <- row_coords[[row_index]]
      score <- data[row_index, col_index]
      
      line <- paste0(c(cc$chr, cc$start - 1, cc$end, rc$chr, rc$start - 1, rc$end, score), collapse = "\t")
      writeLines(line, con = con)
    }
  }
}


load_distance_matrix <- function(file_name, distance_cap = NULL, col_row_names = FALSE, chr = NULL, range_start = NULL, bin_size = NULL) {
  # col_row_names: if TRUE, colnames and rownames will be calculated from range_start and bin_size
  dm <- as.matrix(read.table(file_name, as.is = TRUE))
  stopifnot(nrow(dm) == ncol(dm))
  
  # The column names of distance matrices loaded from the file will have different form,
  # such as X14.19000001, which actually should be 14:19000001.
  # The row names are always correct.
  colnames(dm) <- rownames(dm)
  
  # Clip entries with extremely large distance
  if (!is.null(distance_cap))
    dm[dm > distance_cap] <- distance_cap
  
  if (col_row_names) {
    colnames(dm) <- paste0(chr, ":", sapply(1:nrow(dm), function(v) (v - 1) * bin_size + 1))
    rownames(dm) <- colnames(dm)
  }
  
  dm
}


# Load Hi-C "obs" data from file
# Data format: Juicer Tools dump command
# Return a genomic_matrix object
load_hic_obs <- function(file_name, chr, bin_size) {
  hic_data <- read_tsv(file_name, col_names = c("start1", "start2", "score"), col_types = list(
    col_integer(),
    col_integer(),
    col_double()
  ))
  
  coords <- c(hic_data$start1, hic_data$start2)
  gr_start <- min(coords)
  gr_end <- max(coords) + bin_size
  
  
  source(here::here("src/genomic_matrix.R"), local = TRUE)
  genomic_matrix(hic_data, chr = chr, gr_start = gr_start, gr_end = gr_end, bin_size = bin_size)
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
  
  fileConn <- if (tools::file_ext(output_file) == "gz") gzfile(output_file) else file(output_file)
  writeLines(output_lines, fileConn)
  close(fileConn)
}
