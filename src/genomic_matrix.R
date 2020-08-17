# Under the hood, a genomic matrix is a BEDPE data frame representing data 
# associated with paired genomic intervals, such as HiC data.
#
# Attributes:
# row_coords and col_coords: lists which indicate the genomic coordinates.
# For example: list(chr = "chr22", start = 1000, end = 2000)

library(tidyverse)


.new_genomic_matrix <- function(data, chr, gr_start, gr_end, bin_size) {
  structure(data %>% select(start1, start2, score), 
            class = c("genomic_matrix", class(data)),
            chr = chr, gr_start = gr_start, gr_end = gr_end, bin_size = bin_size)
}


genomic_matrix <- .new_genomic_matrix


# Serialize a genomic matrix (gm) to BEDPE format (dump order: by column)
# max_records: limit the number of output records
# con: a connection for the serialization
# allow_comments: indicate whether include comment lines in the output
.serialize_bedpe <- function(gm, allow_comments, max_records, to_short = FALSE, conn = NULL) {
  gr_start <- attr(gm, "gr_start")
  gr_end <- attr(gm, "gr_end")
  bin_size <- attr(gm, "bin_size")
  chr <- attr(gm, "chr")
  
  gm_size <- (gr_end - gr_start) %/% bin_size
  
  repr <- list()
  
  if (allow_comments) {
    repr <- append(repr, c(
      str_interp("# A genomic matrix: $[d]{gm_size} x $[d]{gm_size}"),
      str_interp("# Genomic range: ${chr}:$[d]{gr_start + 1}-$[d]{gr_end}")
    ))
  }
  
  max_records <- ifelse(max_records <= nrow(gm), max_records, nrow(gm))
  
  contents <- 1:max_records %>% map_chr(function(index) {
    start1 <- gm[index,]$start1
    start2 <- gm[index,]$start2
    score <- gm[index,]$score
    
    if (to_short)
      str_interp("0\t${chr}\t$[d]{start1}\t0\t0\t${chr}\t$[d]{start2}\t1\t${score}")
    else
      str_interp("${chr}\t$[d]{start1}\t$[d]{start1+bin_size}\t${chr}\t$[d]{start2}\t$[d]{start2+bin_size}\t${score}")
  }) 
  contents <- contents %>% paste0(collapse = "\n")
  repr <- append(repr, contents)
  
  if (allow_comments && max_records < nrow(gm)) {
    repr <- append(repr, str_interp("# ... with $[d]{nrow(gm) - max_records} more records"))
  }
  
  paste(repr, collapse = "\n")
}


str.genomic_matrix <- function(gm) {
  cat(.serialize_bedpe(gm, TRUE, 10))
}


print.genomic_matrix <- function(gm) {
  cat(.serialize_bedpe(gm, TRUE, 10))
}


# Convert the genomic matrix to Juicer Tools SHORT format
genomic_matrix_to_short <- function(gm, conn = stdout()) {
  writeLines(.serialize_bedpe(gm, FALSE, nrow(gm), to_short = TRUE), con = conn)
}


# Output the full genomic matrix to connection
dump_genomic_matrix <- function(gm, conn = stdout()) {
  writeLines(.serialize_bedpe(gm, FALSE, nrow(gm)), con = conn)
}


load_genomic_matrix <- function(conn = stdin()) {
  gm_bedpe <-
    read_delim(
      conn,
      delim = "\t",
      col_names = c("chr1", "start1", "end1", "chr2", "start2", "end2", "score"),
      col_types = cols(
        col_factor(),
        col_integer(),
        col_integer(),
        col_factor(),
        col_integer(),
        col_integer(),
        col_double()
      )
    )
  
  chr <- gm_bedpe[1,]$chr1
  bin_size <- with(gm_bedpe[1,], end1 - start1)
  gr_start <- min(gm_bedpe$start1)
  gr_end <- max(gm_bedpe$start1) + bin_size
  
  genomic_matrix(gm_bedpe, chr, gr_start, gr_end, bin_size)
}

convert_to_matrix <- function(gm) {
  gr_start <- attr(gm, "gr_start")
  gr_end <- attr(gm, "gr_end")
  bin_size <- attr(gm, "bin_size")
  
  gm <- as_tibble(gm) %>% transmute(
    bin1 = (start1 - gr_start) %/% bin_size + 1,
    bin2 = (start2 - gr_start) %/% bin_size + 1,
    score = score)
  
  gm_size <- (gr_end - gr_start) %/% bin_size
  m <- matrix(nrow = gm_size, ncol = gm_size)
  
  for (idx in 1:nrow(gm)) {
    bin1 <- gm[idx,]$bin1
    bin2 <- gm[idx,]$bin2
    score <- gm[idx,]$score
    m[bin1, bin2] <- score
    if (bin1 != bin2) m[bin2, bin1] <- score
  }
  
  m
}
