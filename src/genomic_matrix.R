# Under the hood, a genomic matrix is a BEDPE data frame representing data 
# associated with paired genomic intervals, such as HiC data.
#
# Attributes:
# row_coords and col_coords: lists which indicate the genomic coordinates.
# For example: list(chr = "chr22", start = 1000, end = 2000)

library(tidyverse)


.new_genomic_matrix <- function(data, gr, bin_size) {
  stopifnot(length(gr) == 2)
  gr_width <- GenomicRanges::width(gr)
  stopifnot(gr_width[1] == gr_width[2])
  stopifnot(gr_width %% bin_size == 0)
  
  structure(data, # %>% select(start1, start2, score),
            class = c("genomic_matrix", class(data)),
            gr = gr,
            bin_size = bin_size)
}


genomic_matrix <- .new_genomic_matrix


# Serialize a genomic matrix (gm) to BEDPE format (dump order: by column)
# max_records: limit the number of output records
# con: a connection for the serialization
# allow_comments: indicate whether include comment lines in the output
.serialize_bedpe <- function(gm, allow_comments, max_records, to_short = FALSE, conn = NULL) {
  gr <- attr(gm, "gr")
  chr1 <- as.character(GenomicRanges::seqnames(gr)[1])
  chr2 <- as.character(GenomicRanges::seqnames(gr)[2])
  bin_size <- attr(gm, "bin_size")
  
  gm_size <- GenomicRanges::width(gr)[1] %/% bin_size
  
  repr <- list()
  
  if (allow_comments) {
    gr_start1 <- GenomicRanges::start(gr)[1] - 1
    gr_start2 <- GenomicRanges::start(gr)[2] - 1
    gr_end1 <- GenomicRanges::end(gr)[1]
    gr_end2 <- GenomicRanges::end(gr)[2]
    repr <- append(repr, c(
      str_interp("# A genomic matrix: $[d]{gm_size} x $[d]{gm_size}"),
      str_interp("# Genomic range: ${chr1}:$[d]{gr_start1 + 1}-$[d]{gr_end1} vs. ${chr2}:$[d]{gr_start2 + 1}-$[d]{gr_end2}")
    ))
  }
  
  max_records <- ifelse(max_records <= nrow(gm), max_records, nrow(gm))
  
  contents <- 1:max_records %>% map_chr(function(index) {
    start1 <- gm[index,]$start1
    start2 <- gm[index,]$start2
    score <- gm[index,]$score
    
    if (to_short)
      str_interp("0\t${chr1}\t$[d]{start1}\t0\t0\t${chr2}\t$[d]{start2}\t1\t${score}")
    else {
      line_out <- 
        str_interp("${chr1}\t$[d]{start1}\t$[d]{start1+bin_size}\t${chr2}\t$[d]{start2}\t$[d]{start2+bin_size}\t${score}")
      # If there are additional columns
      if (ncol(gm) > 3) {
        line_out <- paste(
          line_out,
          "\t",
          colnames(gm) %>% tail(n = -3) %>%
            map_chr(function(name) {
              as.character(gm[index, ][[name]])
            }) %>%
            paste0(collapse = "\t")
        )
      }
      line_out
    }
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
  gm1 <- gm %>% select(start1, start2, score) %>% filter(is.numeric(score) & is.finite(score))
  writeLines(.serialize_bedpe(gm1, FALSE, nrow(gm1), to_short = TRUE), con = conn)
}


# Output the full genomic matrix to connection
dump_genomic_matrix <- function(gm, conn = stdout()) {
  writeLines(.serialize_bedpe(gm, FALSE, nrow(gm)), con = conn)
}


load_genomic_matrix <- function(conn = stdin(), additional_col_names = NULL, gr = NULL, bin_size = NULL) {
  full_col_names <- c(
    c("chr1", "start1", "end1", "chr2", "start2", "end2", "score"),
    if (is.null(additional_col_names))
      character()
    else
      additional_col_names)
  full_col_types <- c(
    list(
      col_factor(),
      col_integer(),
      col_integer(),
      col_factor(),
      col_integer(),
      col_integer(),
      col_guess()
    ),
    if (is.null(additional_col_names))
      list()
    else
      1:length(additional_col_names) %>% map(function(v)
        col_guess())
  )
  gm_bedpe <-
    read_delim(
      conn,
      delim = "\t",
      col_names = full_col_names,
      col_types = do.call(cols, full_col_types)
    ) %>% mutate(score = as.numeric(score))
  
  # Determine chr
  stopifnot(length(levels(gm_bedpe$chr1)) == 1)
  stopifnot(length(levels(gm_bedpe$chr2)) == 1)
  chr1 <- as.character(gm_bedpe$chr1[1])
  chr2 <- as.character(gm_bedpe$chr2[1])
  
  if (is.null(gr)) {
    gr_start1 <- min(gm_bedpe$start1)
    gr_end1 <- max(gm_bedpe$end1)
    gr_start2 <- min(gm_bedpe$start2)
    gr_end2 <- max(gm_bedpe$end2)
    
    gr <- GenomicRanges::GRanges(c(chr1, chr2),
                                 c(
                                   IRanges::IRanges(gr_start1 + 1, gr_end1),
                                   IRanges::IRanges(gr_start2 + 1, gr_end2)
                                 ))
  }
  
  if (is.null(bin_size))
    bin_size <- with(gm_bedpe[1,], end1 - start1)
  
  genomic_matrix(gm_bedpe %>% select(-chr1, -end1, -chr2, -end2), gr, bin_size)
}

convert_to_matrix <- function(gm) {
  bin_size <- attr(gm, "bin_size")
  gr <- attr(gm, "gr")
  gr_start1 <- GenomicRanges::start(gr)[1] - 1
  gr_start2 <- GenomicRanges::start(gr)[2] - 1
  
  gm_size <- GenomicRanges::width(gr)[1] %/% bin_size
  m <- matrix(NA, nrow = gm_size, ncol = gm_size)
  
  data <- as_tibble(gm) %>% 
    mutate(idx = (start2 - gr_start2) %/% bin_size * gm_size + (start1 - gr_start1) %/% bin_size + 1)
  
  data2 <- as_tibble(gm) %>% 
    mutate(s1 = start1 + start2,
           start1 = s1 - start1,
           start2 = s1 - start2,
           idx = (start2 - gr_start2) %/% bin_size * gm_size + (start1 - gr_start1) %/% bin_size + 1)
  
  m[data$idx] <- data$score
  m[data2$idx] <- data2$score
  m
}

trim_genomic_matrix <- function(gm, gr) {
  bin_size <- attr(gm, "bin_size")
  gr1 <- gr[1]
  gr2 <- if (length(gr) >= 2) gr[2] else gr1
  stopifnot(all(list(gr1, gr2) %>% map_lgl(function(v) {
    GenomicRanges::width(v) %% bin_size == 0
    })))
  stopifnot(GenomicRanges::width(gr1) == GenomicRanges::width(gr2))
  
  limits <- list(gr1, gr2) %>% map(function(v) c(GenomicRanges::start(v) - 1, GenomicRanges::end(v)))
  
  gm %>% as_tibble() %>% 
    filter(between(start1, limits[[1]][1], limits[[1]][2]) & between(start2, limits[[2]][1], limits[[2]][2])) %>%
    genomic_matrix(gr = GenomicRanges::GRanges(list(gr1, gr2) %>% map_chr(as.character)), bin_size = bin_size)
}
