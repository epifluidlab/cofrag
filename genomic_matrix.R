# Under the hood, a genomic matrix is a matrix representing data associated with
# paired genomic intervals, such as HiC data.
#
# Attributes:
# row_coords and col_coords: lists which indicate the genomic coordinates.
# For example: list(chr = "chr22", start = 1000, end = 2000)

library(stringr)

.new_genomic_matrix <-
  function(data, row_intervals, col_intervals) {
    if (nrow(data) != length(row_intervals) ||
        ncol(data) != length(col_intervals))
      stop("Mismatched data and interval specification")
    
    structure(
      data,
      row_intervals = row_intervals,
      col_intervals = col_intervals,
      class = "genomic_matrix"
    )
  }

genomic_matrix <- function(data, row_intervals, col_intervals) {
  .new_genomic_matrix(data, row_intervals = row_intervals, col_intervals = col_intervals)
}

# Serialize a genomic matrix (gm) to BEDPE format (dump order: by column)
# max_records: limit the number of output records
# con: a connection for the serialization
# allow_comments: indicate whether include comment lines in the output
.serialize_bedpe <- function(gm, con, allow_comments, max_records) {
  row_intervals <- attr(gm, "row_intervals")
  col_intervals <- attr(gm, "col_intervals")
  
  cnt <- 0
  
  if (allow_comments)
    cat(str_interp("# A genomic matrix: ${nrow(gm)} x ${ncol(gm)}\n"), file = con)
  
  for (col_index in 1:ncol(gm)) {
    cinterval <- col_intervals[[col_index]]
    for (row_index in 1:nrow(gm)) {
      rinterval <- row_intervals[[row_index]]
      score <- gm[row_index, col_index]
      
      line <- paste(
        cinterval$chr,
        cinterval$start,
        cinterval$end,
        rinterval$chr,
        rinterval$start,
        rinterval$end,
        score,
        sep = "\t"
      )
      cat(line, file = con)
      cat("\n", file = con)
      
      cnt <- cnt + 1
      if (cnt >= max_records) break
    }
    if (cnt >= max_records) break
  }
  if (allow_comments && cnt < length(gm)) {
    cat(str_interp("# ... with ${length(gm) - cnt} more records\n"), file = con)
  }
}

str.genomic_matrix <- function(gm) {
  .serialize_bedpe(gm, stdout(), TRUE, 10)
}

print.genomic_matrix <- function(gm) {
  .serialize_bedpe(gm, stdout(), TRUE, 10)
}

# Output the full genomic matrix to connection
write_genomic_matrix <- function(gm, con = stdout()) {
  .serialize_bedpe(gm, con, FALSE, length(gm))
}
