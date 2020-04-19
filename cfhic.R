#!/usr/bin/env Rscript

process_args <- function() {
  library(optparse)
  
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0)
    stop("No subcommand found")
  
  subcommand <- args[[1]]
  args <- tail(args, length(args) - 1)
  
  if (subcommand == "distance") {
    list(subcommand = subcommand, args = process_dist_args(args))
  } else {
    stop(sprintf("Unknown subcommand: %s", subcommand))
  }
}


main <- function() {
  result <- process_args()
  subcommand <- result$subcommand
  args <- result$args
  
  if (subcommand == "distance") {
    source(here::here("calc_distance.R"))
    dm <- calc_distance(
      bam_file = args$bam_file,
      gr = args$range,
      bin_size = args$bin_size,
      block_size = args$block_size,
      nthreads = args$ncores,
      metrics = args$metrics,
      opts = list(min_samples = args$min_samples, max_frag_size = args$max_frag_size)
    )
    if (!is.null(args$output_file)) {
      write.table(
        dm,
        file = args$output_file,
        sep = "\t",
        row.names = T,
        col.names = T
      )
    }
  }
}

process_dist_args <- function(args) {
  parser <- OptionParser(
    option_list = list(
      make_option(c("--range"), help = "Genomic range of the region under study"),
      make_option(c("-s", "--bin-size"), help = "Size of each bin in base pairs", type = "integer"),
      make_option(c("-b", "--block-size"), help = "Size of each block in base pairs", type = "integer"),
      make_option(c("-o", "--output-file"), help = "File name for storing the calculated distance matrix", type = "character"),
      make_option(
        c("-m", "--metrics"), 
        help = "The statistical distance metrics. Currently only 'ks' (Kolmogorov–Smirnov test) is supported.",
        default = "ks"
      ),
      make_option(
        c("--min-samples"),
        help = "The minimul number of samples from each bin required to perform the calculation.",
        default = 1000,
        type = "integer"
      ),
      make_option(
        c("--max-frag-size"),
        help = "Filter out any fragments with larger sizes",
        default = 240,
        type = "integer"
      ),
      make_option(
        c("-n", "--number-of-cores"),
        help = "The number of cores to be used in the computing",
        default = 1,
        dest = "ncores"
      )
    )
  )
  
  args <-
    parse_args(parser, args = args, positional_arguments = TRUE)
  
  if (length(args$args) != 1)
    stop("Please specify one (1) filtered BAM file to process")
  
  options <- args$options
  stopifnot(options$metrics == "ks" || options$metrics == "cucconi")
  list(
    range = options$range,
    bin_size = options$`bin-size`,
    block_size = options$`block-size`,
    output_file = options$`output-file`,
    ncores = options$ncores,
    metrics = options$metrics,
    min_samples = options$`min-samples`,
    max_frag_size = options$`max-frag-size`,
    bam_file = args$args[1]
  )
}

main()