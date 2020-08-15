#!/usr/bin/env Rscript

# Determine the project home directory through the package "here"
# This results in consistent source/import in the future
get_script_path <- function() {
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl = TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if (length(script.dir) == 0)
    stop("can't determine script dir: please call the script with Rscript")
  if (length(script.dir) > 1)
    stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

original_wd <- getwd()
setwd(get_script_path())
requireNamespace("here")
setwd(original_wd)

library(logging)
library(tidyverse)
requireNamespace("optparse")

source(here::here("src/customized_logging.R"))
source(here::here("src/compartment.R"))


logReset()
writeToConn <- writeToConnFactory(stderr())
addHandler(writeToConn)

# Parse common options
parser <- optparse::OptionParser(option_list = list(
  optparse::make_option(
    c("-v", "--verbose"),
    help = "Show verbose messagegs",
    type = "logical",
    default = FALSE,
    action = "store_true"
  )
))
args <-
  optparse::parse_args(parser,
             args = commandArgs(trailingOnly = TRUE),
             positional_arguments = TRUE)
if (args$options$verbose)
  setLevel("DEBUG")

args <- args$args
subcommand <- args[1]

# Parse a numeric string, e.g. 50k, 500K, 10M, 10m, 25.5k, etc.
parse_num_mk <- function(s) {
  s <- str_replace_all(str_trim(s), ",", "")
  
  num_pattern <- "([[:digit:]\\.]+)"
  postfix <- list(k = 1e3L, m = 1e6L)
  patterns <- c(
    sprintf("^%s$", num_pattern),
    sprintf("^%s(k|K)$", num_pattern),
    sprintf("^%s(m|M)$", num_pattern),
    sprintf("^%s(g|G)$", num_pattern))
  
  for (pattern in patterns) {
    match <- str_match(s, pattern)
    if (is.na(match[1])) next
    
    num <- as.numeric(match[2])
    base <- if (length(match) == 3) postfix[[tolower(match[3])]] else 1L
    num <- as.integer(num * base)
    if (is.integer(num)) return(num)
  }
  stop(str_interp("Invalid argument: ${s}"))
}


# Parse the command line arguments
parse_contact_args <- function(args) {
  parser <- optparse::OptionParser(
    option_list = list(
      # optparse::make_option(c("--range"), help = "Genomic range of the region under study"),
      optparse::make_option(
        c("-s", "--bin-size"),
        help = "Size of each bin in base pairs",
        default = "500k"
      ),
      optparse::make_option(
        c("-b", "--block-size"),
        help = "Size of each block in base pairs",
        default = "5m"
      ),
      # optparse::make_option(c("-o", "--output-file"), help = "File name for storing the calculated distance matrix", type = "character"),
      optparse::make_option(c("-m", "--metrics"),
                  help = "The statistical distance metrics. The default is 'ks' for Kolmogorov–Smirnov test.",
                  default = "ks"),
      optparse::make_option(c("-q", "--min-mapq"),
                  help = "Exclude all fragments with smaller MAPQ scores. Default is 0",
                  type = "integer",
                  default = 0),
      # optparse::make_option(c("--min-samples"),
      #             help = "The minimal number of samples from each bin required to perform the calculation.",
      #             type = "integer"),
      optparse::make_option(c("--max-frag-size"),
                  help = "Exclude fragments with larger sizes. Default is 1000",
                  default = 1000,
                  type = "integer"),
      optparse::make_option(
        c("-n", "--num-cores"),
        help = "The number of cores to be used in the computing",
        default = 1
      )
    )
  )
  
  args <- optparse::parse_args(parser, args = args, positional_arguments = TRUE)
  options <- args$options
  
  # Supported distance metrics: ks and cucconi
  stopifnot(options$metrics == "ks" || options$metrics == "cucconi")
  
  # Parse the numerical arguments
  options$`bin-size` <- parse_num_mk(options$`bin-size`)
  options$`block-size` <- parse_num_mk(options$`block-size`)
  
  # Parse the genomic range (example: chr22:100k-150k)
  if (!is.null(options$range)) {
    range_str <- str_trim(options$range)
    match <- str_match(range_str, "^([^:]+):([^-]+)-([^-]+)")
    
    start <- .parse_num_mk(match[3])
    end <- .parse_num_mk(match[4])
    chr <- match[2]
    
    if (is.na(match[1]) || !is.integer(start) || !is.integer(end))
      stop(str_interp("Invalid genomic range specification: ${range_str}"))
    
    options$`range` <- str_interp("${chr}:${start}-${end}")
  }
  
  # if (length(args$args) != 1)
  #   stop("BAM file not specified")
  # options$bam_file <- args$args[1]
  
  return(options)
}



parse_compartment_args <- function(args) {
  parser <- optparse::OptionParser(
    option_list = list(
      optparse::make_option(c("-g", "--genes"), help = "A BED file of gene annotations. Can read from .gz and .bzip2 compresed files.")
    )
  )
  
  optparse::parse_args(parser, args = args)
}




switch (subcommand,
        
        # "distance" = {
        #   logdebug("Calculating distance matrix...")
        #   source(here::here("src/calc_distance.R"))
        #   calc_distance_cli(parse_contact_args(tail(args, n = -1)))
        # },
        
        "contact" = {
          logdebug("Calculating 3D genome contact matrix...")
          source(here::here("src/contact_matrix.R"))
          call_contact_matrix_cli(parse_contact_args(tail(args, n = -1)))
        },
        
        "compartment" = {
          logdebug("Calling A/B compartments...")
          
          options <- parse_compartment_args(tail(args, n = -1))
          call_compartments_cli(options)
        },
        
        {
          logerror("Invalid subcommand")
        })