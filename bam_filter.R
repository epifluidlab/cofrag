# Filter the BAM file based on flags, mapping quality, and mappability.
# Usage:
# Rscript bam_filter.R --chr 22 --maps 0.75 --mapq 30 --map-track m.bw n.bam

library(logging)
basicConfig()

library(dplyr)
library(doParallel)

bam_filter <-
  function(bam_file,
           map_track,
           output_file,
           grange,
           maps_threshold = 0.75,
           mapq_threshold = 30,
           nthreads = 1) {
    # Only regions with high-mappability scores matter.
    n_map_track_1 <- nrow(map_track)
    map_track <- map_track %>% filter(score >= maps_threshold)
    n_map_track_2 <- nrow(map_track)
    print(sprintf(
      "Filtering the mappability track, from %d to %d obs",
      n_map_track_1,
      n_map_track_2
    ))
    
    convert_map_track <- function(start, end) {
      # Convert the mappability track from BED-like format to mundane coord-series:
      # [(coord, score), ...]
      # Focus on the region between start and end
      series <- rep(0, end - start + 1)
      for (row_idx in 1:nrow(map_track)) {
        row <- map_track[row_idx, ]
        # Scan the map_track and set values accordingly
        if (row$end < start)
          next
        else if (row$start <= start && row$end >= start) {
          idx1 <- 1
          idx2 <- row$end - start + 1
          series[idx1:idx2] <- row$score
        }
        else if (row$start > start && row$start <= end) {
          idx1 <- row$start - start + 1
          if (row$end <= end)
            idx2 <- row$end - start + 1
          else
            idx2 <- end - start + 1
          series[idx1:idx2] <- row$score
        } else
          # row$start > end
          break
      }
      print(sprintf("Qualified mappability coverage: %f", mean(series > 0)))
      series
    }
    
    construct_map_track <- function(start_pos, end_pos, nthreads) {
      # Construct the mappability track in "flat" format.
      # Originally, the track is in BED-like format. For the downstream analysis
      # purpose, we need to convert the track to the following format:
      # c(score1, score2, ...), which is simply a numeric vector reflecting the
      # mappability scores of the region start_pos:end_pos
      
      # Split the task for parallel computing
      stop_points <-
        round(seq(start_pos, end_pos, by = (end_pos - start_pos) / nthreads))
      
      cl <- makeCluster(nthreads)
      registerDoParallel(cl)
      map_track_roi <-
        foreach (
          chunk_idx = 1:(length(stop_points) - 1),
          .combine = c,
          .export = "convert_map_track"
        ) %dopar% {
          n1 = stop_points[chunk_idx]
          if (chunk_idx == length(stop_points) - 1)
            n2 = stop_points[chunk_idx + 1]
          else
            n2 = stop_points[chunk_idx + 1] - 1
          convert_map_track(n1, n2)
        }
      
      stopCluster(cl)
      map_track_roi
    }
    
    # Construct the BAM reads filter
    filter_func <- function(bam_records) {
      # Filter #1: filter reads according to the flag and mapping score
      filter_flags <-
        bam_records$flag %in% c(83, 99) &
        bam_records$mapq >= mapq_threshold
      filtered_indexes <- which(filter_flags)
      
      print(sprintf(
        "Filtering BAM chunks between coords: %d-%d.",
        bam_records$pos[1],
        tail(bam_records$pos, n = 1)
      ))
      print(sprintf("Original total reads: %d", length(bam_records$pos)))
      print(sprintf("Post Filter #1 reads: %d", length(filtered_indexes)))
      
      # Among reads that passes Filter #1, apply Filter #2 based on the mappability score
      
      # Extract the region-of-interest from the mappability track
      start_pos <- bam_records$pos[1]
      end_pos <- tail(bam_records$pos, n = 1)
      map_track_roi <- construct_map_track(start_pos, end_pos, nthreads)
      
      filter_flags <-
        filter_flags &
        sapply(bam_records$pos, function(x)
          map_track_roi[x - start_pos + 1] >= maps_threshold)
      
      print(sprintf("Post Filter #2 reads: %d", sum(filter_flags)))
      filter_flags
    }
    
    which <- GRanges(grange)
    param <- ScanBamParam(which = which,
                          what = c("flag", "pos", "seq", "isize", "qwidth", "mapq"))
    
    print(sprintf("Start filtering BAM: %s", grange))
    
    filterBam(
      file = bam_file,
      destination = output_file,
      param = param,
      filter = FilterRules(list(Func = filter_func))
    )
  }

main <- function(arguments) {
  # Load the mappability track
  library(rtracklayer)
  library(Repitools)
  library(Rsamtools)
  library(stringr)
  
  match_result <- str_match(arguments$grange, "(chr)?([^.]+):(.+)")
  chr <- match_result[1, 3]
  chr_name <- sprintf("chr%s", chr)
  arguments$grange <- sprintf("%s:%s", chr, match_result[1, 4])
  print(sprintf("Loading the mappability track for %s...", chr_name))
  
  map_track <-
    annoGR2DF(import(arguments$map_track)) %>% filter(chr == chr_name) %>% select(start, end, score)
  print("Completed loading the mappability track.")
  
  print(sprintf("Start filtering the BAM file: %s", arguments$bam_file))
  print(sprintf("Output file: %s", arguments$output_file))
  print(
    sprintf(
      "With parameters: maps_threshold: %f, mapq_threshold: %f",
      arguments$maps,
      arguments$mapq
    )
  )
  print(sprintf("Mappability track: %s", arguments$map_track))
  
  bam_filter(
    bam_file =arguments$bam_file,
    map_track = map_track,
    output_file = arguments$output_file,
    grange = arguments$grange,
    maps_threshold = arguments$maps,
    mapq_threshold = arguments$mapq,
    nthreads = arguments$nthreads
  )
}

process_arguments <- function() {
  # Default value
  arguments <- list(maps = 0.75, mapq = 30, nthreads = 1)
  args_list <- commandArgs(trailingOnly = TRUE)
  idx <- 1
  while (idx <= length(args_list)) {
    action <- args_list[idx]
    if (action == "--grange") {
      idx <- idx + 1
      arguments$grange <- args_list[idx]
    } else if (action == "--maps") {
      idx <- idx + 1
      arguments$maps <- as.numeric(args_list[idx])
    } else if (action == "--mapq") {
      idx <- idx + 1
      arguments$mapq <- as.numeric(args_list[idx])
    } else if (action == "--map-track") {
      idx <- idx + 1
      arguments$map_track <- args_list[idx]
    } else if (action == "-o") {
      idx <- idx + 1
      arguments$output_file <- args_list[idx]
    } else if (action == "-n") {
      idx <- idx + 1
      arguments$nthreads <- as.integer(args_list[idx])
    } else if (idx == length(args_list)) {
      arguments$bam_file <- args_list[idx]
      break
    }
    idx <- idx + 1
  }
  arguments
}

main(process_arguments())
# main(
#   list(
#     map_track = "./hg19_mappability/hg19_45_mer.chr22.bw",
#     # map_track = "./hg19_mappability/hl_chr21_45.bw",
#     bam_file = "./SRR2129993.chr22.part1.sorted.bam",
#     output_file = "./filter.bam",
#     grange = "chr22:1-51304566",
#     maps = 0.75,
#     mapq = 30,
#     nthreads = 1
#   )
# )
