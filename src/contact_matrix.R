# This is the main module for the contact matrix analysis, based on fragmentation data

library(foreach)
library(doParallel)
library(logging)
library(tidyverse)
requireNamespace("GenomicRanges")


source(here::here("src/genomic_matrix.R"))
source(here::here("src/customized_logging.R"))


# Read fragments data from connection
# Return a
import_fragments <- function(conn) {
  frag_data <-
    read_delim(
      conn,
      delim = "\t",
      col_names = c("chr", "start", "end", "mapq", "strand"),
      col_types = cols(
        col_factor(),
        col_integer(),
        col_integer(),
        col_integer(),
        col_character()
      )
    )
  
  # Only work for intra-chromosome cases
  chr_names <- unique(frag_data$chr)
  stopifnot(length(chr_names) == 1)
  chr <- as.character(unique(frag_data$chr))[1]
  
  loginfo("Finished importing fragments")
  loginfo(str_interp("Total fragments #: ${nrow(frag_data)}"))
  
  range_start <- head(frag_data, n = 1)$start
  range_end <- tail(frag_data, n = 1)$end
  loginfo(str_interp("Range: ${chr}:${range_start + 1}-${range_end}"))
  
  frag_data %>% mutate(length = end - start, mid = as.integer((end + start) / 2))
}


# In order to calculate the distance matrix, we first need to infer the genomic
# range from the fragments, and make paddings as necessary so that: the entire
# genomic range is multiple of bin_size
infer_genomic_range <- function(frag_data, bin_size) {
  # Infer the genomic range
  chr_name <- frag_data[1,]$chr
  range_start <- min(frag_data$start)
  range_end <- max(frag_data$end)
  
  # Start is round to multiples of bin_size. Notice that range_start is 0-based and semi-inclusive
  # For example: padded range: 1000-2500 (BED coordinates, bin_size=500)
  logdebug(str_interp("range_start: $[d]{range_start}"))
  padded_start <- range_start %/% bin_size * bin_size
  padded_end <-
    (range_end - padded_start) %/% bin_size * bin_size + padded_start
  logdebug(
    str_interp(
      "padded_start: $[d]{padded_start}, bin_size: $[d]{bin_size}, padded_end: $[d]{padded_end}"
    )
  )
  if (padded_end < range_end)
    padded_end <- padded_end + bin_size
  
  GenomicRanges::GRanges(str_interp("${chr_name}:$[d]{padded_start + 1}-$[d]{padded_end}"))
}


# Calculate binned fragmentation profile
# Parameters:
#   * frag_data: BED-format data frame representing fragments
#   * gr: a GRanges object indicating the region of interest
#   * bin_size
#
# Return a data frame containing fragment lengths for each bin:
# # A tibble: 704 x 3
#     bin_idx bin_start frag
#       <int>     <int> <list>
#   1       1  16100000 <int [18,163]>
#   2       2  16150000 <int [39,340]>
calc_bfp <-
  function(frag_data, gr, bin_size) {
    stopifnot(GenomicRanges::width(gr) %% bin_size == 0)
    
    # Build the bin layout
    gr_start <- GenomicRanges::start(gr) - 1
    bin_layout <-
      sapply(1:(GenomicRanges::width(gr) %/% bin_size), function(v) {
        as.integer(gr_start + (v - 1) * bin_size)
      })
    
    # Calculate the bin index for a certain coordinate (pos)
    .calc_bin_index <- function(pos, gr_start, bin_size) {
      as.integer((pos - gr_start) %/% bin_size + 1)
    }
    
    # bin1: index of the bin which contains the fragment start
    # bin2: index of the bin which contains the fragment end
    # cross_bin: logical value indicating whether the fragment straddles across adjacent bins
    frag_data <- frag_data %>% mutate(
      bin = as.integer(.calc_bin_index(mid, gr_start, bin_size))
      # bin1 = as.integer(.calc_bin_index(start, gr_start, bin_size)),
      # bin2 = as.integer(.calc_bin_index(end - 1, gr_start, bin_size)),
      # cross_bin = (bin1 != bin2)
    )
    
    frag_coll <- 
      frag_data %>% group_by(bin) %>% summarize(frag_coll = list(list(length, start)))
    
    # frag_coll_1 <-
    #   frag_data %>% group_by(bin1) %>% summarize(frag_coll = list(list(length, start)))
    # # For cross-bin fragments, group them by bins. When retrieving fragments in
    # # a certain bin, we should not only consider the fragments starting in the
    # # bin, but also those ending in the bin. Thus, we need to combine
    # # frag_coll_1 and frag_coll_2
    # frag_coll_2 <-
    #   frag_data %>% filter(cross_bin == TRUE) %>% group_by(bin2) %>% summarize(frag_coll = list(list(length, start)))
    
    # A helper function which retrieve all fragments for a specified bin
    # Properly deal with cases where the bin contains no data
    retrieve_frag_coll <- function(data, bin_name, bin_idx) {
      results <- data %>% filter(.data[[bin_name]] == bin_idx)
      if (nrow(results) > 0)
        results$frag_coll[[1]]
      else
        list(integer(0), integer(0))
    }
    
    as_tibble(list(
      bin_idx = seq_along(bin_layout),
      bin_start = bin_layout,
      frag = lapply(seq_along(bin_layout), function(v) {
        retrieve_frag_coll(frag_coll, "bin", v)
        # coll1 <- retrieve_frag_coll(frag_coll_1, "bin1", v)
        # coll2 <- retrieve_frag_coll(frag_coll_2, "bin2", v)
        # # Concatenate these two by row-binding
        # 
        # list(c(coll1[[1]], coll2[[1]]), c(coll1[[2]], coll2[[2]]))
      })
    ))
  }



.build_block_layout <- function(gr, block_size) {
  gr_start <- GenomicRanges::start(gr) - 1
  gr_end <- GenomicRanges::end(gr)
  block_start <-
    as.integer(seq(gr_start, gr_end, by = block_size) %>% keep(function(x)
      x < gr_end))
  block_end <- (block_start + block_size) %>% map_int(function(x) {
    as.integer(ifelse(x < gr_end, x, gr_end))
  })
  as_tibble(list(start = block_start, end = block_end))
}


# Divide the genomic range into blocks, and return a list of paired blocks from the upper-right half of the matrix
# Example:
#       i     j start1  end1 start2  end2
#   <int> <int>  <int> <int>  <int> <int>
# 1     1     1      0   300      0   300
# 2     1     2      0   300    300   600
# 3     2     2    300   600    300   600
.build_block_pairs <- function(gr, block_size) {
  gr_start <- GenomicRanges::start(gr) - 1
  gr_end <- GenomicRanges::end(gr)
  num_blocks <-
    length(seq(gr_start, gr_end, by = block_size) %>% keep(function(x)
      x < gr_end))
  
  block_layout <- .build_block_layout(gr, block_size)
  
  as_tibble(expand.grid(1:num_blocks, 1:num_blocks) %>% filter(Var2 >= Var1)) %>%
    transmute(
      i = Var1,
      j = Var2,
      start1 = block_layout$start[i],
      end1 = block_layout$end[i],
      start2 = block_layout$start[j],
      end2 = block_layout$end[j]
    )
}


consolidate_models <- function(model_results) {
  # For fraglen models: flip
  m <- model_results[["fraglen"]]
  
  max_score <- max(m$score)
  m %>% mutate(score = max_score - score)
}


# IMPORTANT: define a named list of models and call contact matrix based on each of them
init_models <- function(...) {
  fraglen_env <- new.env()
  source(here::here("src/fraglen_model.R"), local = fraglen_env)
  
  # fraglen_model <- function(bfp, gr, bin_size, interval1, interval2, metrics) {
  list(fraglen = function(metrics = "ks", ...) {
    fraglen_env$fraglen_model(metrics = "ks", ...)
  })
}


# Each model is a function(bfp, gr, bin_size, interval1, interval2, ...) {
call_contact_matrix <-
  function(frag_data,
           model_envs,
           genome_tracks = NULL,
           bin_size = 500000,
           block_size = 5000000,
           ncores = 1,
           rng_seed = NULL,
           ...) {
    # Calculate the distance matrix
    # Parameters:
    # * frag_data: a BED-format data frame containing fragments
    # * bin_zize
    # * block_size: the calculation is performed based on "blocks". Notice that
    #     block_size should be multiples of bin_size
    # * ncores: multi-processing
    # * opts: other options for the calculation
    
    # Argument check
    # block_size should be multiple of bin_size
    if (block_size %% bin_size != 0)
      stop(str_interp(
        "Block size ${block_size} should be multiple of bin size ${bin_size}"
      ))
    
    gr <- infer_genomic_range(frag_data, bin_size)
    
    # Calculate the binned fragment profile (BFP)
    bfp <- calc_bfp(frag_data, gr, bin_size)
    
    # Divide the genomic range into blocks
    block_pairs <- .build_block_pairs(gr, block_size)
    
    interval_pairs <- 1:nrow(block_pairs) %>% map(function(idx) {
      pair <- block_pairs[idx, ]
      interval1 <-
        GenomicRanges::GRanges(
          GenomicRanges::seqnames(gr),
          str_interp("${pair[['start1']] + 1}-${pair[['end1']]}")
        )
      interval2 <-
        GenomicRanges::GRanges(
          GenomicRanges::seqnames(gr),
          str_interp("${pair[['start2']] + 1}-${pair[['end2']]}")
        )
      list(interval1, interval2)
    })
    
    args <- list(...)
    
    models <- model_envs %>% map(function(model_env) {
      # Clone and create the working model environment
      parent.env(model_env) <- environment()
      working_env <- new.env(parent = model_env)
      with(working_env, {
        model_args <- c(list(bfp = bfp, gr = gr, bin_size = bin_size), args)
        entry <- function(interval1, interval2) {
          model_args <- c(list(interval1 = interval1, interval2 = interval2), model_args)
          do.call(model_func, model_args)
        }
      })
      working_env
    })
    names(models) <- names(model_envs)
    # 
    # # Init models
    # fraglen_env <- new.env()
    # source(here::here("src/fraglen_model.R"), local = fraglen_env)
    # with(fraglen_env, {
    #   metrics <- ifelse(is.null(args$metrics), "ks", args$metrics)
    #   entryfunc <- function(interval1, interval2) {
    #     fraglen_model(bfp, gr, bin_size, interval1, interval2, metrics)
    #   }
    # })
    # models <- list(fraglen = fraglen_env)

    
    raw_matrices <- lapply(names(models), function(model_name) {
      loginfo(str_interp("Calling contact matrix baed on model: ${model_name}"))
      model <- models[[model_name]]
      
      loginfo(str_interp("Registering cluster for ${ncores} cores"))
      logger_level <- getLogger()[['level']]
      
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      results <-
        foreach (
          pair_index = 1:length(interval_pairs),
          .combine = "rbind",
          .export = c("interval_pairs", "rng_seed")
        ) %dopar% {
          library(logging)
          library(tidyverse)
          
          if (!is.null(rng_seed)) {
            set.seed(rng_seed + pair_index)
          }

          # Create a logger for the worker
          logger <- getLogger('worker')
          addHandler(writeToFile, logger = "worker", file = "/dev/stderr")
          setLevel(logger_level, container = "worker")

          pair <- interval_pairs[[pair_index]]
          interval1 <- pair[[1]]
          interval2 <- pair[[2]]
          loginfo(
            str_interp(
              "Processing block #${pair_index} / ${length(interval_pairs)}: ${interval1} vs. ${interval2}"
            ),
            logger = "worker"
          )
          
          model$entry(interval1, interval2)
        }
      stopCluster(cl)
      
      chr_name <- as.character(GenomicRanges::seqnames(gr))
      gr_start <- GenomicRanges::start(gr) - 1
      gr_end <- GenomicRanges::end(gr)
      
      gm_gr <- GenomicRanges::GRanges(c(chr_name, chr_name),
                                      c(IRanges::IRanges(gr_start + 1, gr_end),
                                        IRanges::IRanges(gr_start + 1, gr_end)))
      
      genomic_matrix(results, gm_gr, bin_size)
    })
    names(raw_matrices) <- names(models)
    
    consolidate_models(raw_matrices)
  }


# Entrypoint for distance calculation
call_contact_matrix_cli <- function(options) {
  loginfo("Command line arguments:")
  loginfo(paste(
    names(options),
    options,
    sep = "=",
    collapse = "; "
  ))
  
  # Import fragmentation data from stdin
  frag_data <-
    import_fragments(conn = fifo("/dev/stdin")) %>% filter(mapq >= options$`min-mapq` &
                                                              length <= options$`max-frag-size`)
  
  model_envs <- list(fraglen = new.env())
  source(here::here("src/fraglen_model.R"), local = model_envs$fraglen)
  model_envs$fraglen$model_func <- model_envs$fraglen$fraglen_model
  
  contact_matrix <- call_contact_matrix(
    frag_data,
    model_envs,
    genome_tracks = NULL,
    bin_size = options$`bin-size`,
    block_size = options$`block-size`,
    ncores = options$`num-cores`,
    metrics = options$metrics,
    rng_seed = options$seed,
    bootstrap = options$bootstrap,
    subsample = options$subsample
  )
  
  dump_genomic_matrix(contact_matrix)
}