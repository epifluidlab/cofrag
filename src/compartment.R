source(here::here("src/genomic_matrix.R"))


# Calculate the number of genes in a genomic interval
.calc_num_genes <- function(gr, gene_density) {
  helper <- function(interv_chr, interv_start, interv_end) {
    gene_density %>%
      filter(chr == as.character(interv_chr) &
               end > interv_start & start < interv_end) %>% nrow()
  }
  
  sapply(1:length(gr), function(idx) {
    interval <- gr[idx]
    helper(
      as.character(GenomicRanges::seqnames(interval)),
      as.integer(GenomicRanges::start(interval) - 1),
      as.integer(GenomicRanges::end(interval))
    )
  })
}


# Perform the A/B compartment analysis
# Parameters:
#   * gm: a genomic_matrix object
#   * gene_density: a BED-format data frame of gene annotation
#   * pc_compartment: indicates which component is related to compartment
#   * n_component: return the first n_component components
# Return: a vector indicating the compartment level
call_compartments <- function(gm, 
                              gene_density = NULL,
                              valid_strips_pct = 0.9,
                              center = TRUE,
                              scale = FALSE,
                              pc_compartment = 1,
                              n_component = 1,
                              rollmean_width = NULL) {
  m <- convert_to_matrix(gm)
  # invalid_values <- is.infinite(m) | is.na(m)
  # if (sum(invalid_values) > 0 & !is.null(nuance))
  #   m[invalid_values] <- rnorm(sum(invalid_values), sd = nuance)
  
  # Identify all-NA strips and remove them
  # all_na <- 1:ncol(m) %>% map_lgl(~ m[, .] %>% is.na() %>% all())
  
  # Identify strips of which more than 20% are missing
  invalid_strips <- 1:ncol(m) %>% map_lgl(~ (m[, .] %>% is.na() %>% mean) >= valid_strips_pct)
  valid_m <- m[!invalid_strips, !invalid_strips]
  
  eigenvectors <-
    prcomp(
      cor(valid_m, use = "pairwise.complete.obs", method = "pearson"), 
      center = TRUE,
      scale. = scale,
      )$rotation  #[, 1]
  
  comp <- rep(NA, ncol(m))
  comp[!invalid_strips] <- eigenvectors[, pc_compartment]
  
  bin_size <- attr(gm, "bin_size")
  gr <- attr(gm, "gr")
  stopifnot(as.character(gr[1]) == as.character(gr[2]))
  chr <- as.character(GenomicRanges::seqnames(gr))[1]
  gr_start <- GenomicRanges::start(gr)[1] - 1
  gr_end <- GenomicRanges::end(gr)[1]

  gr <- 1:((gr_end - gr_start) %/% bin_size) %>%
    map_chr(function(idx) {
      start <- gr_start + (idx - 1) * bin_size
      str_interp("${chr}:$[d]{start + 1}-$[d]{start + bin_size}")
    }) %>%
    GenomicRanges::GRanges()
  
  # Compare the compartment level with gene density and determine whether to flip
  if (!is.null(gene_density)) {
    gd_cnt <- .calc_num_genes(gr, gene_density)
    if (cor(comp, gd_cnt, use = "complete.obs") < 0)
      comp <- -comp
  }
  
  result <- list(
    start = GenomicRanges::start(gr) - 1, 
    end = GenomicRanges::end(gr), 
    score = comp) %>%
    as_tibble()
  
  if (center) {
    result <- result %>% mutate(score = score - mean(score, na.rm = TRUE))
  }
  
  if (!is.null(rollmean_width)) {
    result <- result %>% 
      mutate(score = zoo::rollmean(score, rollmean_width, fill = NA))
  }
  
  # Extract the first n_component eigenvectors
  pc <- matrix(rep(NA, ncol(m) * n_component), ncol = n_component)
  pc[!invalid_strips,] <- eigenvectors[, 1:n_component]
  pc %<>% data.frame() %>% as_tibble() %>% set_colnames(1:n_component %>% map_chr(~ paste0("PC", .)))
  
  return(bind_cols(result, pc))
}


load_gene_density <- function(conn) {
  read_tsv(
    conn,
    col_names = c("chr", "start", "end", "gene"),
    col_types = list(col_factor(), col_integer(), col_integer(), col_factor())
  )
}


call_compartments_cli <- function(options) {
  gene_annot <- options$genes
  logdebug(str_interp("Loading gene annotation: ${gene_annot}"))
  
  conn <- if (endsWith(gene_annot, ".gz")) {
    gzfile(gene_annot)
  } else if (endsWith(gene_annot, ".bzip2")) {
    bzfile(gene_annot)
  } else {
    gene_annot
    # file(gene_annot)
  }
  
  logdebug(str_interp("Loading gene density file from ${gene_annot}..."))
  gene_density <- load_gene_density(conn)
  
  logdebug("Loading genomic matrix...")
  gm <- load_genomic_matrix(conn = fifo("/dev/stdin"))
  
  logdebug("Calling compartments...")
  comp <- call_compartments(gm = gm, gene_density = gene_density)
  
  chr <- attr(gm, "chr")
  
  contents <- 1:nrow(comp) %>% 
    map_chr(function(idx) {
      line <- comp[idx,]
      start <- line$start + 1
      end <- line$end
      score <- line$score
      str_interp("${chr}\t$[d]{start}\t$[d]{end}\t${score}")
    })
  writeLines(contents)
}


calc_pearson_correlation <- function(comp1, comp2, use = "everything", method = "pearson") {
  # Check validity
  bin_size1 <- (comp1$end - comp1$start) %>% unique()
  bin_size2 <- (comp2$end - comp2$start) %>% unique()
  stopifnot(length(bin_size1) == 1 && length(bin_size2) == 1 && bin_size1[1] == bin_size2[1])
  bin_size = bin_size1[1]
  
  # Align the two compartment dataset
  comp1 <- comp1 %>% mutate(bin_idx = start %/% bin_size) %>% 
    filter(!is.na(score))
  comp2 <- comp2 %>% mutate(bin_idx = start %/% bin_size) %>%
    filter(!is.na(score))
  
  common_bins <- intersect(comp1$bin_idx, comp2$bin_idx)
  comp1 <- comp1 %>% filter(bin_idx %in% common_bins)
  comp2 <- comp2 %>% filter(bin_idx %in% common_bins)
  
  cor(comp1$score, comp2$score, use = use, method = method)
}


