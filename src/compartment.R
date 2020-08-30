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
#   * input_as_cor: if TRUE, the input gm is already a correlation matrix, so that we can 
#                   skip calling cor()
# Return: a vector indicating the compartment level
call_compartments <- function(gm, 
                              input_as_cor = FALSE, 
                              gene_density = NULL, 
                              nuance = NULL,
                              center = TRUE) {
  m <- convert_to_matrix(gm)
  invalid_values <- is.infinite(m) | is.na(m)
  if (sum(invalid_values) > 0 & !is.null(nuance))
    m[invalid_values] <- rnorm(sum(invalid_values), sd = nuance)
  
  comp <-
    prcomp(
      # if (input_as_cor) m else cor(m, use = "pairwise.complete.obs"),
      if (input_as_cor) m else cor(m, use = "pairwise.complete.obs", method = "pearson"),
           center = TRUE
           # scale. = TRUE
      )$rotation[, 1]
  
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
    if (cor(comp, gd_cnt) < 0)
      comp <- -comp
  }
  
  result <- list(
    start = GenomicRanges::start(gr) - 1, 
    end = GenomicRanges::end(gr), 
    score = comp) %>%
    as_tibble()
  
  if (center) {
    return(result %>% mutate(score = score - mean(score, na.rm = TRUE)))
  } else {
    return(result)
  }
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