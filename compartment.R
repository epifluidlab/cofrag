source(here::here("genomic_matrix.R"))


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
      GenomicRanges::seqnames(interval),
      GenomicRanges::start(interval) - 1,
      GenomicRanges::end(interval)
    )
  })
}


# Perform the A/B compartment analysis
# Parameters:
#   * gm: a genomic_matrix object
#   * gene_density: a BED-format data frame of gene annotation
# Return: a vector indicating the compartment level
call_compartments <- function(gm, gene_density = NULL) {
  m <- convert_to_matrix(gm)
  comp <-
    prcomp(cor(m, use = "pairwise.complete.obs"),
           center = TRUE,
           scale. = TRUE)$rotation[, 1]
  
  # Compare the compartment level with gene density and determine whether to flip
  if (!is.null(gene_density)) {
    gr_start <- attr(gm, "gr_start")
    gr_end <- attr(gm, "gr_end")
    bin_size <- attr(gm, "bin_size")
    
    gr <- 1:((gr_end - gr_start) %/% bin_size) %>%
      map_chr(function(idx) {
        chr <- as.character(attr(gm, "chr"))
        start <- gr_start + (idx - 1) * bin_size
        str_interp("${chr}:$[d]{start + 1}-$[d]{start + bin_size}")
      }) %>%
      GenomicRanges::GRanges()
    
    gd_cnt <- .calc_num_genes(gr, gene_density)
    
    if (cor(comp, gd_cnt) < 0)
      comp <- -comp
  }
  
  comp
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
  gene_density <- read_tsv(
    conn,
    col_names = c("chr", "start", "end", "gene"),
    col_types = list(col_factor(), col_integer(), col_integer(), col_factor())
  )
  
  logdebug("Loading genomic matrix...")
  gm <- load_genomic_matrix(conn = fifo("/dev/stdin"))
  
  logdebug("Calling compartments...")
  comp <- call_compartments(gm = gm, gene_density = gene_density)
  
  # compartment output as BED
  gr_start <- attr(gm, "gr_start")
  gr_end <- attr(gm, "gr_end")
  chr <- attr(gm, "chr")
  bin_size <- attr(gm, "bin_size")
  
  gm_size <- (gr_end - gr_start) %/% bin_size
  contents <- 1:gm_size %>%
    map_chr(function(idx) {
      start <- gr_start + bin_size * (idx - 1) + 1
      end <- start + bin_size
      str_interp("${chr}\t$[d]{start}\t$[d]{end}\t${comp[idx]}")
    })
  writeLines(contents)
}