plot_genomic_matrix <- function(
  gm, scale_factor = c(1, 1), symfill = TRUE, missing_value = NULL, color_palette = "viridis",
  n.breaks = NULL,
  xlab = NULL,
  ylab = NULL,
  ggtitle = NULL) {
  # Replace missing values by ...
  if (is.null(missing_value)) {
    new_gm <- gm
  } else {
    gr <- attr(gm, "gr")
    bin_size <- attr(gm, "bin_size")
    gr_start1 <- (GenomicRanges::start(gr) - 1)[1]
    gr_start2 <- (GenomicRanges::start(gr) - 1)[2]
    n_bins <- GenomicRanges::width(gr)[1] %/% bin_size
    
    new_gm <- expand_grid(start1 = seq.int(gr_start1, by = bin_size, length.out = n_bins),
                          start2 = seq.int(gr_start1, by = bin_size, length.out = n_bins)) %>%
      filter(start1 <= start2)
    new_gm$score <- new_gm %>% apply(1, function(v) {
      res <- gm %>% filter(start1 == v[1] & start2 == v[2])
      if (nrow(res) == 0) {
        missing_value
      } else {
        res[1,]$score
      }
    })
  }
  
  if (symfill) {
    new_gm <- rbind(
      as_tibble(new_gm),
      as_tibble(new_gm) %>% filter(start1 != start2) %>%
        mutate(
          s = start1 + start2,
          start1 = (s - start1),
          start2 = (s - start2)
        ) %>% select(-s)
    )
  }
  
  ggplot(new_gm %>% mutate(start1 = start1 / scale_factor[1], start2 = start2 / scale_factor[2]),
    aes(x = start1, y = start2, fill = score)) +
    geom_tile(color = "black") + 
    scale_x_continuous(labels = scales::number_format(accuracy = 1), n.breaks = n.breaks, position = "top") + 
    scale_y_reverse(labels = scales::number_format(accuracy = 1), n.breaks = n.breaks) + 
    scale_fill_viridis_c(option = color_palette) + xlab + ylab + ggtitle
}