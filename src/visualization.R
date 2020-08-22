plot_genomic_matrix <- function(
  gm, scale_factor = c(1, 1), symfill = TRUE, missing_value = NULL, 
  n.breaks = NULL,
  color_palette = "viridis",
  gamma = 1,
  tile_outline = NULL
  ) {
  gr <- attr(gm, "gr")
  bin_size <- attr(gm, "bin_size")
  
  # Replace missing values by ...
  if (is.null(missing_value)) {
    new_gm <- gm
  } else {
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
  
  hic_colors <- list(
    colors = c("#FFFFFF", "#FFF2F2", "#FFE8E8", "#FFCBCB", "#FFB3B3", "#FFA4A4", "#FF6565", "#FF0402"), 
    values = (c(0, 56, 95, 218, 265, 369, 603, 1033)/1033) ** 0.45)
  
  ggplot(new_gm %>% mutate(start1 = start1 / scale_factor[1], start2 = start2 / scale_factor[2]),
    aes(x = start1, y = start2, fill = score)) +
    (if (is.null(tile_outline)) geom_tile() else geom_tile(color = tile_outline)) +
    # ifelse(is.null(tile_outline), geom_tile(), geom_tile(color = tile_outline)) +
    # geom_tile(color = tile_outline) + 
    # geom_tile(color = "black") + 
    scale_x_continuous(labels = scales::number_format(accuracy = 1), n.breaks = n.breaks, position = "top") + 
    scale_y_reverse(labels = scales::number_format(accuracy = 1), n.breaks = n.breaks) + 
    scale_fill_gradientn(colors = hic_colors$colors, values = hic_colors$values ** gamma) +
    # scale_fill_viridis_c(values = (seq(0, 10) / 10) ** gamma, option = color_palette) +
    ylab(as.character(gr)[1]) + xlab(as.character(gr)[2])
}