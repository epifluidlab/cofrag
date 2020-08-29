plot_genomic_matrix <- function(gm,
                                control_gm = NULL,
                                scale_factor = c(1, 1),
                                # symfill = TRUE,
                                missing_value = NULL,
                                n.breaks = NULL,
                                color_palette = "viridis",
                                gamma = 1,
                                tile_outline = NULL) {
  gr <- attr(gm, "gr")
  bin_size <- attr(gm, "bin_size")
  
  # Replace missing values by ...
  if (is.null(missing_value)) {
    new_gm <- gm
  } else {
    gr_start1 <- (GenomicRanges::start(gr) - 1)[1]
    gr_start2 <- (GenomicRanges::start(gr) - 1)[2]
    n_bins <- GenomicRanges::width(gr)[1] %/% bin_size
    
    new_gm <-
      expand_grid(
        start1 = seq.int(gr_start1, by = bin_size, length.out = n_bins),
        start2 = seq.int(gr_start1, by = bin_size, length.out = n_bins)
      ) %>%
      filter(start1 <= start2)
    new_gm$score <- new_gm %>% apply(1, function(v) {
      res <- gm %>% filter(start1 == v[1] & start2 == v[2])
      if (nrow(res) == 0) {
        missing_value
      } else {
        res[1, ]$score
      }
    })
  }
  
  build_full_gm <- function(gm1, gm2) {
    rbind(
      as_tibble(gm1),
      as_tibble(gm2) %>% filter(start1 != start2) %>%
        mutate(
          s = start1 + start2,
          start1 = (s - start1),
          start2 = (s - start2)
        ) %>% select(-s)
    )
  }
  
  if (is.null(control_gm)) {
    new_gm <- build_full_gm(gm1 = gm, gm2 = new_gm)
  } else {
    # Ensure the dimensions are compatible
    # stopifnot(all(attr(gm, "gr") == attr(control_gm, "gr")))
    stopifnot(attr(gm, "bin_size") == attr(control_gm, "bin_size"))
    new_gm <- build_full_gm(gm1 = control_gm, gm2 = gm)
  }
  
  hic_colors <- list(
    colors = c(
      "#FFFFFF",
      "#FFF2F2",
      "#FFE8E8",
      "#FFCBCB",
      "#FFB3B3",
      "#FFA4A4",
      "#FF6565",
      "#FF0402"
    ),
    values = (c(0, 56, 95, 218, 265, 369, 603, 1033) / 1033) ** 0.45
  )
  
  ggplot(
    new_gm %>% mutate(
      start1 = start1 / scale_factor[1],
      start2 = start2 / scale_factor[2]
    ),
    aes(x = start1, y = start2, fill = score)
  ) +
    (if (is.null(tile_outline))
      geom_tile()
     else
       geom_tile(color = tile_outline)) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 1),
      n.breaks = n.breaks,
      position = "top"
    ) +
    scale_y_reverse(labels = scales::number_format(accuracy = 1),
                    n.breaks = n.breaks) +
    scale_fill_gradientn(colors = hic_colors$colors, values = hic_colors$values ** gamma) +
    ylab(paste0(
      ifelse(is.null(control_gm), "", "Control: "), as.character(gr)[1])) + 
    xlab(as.character(gr)[2])
}


# Plot the A/B compartment profile
plot_compartments <- function(data) {
  resolution <- data %>% .$start %>% diff() %>% min()
  data %>%
    mutate(compartment = factor(
      ifelse(score > 0, "open", "close"), levels = c("open", "close"))) %>%
    
    ggplot(aes(x = start, y = score, fill = compartment)) +
    geom_col(width = resolution * .75) +
    scale_x_continuous(labels = scales::comma) +
    xlab("coordinate")
}