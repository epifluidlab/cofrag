plot_genomic_matrix <- function(
  gm, scale_factor = c(1, 1), symfill = TRUE, color_palette = "viridis",
  n.breaks = NULL,
  xlab = NULL,
  ylab = NULL,
  ggtitle = NULL) {
  if (symfill) {
    gm <- rbind(
      as_tibble(contact_matrix),
      as_tibble(contact_matrix) %>% filter(start1 != start2) %>%
        mutate(
          s = start1 + start2,
          start1 = (s - start1),
          start2 = (s - start2)
        ) %>% select(-s)
    )
  }
  ggplot(gm %>% mutate(start1 = start1 / scale_factor[1], start2 = start2 / scale_factor[2]),
    aes(x = start1, y = start2, fill = score)) +
    geom_tile(color = "black") + 
    scale_x_continuous(labels = scales::number_format(accuracy = 1), n.breaks = n.breaks, position = "top") + 
    scale_y_reverse(labels = scales::number_format(accuracy = 1), n.breaks = n.breaks) + 
    scale_fill_viridis_c(option = color_palette) + xlab + ylab + ggtitle
}