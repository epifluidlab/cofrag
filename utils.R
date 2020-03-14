library(gplots)

default_tick_func <- function(ticks)
  sapply(round(ticks / 1000), function(v)
    sprintf("%s kb", format(v, big.mark = ",")))

plot_hic_matrix <-
  function(m,
           start,
           end,
           tick_func = default_tick_func,
           color = "hic",
           main_title = NULL,
           ...) {
    my.image <- function(figData,  zlim, col, na.color = 'gray', ...)
    {
      newz.na <-
        zlim[2] + (zlim[2] - zlim[1]) / length(col) # new z for NA
      figData[which(is.na(figData))] <-
        newz.na # we affect newz.outside
      zlim[2] <-
        zlim[2] + (zlim[2] - zlim[1]) / length(col) # we finally extend the z limits to include the two new values
      col <-
        c(col, na.color) # we construct the new color range by including: na.color and outside.color
      image(figData,  zlim = zlim, col = col, ...) # we finally call image(...)
    }
    
    if (identical(color, "hic"))
      cp <- colorpanel(64, "#FB040A", "#FF7F82", "white")
    else
      cp <- color

    min_val <- min(as.vector(dm[!is.na(dm)]))
    max_val <- max(as.vector(dm[!is.na(dm)]))
    zlim <- c(min_val, max_val)
    my.image(
      m[, nrow(m):1],
      zlim,
      cp,
      na.color = "gray",
      useRaster = TRUE,
      axes = FALSE,
      ...
    )
    
    ticks <- tick_func(seq(start, end, length.out = 6))
    axis(
      2,
      at = seq(0, 1, length.out = 6),
      labels = rev(ticks),
      srt = 45,
      tick = TRUE
    )
    axis(
      3,
      at = seq(0, 1, length.out = 6),
      labels = ticks,
      srt = 45,
      tick = TRUE
    )
    if (is.null(main_title))
      main_title <-
      sprintf(
        "Hi-C heatmap for %s-%s",
        format(start, digits = 15, big.mark = ","),
        format(end, digits = 15, big.mark = ",")
      )
    title(main_title, line = 2.5)
  }
