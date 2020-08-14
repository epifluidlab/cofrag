.build_msg_coloring <- function() {
  crayon_env <- tryCatch(asNamespace("crayon"),
                         error = function(e) NULL)
  
  default_color_msg <- function(msg, level_name) msg
  if (is.null(crayon_env)) {
    return(default_color_msg)
  }
  
  if (is.null(crayon_env$make_style) ||
      is.null(crayon_env$combine_styles) ||
      is.null(crayon_env$reset)) {
    return(default_color_msg)
  }
  
  color_msg <- function(msg, level_name) {
    style <- switch(level_name,
                    "FINEST" = crayon_env$make_style("gray80"),
                    "FINER" = crayon_env$make_style("gray60"),
                    "FINE" = crayon_env$make_style("gray60"),
                    "DEBUG" = crayon_env$make_style("deepskyblue4"),
                    "INFO" = crayon_env$reset,
                    "WARNING" = crayon_env$make_style("darkorange"),
                    "ERROR" = crayon_env$make_style("red4"),
                    "CRITICAL" =
                      crayon_env$combine_styles(crayon_env$bold,
                                                crayon_env$make_style("red1")),
                    crayon_env$make_style("gray100"))
    res <- paste0(style(msg), crayon_env$reset(""))
    return(res)
  }
  return(color_msg)
}

writeToConnFactory <- function(conn) {
  writeToConn <- function(msg, handler, ...) {
    if (length(list(...)) && "dry" %in% names(list(...))) {
      if (!is.null(handler$color_output) &&
          handler$color_output == FALSE) {
        handler$color_msg <- function(msg, level_name)
          msg
      } else {
        handler$color_msg <- .build_msg_coloring()
      }
      return(TRUE)
    }
    
    stopifnot(length(list(...)) > 0)
    
    level_name <- list(...)[[1]]$levelname
    msg <- handler$color_msg(msg, level_name)
    cat(paste0(msg, "\n"), file = conn)
  }
  return(writeToConn)
}