

plot_fits <- function(model,level = c("ind","group","global"), show_ci = TRUE, facet = show_ci, show_events = TRUE,
                      palette = "Purples", ...) {
  level <- match(level)
  tidy_model <- reconstruct_traces(model,level = level, ...)
  if (any(c("one_ind","one_group")) %in% class(model) && level == "global") {
    stop("Global trace not present for models with one individual or group")
  }

  if (!show_ci) {
    p <- ggplot(tidy_model, aes(x = x, y = y, ymin = ymin, ymax = ymax, fill = factor(.width),group = ind))
  } else {
    p <- p + ggdist::geom_lineribbon(alpha=0.6,linewidth=0.5) +
      scale_fill_brewer(palette = palette, direction = -1)
  }

  if (!is.null(model$traces) && show_traces) {
      p <- p + geom_line(data=model$traces,aes(x=x,y=y+0.5,group=.data[[level]]),
                         inherit.aes=FALSE,linewidth=1)
    }

  if (facet) {
    p <- p + facet_wrap(~.data[[level]])
  }


}

    #facet_wrap(~ind)


