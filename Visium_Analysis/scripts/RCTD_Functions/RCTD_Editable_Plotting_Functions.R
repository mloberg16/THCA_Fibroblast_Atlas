### Author: Matthew Aaron Loberg

### Function Notes:
# Creating plot_puck_continuous_edited as a function
# The major change here is to replace size with size = size
# When running RCTD, I will call "plot_puck_continuous_edited" instead of "plot_puck_continuous"

plot_puck_continuous_edited <- function (puck, barcodes, plot_val, ylimit = c(0, 1), title = NULL, 
                                         counter_barcodes = NULL, label = F, my_pal = NULL, xlim = NULL, 
                                         ylim = NULL, size = 2, alpha = 1, small_point = F) 
{
  if (is.null(my_pal)) 
    my_pal = pals::kovesi.rainbow(20)
  my_table = puck@coords[barcodes, ]
  plot_val <- pmax(pmin(plot_val, ylimit[2] - 1e-08), ylimit[1] + 
                     1e-08)
  my_table$value = plot_val[barcodes]
  if (!is.null(ylimit)) 
    sc <- ggplot2::scale_colour_gradientn(colors = my_pal, 
                                          limits = ylimit)
  else sc <- ggplot2::scale_colour_gradientn(colors = my_pal)
  plot <- ggplot2::ggplot(my_table, ggplot2::aes(x = x, y = y))
  if (small_point) 
    plot <- plot + ggplot2::geom_point(ggplot2::aes(size = size, 
                                                    shape = 16, color = value, stroke = 0), alpha = alpha)
  else plot <- plot + ggplot2::geom_point(ggplot2::aes(size = size, 
                                                       shape = 19, color = value), alpha = alpha)
  plot <- plot + sc + ggplot2::scale_shape_identity() + ggplot2::theme_classic() + 
    ggplot2::scale_size_identity()
  if (label) 
    plot <- plot + ggplot2::aes(label = which(rownames(puck@coords) %in% 
                                                inter_barcodes[my_class == i])) + ggplot2::geom_text()
  if (!is.null(counter_barcodes)) {
    my_table = puck@coords[counter_barcodes, ]
    plot <- plot + ggplot2::geom_point(data = my_table, 
                                       ggplot2::aes(x = x, y = y, size = size), alpha = 0.1)
  }
  if (is.null(xlim)) 
    xlim <- c(min(puck@coords$x) - 1, max(puck@coords$x) + 
                1)
  if (is.null(ylim)) 
    ylim <- c(min(puck@coords$y) - 1, max(puck@coords$y) + 
                1)
  plot <- plot + ggplot2::coord_fixed() + ggplot2::xlim(xlim) + 
    ggplot2::ylim(ylim)
  if (!is.null(title)) 
    plot <- plot + ggplot2::ggtitle(title)
  plot
}

# Editable plot_puck_wrapper
plot_puck_wrapper_edited <- function (puck, plot_val, cell_type = NULL, minUMI = 0, maxUMI = 2e+05, 
          min_val = NULL, max_val = NULL, title = NULL, my_cond = NULL) 
{
  UMI_filter = (puck@nUMI > minUMI) & (puck@nUMI < maxUMI)
  ylimit = NULL
  if (!is.null(my_cond)) 
    my_cond = UMI_filter & my_cond
  else my_cond = UMI_filter
  if (!is.null(cell_type)) 
    my_cond = my_cond & (puck@cell_labels == cell_type)
  if (!is.null(min_val)) 
    my_cond = my_cond & (plot_val > min_val)
  if (!is.null(max_val)) {
    epsilon = 1e-05
    plot_val[plot_val >= max_val - epsilon] = max_val - 
      epsilon
    if (!is.null(min_val)) 
      ylimit = c(min_val, max_val)
  }
  plot_puck_continuous_edited(puck, names(which(my_cond)), plot_val, 
                       title = title, ylimit = ylimit)
}

# Plot weights unthreshold edited
plot_weights_unthreshold_edited <- function (cell_type_names, puck, resultsdir, weights) 
{
  plots <- vector(mode = "list", length = length(cell_type_names))
  for (i in 1:length(cell_type_names)) {
    cell_type = cell_type_names[i]
    plot_var <- weights[, cell_type]
    names(plot_var) = rownames(weights)
    if (sum(weights[, cell_type]) > 0) 
      plots[[i]] <- plot_puck_wrapper_edited(puck, plot_var, 
                                      NULL, minUMI = 100, maxUMI = 2e+05, min_val = 0, 
                                      max_val = 1, title = cell_type)
  }
  pdf(file.path(resultsdir, "cell_type_weights_unthreshold.pdf"))
  invisible(lapply(plots, print))
  dev.off()
}
