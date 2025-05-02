#' 2D Scatter Plot of Sender and Receiver Pathway-Level Changes Across Multiple Groups
#'
#' @description
#' Creates a 2D scatter plot showing how cell types change in both sender and receiver
#' dimensions across multiple conditions, using pathway-level information from the netP slot.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (NULL for all objects)
#' @param reference Integer or character, index or name of the reference object to compare against.
#'   If NULL (default), the first object in comparison will be used as reference.
#' @param pathways Vector of pathway names to include (NULL for all)
#' @param thresh P-value threshold for significant interactions (default: 0.05)
#' @param log.scale Whether to use log scale for axes (default: TRUE)
#' @param arrow.size Size of trajectory arrows (default: 1)
#' @param arrow.alpha Transparency of trajectory arrows (default: 0.8)
#' @param arrow.color Color of arrows (default: NULL, uses cell type colors)
#' @param arrow.type Type of arrow head, 'open' or 'closed' (default: "closed")
#' @param label.cell Whether to label cell types (default: TRUE)
#' @param label.size Size of cell type labels (default: 3)
#' @param label.color Color of cell type labels (default: "black")
#' @param label.repel Whether to use ggrepel for label placement (default: TRUE)
#' @param label.box Whether to add box around labels (default: FALSE)
#' @param label.min.segment.length Minimum segment length for repelled labels (default: 0)
#' @param label.max.overlaps Maximum allowed label overlaps (default: 10)
#' @param add.quadrants Whether to add quadrant lines and labels (default: TRUE)
#' @param quadrant.line.type Type of quadrant dividing lines (default: "dashed")
#' @param quadrant.line.color Color of quadrant dividing lines (default: "grey50")
#' @param quadrant.line.alpha Alpha of quadrant dividing lines (default: 0.5)
#' @param quadrant.label.size Size of quadrant labels (default: 3)
#' @param quadrant.label.alpha Alpha of quadrant labels (default: 0.7)
#' @param point.size.range Range of point sizes (default: c(2, 6))
#' @param point.alpha.values Vector of alpha values for conditions (default: c(0.7, 1))
#' @param margin.factor.x Horizontal expansion factor for plot margins (default: 0.3)
#' @param margin.factor.y Vertical expansion factor for plot margins (default: 0.3)
#' @param colors Vector of colors for cell types (default: NULL, uses global_colors)
#' @param group.colors Vector of colors for groups/conditions (default: NULL)
#' @param group.shapes Vector of shapes for different conditions (default: NULL, auto-generated)
#' @param use.group.colors Whether to color points by group instead of cell type (default: FALSE)
#' @param show.group.legend Whether to show the group legend (default: TRUE)
#' @param convex.hull Whether to draw convex hulls around points of the same group (default: FALSE)
#' @param hull.alpha Transparency of convex hulls (default: 0.2)
#' @param group.label.size Size of group labels (default: 3.5)
#' @param show.group.labels Whether to add labels for each group (default: FALSE)
#' @param title Plot title (default: "Changes in Sender-Receiver Roles")
#' @param title.size Size of plot title (default: 14)
#' @param title.face Font face of plot title (default: "bold")
#' @param axis.text.size Size of axis text (default: 10)
#' @param axis.title.size Size of axis titles (default: 12)
#' @param legend.position Position of the legend (default: "right")
#' @param base.theme Base theme for the plot (default: ggplot2::theme_bw())
#' @param plot.type Type of plot: "combined" (all trajectories on one plot),
#'   "facet" (separate plot for each condition vs reference),
#'   "grid" (separate plots arranged in a grid), or
#'   "group_colored" (all trajectories on one plot with group coloring)
#' @param facet.ncol Number of columns when using faceted plot
#' @param highlight.cells Vector of cell types to highlight
#' @param highlight.color Color for highlighted cells
#' @param comparison_method Method for comparing groups: "all_vs_ref" (each group vs reference)
#'   or "sequential" (each group vs next) or "custom_pairs" (specified pairs)
#' @param custom_comparisons List of custom comparison pairs
#' @param common.scale Whether to use a common scale for all plots
#' @param save_plot Whether to save the plot to file
#' @param save_name Filename for saving the plot
#' @param save_width Width of the saved plot in inches
#' @param save_height Height of the saved plot in inches
#' @param return_data Whether to return the data along with the plot
#'
#' @return A ggplot object or a list with plot and data if return_data=TRUE
#' @importFrom ggplot2 ggplot aes geom_point geom_segment theme_bw labs scale_color_manual
#'   scale_shape_manual scale_alpha_manual scale_size annotate xlim ylim theme element_text
#'   facet_wrap ggtitle
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by arrange slice filter
#' @importFrom gridExtra grid.arrange
#' @export
scatterDiff2DM <- function(object.list, comparison = NULL, reference = NULL,
                           pathways = NULL, thresh = 0.05,
                           log.scale = TRUE,
                           arrow.size = 1, arrow.alpha = 0.8,
                           arrow.color = NULL, arrow.type = "closed",
                           label.cell = TRUE, label.size = 3,
                           label.color = "black", label.repel = TRUE,
                           label.box = FALSE,
                           label.min.segment.length = 0,
                           label.max.overlaps = 10,
                           add.quadrants = TRUE,
                           quadrant.line.type = "dashed",
                           quadrant.line.color = "grey50",
                           quadrant.line.alpha = 0.5,
                           quadrant.label.size = 3,
                           quadrant.label.alpha = 0.7,
                           point.size.range = c(2, 6),
                           point.alpha.values = c(0.7, 1),
                           margin.factor.x = 0.3,
                           margin.factor.y = 0.3,
                           colors = NULL,
                           group.colors = NULL,
                           group.shapes = NULL,
                           use.group.colors = FALSE,
                           show.group.legend = TRUE,
                           convex.hull = FALSE,
                           hull.alpha =  0.2,
                           group.label.size = 3.5,
                           show.group.labels = FALSE,
                           title = "Changes in Sender-Receiver Roles",
                           title.size = 14, title.face = "bold",
                           axis.text.size = 10, axis.title.size = 12,
                           legend.position = "right",
                           base.theme = ggplot2::theme_bw(),
                           plot.type = c("combined", "facet", "grid", "group_colored"),
                           facet.ncol = 2,
                           highlight.cells = NULL,
                           highlight.color = "red",
                           comparison_method = c("all_vs_ref", "sequential", "custom_pairs"),
                           custom_comparisons = NULL,
                           common.scale = TRUE,
                           save_plot = FALSE, save_name = "scatterDiff2DM.pdf",
                           save_width = 10, save_height = 8,
                           return_data = FALSE) {

  # Match arguments
  plot.type <- match.arg(plot.type)
  comparison_method <- match.arg(comparison_method)

  # For backward compatibility - if plot.type is "group_colored", set use.group.colors to TRUE
  if (plot.type == "group_colored") {
    use.group.colors <- TRUE
  }

  # If comparison is NULL, use all objects
  if (is.null(comparison)) {
    comparison <- 1:length(object.list)
  } else if (is.character(comparison)) {
    # If comparison contains names, match them to indices
    if (all(comparison %in% names(object.list))) {
      comparison <- match(comparison, names(object.list))
    } else {
      stop("If comparison contains names, all must match names in object.list")
    }
  }

  # Process reference
  if (is.null(reference)) {
    reference <- comparison[1]
    message("Using first object in comparison as reference")
  } else if (is.character(reference) && length(reference) == 1) {
    if (reference %in% names(object.list)) {
      reference <- match(reference, names(object.list))
    } else {
      warning("Reference name not found in object.list. Using first object as reference.")
      reference <- comparison[1]
    }
  }

  if (!reference %in% comparison) {
    warning("Reference index not in comparison. Using first object in comparison as reference.")
    reference <- comparison[1]
  }

  # Get condition names
  condition_names <- names(object.list)
  if (is.null(condition_names) || any(condition_names == "")) {
    condition_names <- paste0("Condition_", 1:length(object.list))
    names(object.list) <- condition_names
  }

  # Extract cell types
  cell_types_list <- lapply(object.list[comparison], function(obj) {
    rownames(obj@netP$prob)
  })
  common_cell_types <- Reduce(intersect, cell_types_list)

  if (length(common_cell_types) == 0) {
    stop("No common cell types found across all conditions")
  }

  # Get all pathways if not specified
  if (is.null(pathways)) {
    pathways_list <- lapply(object.list[comparison], function(obj) {
      dimnames(obj@netP$prob)[[3]]
    })
    pathways <- Reduce(intersect, pathways_list)

    if (length(pathways) == 0) {
      stop("No common pathways found across all conditions")
    }
  } else {
    # Verify that all specified pathways exist in all objects
    for (idx in comparison) {
      obj_pathways <- dimnames(object.list[[idx]]@netP$prob)[[3]]
      missing_pathways <- setdiff(pathways, obj_pathways)
      if (length(missing_pathways) > 0) {
        warning(paste("Pathways", paste(missing_pathways, collapse = ", "),
                      "not found in condition", condition_names[idx]))
      }
    }
    # Filter to only include pathways present in all objects
    common_pathways <- lapply(object.list[comparison], function(obj) {
      intersect(pathways, dimnames(obj@netP$prob)[[3]])
    })
    pathways <- Reduce(intersect, common_pathways)

    if (length(pathways) == 0) {
      stop("None of the specified pathways were found in all conditions")
    }
  }

  # Determine which comparisons to make
  comparison_pairs <- list()

  if (comparison_method == "all_vs_ref") {
    # All conditions vs reference
    ref_idx <- which(comparison == reference)
    for (i in 1:length(comparison)) {
      if (i != ref_idx) {
        comparison_pairs[[length(comparison_pairs) + 1]] <- c(ref_idx, i)
      }
    }
  } else if (comparison_method == "sequential") {
    # Each condition vs the next
    for (i in 1:(length(comparison)-1)) {
      comparison_pairs[[i]] <- c(i, i+1)
    }
  } else if (comparison_method == "custom_pairs") {
    # Use custom comparison pairs
    if (is.null(custom_comparisons) || !is.list(custom_comparisons)) {
      stop("For custom_pairs comparison method, please provide a list of comparison pairs")
    }

    # Validate custom comparisons
    for (i in 1:length(custom_comparisons)) {
      pair <- custom_comparisons[[i]]
      if (length(pair) != 2) {
        stop(paste("Each comparison pair must have exactly 2 elements. Pair", i, "has",
                   length(pair), "element(s)"))
      }

      # Convert names to indices if needed
      if (is.character(pair)) {
        if (!all(pair %in% names(object.list))) {
          stop(paste("Invalid condition name(s) in pair", i))
        }
        pair <- match(pair, names(object.list))
      }

      # Check if indices are valid
      if (!all(pair %in% comparison)) {
        stop(paste("Indices in pair", i, "must be in the comparison vector"))
      }

      comparison_pairs[[i]] <- pair
    }
  }

  if (length(comparison_pairs) == 0) {
    stop("No valid comparison pairs were generated")
  }

  # Set cell type colors if not provided
  if (is.null(colors)) {
    if (exists("global_colors") && length(global_colors) >= length(common_cell_types)) {
      colors <- global_colors[1:length(common_cell_types)]
    } else {
      colors <- rainbow(length(common_cell_types))
    }
    names(colors) <- common_cell_types
  } else if (length(colors) < length(common_cell_types)) {
    # Recycle colors if needed
    colors <- rep(colors, length.out = length(common_cell_types))
    names(colors) <- common_cell_types
  }

  # Set group colors if not provided
  if (is.null(group.colors)) {
    if (exists("global_colors") && length(global_colors) >= length(comparison)) {
      group.colors <- global_colors[1:length(comparison)]
    } else {
      group.colors <- rainbow(length(comparison))
    }
    names(group.colors) <- condition_names[comparison]
  } else if (length(group.colors) < length(comparison)) {
    # Recycle colors if needed
    group.colors <- rep(group.colors, length.out = length(comparison))
    names(group.colors) <- condition_names[comparison]
  }

  # Set default shapes for groups if not provided
  if (is.null(group.shapes)) {
    # Common distinctive shapes
    shape_options <- c(16, 17, 15, 18, 8, 7, 9, 10, 11, 12, 13, 14)
    group.shapes <- shape_options[1:length(comparison)]
    names(group.shapes) <- condition_names[comparison]
  } else if (length(group.shapes) < length(comparison)) {
    # Recycle shapes if needed
    group.shapes <- rep(group.shapes, length.out = length(comparison))
    names(group.shapes) <- condition_names[comparison]
  }

  # Initialize lists to store results
  all_scores <- list()
  all_data <- list()

  # Get signaling scores for all conditions
  for (i in 1:length(comparison)) {
    idx <- comparison[i]
    obj <- object.list[[idx]]

    # Calculate sender and receiver scores
    scores <- calculateScores(obj, common_cell_types, pathways)
    all_scores[[i]] <- scores

    # Get interaction counts
    counts <- getInteractionCounts(obj, signaling = pathways,
                                   slot.name = "net", thresh = thresh)

    # Create data frame for this condition
    data <- data.frame(
      celltype = common_cell_types,
      sender = scores$sender,
      receiver = scores$receiver,
      condition = condition_names[idx],
      condition_idx = idx,
      count = counts[common_cell_types],
      shape = group.shapes[condition_names[idx]]
    )

    # Apply log transformation if requested
    if (log.scale) {
      data$sender <- log1p(data$sender)
      data$receiver <- log1p(data$receiver)
    }

    all_data[[i]] <- data
  }

  # Combine all condition data
  combined_data <- do.call(rbind, all_data)

  # Find global min/max for consistent scales if requested
  if (common.scale) {
    x_range <- range(combined_data$sender, na.rm = TRUE)
    y_range <- range(combined_data$receiver, na.rm = TRUE)

    # Expand the range by user-defined factors
    x_span <- diff(x_range)
    y_span <- diff(y_range)

    x_min <- x_range[1] - margin.factor.x * x_span
    x_max <- x_range[2] + margin.factor.x * x_span
    y_min <- y_range[1] - margin.factor.y * y_span
    y_max <- y_range[2] + margin.factor.y * y_span
  }

  # Create data for each comparison pair
  all_plot_data <- list()
  all_plots <- list()

  for (p in 1:length(comparison_pairs)) {
    pair <- comparison_pairs[[p]]
    idx1 <- comparison[pair[1]]
    idx2 <- comparison[pair[2]]

    # Get data for both conditions
    data1 <- all_data[[pair[1]]]
    data2 <- all_data[[pair[2]]]

    # Create arrows data showing direction of change
    arrows_data <- data.frame(
      celltype = common_cell_types,
      x = data1$sender,
      y = data1$receiver,
      xend = data2$sender,
      yend = data2$receiver,
      condition1 = data1$condition,
      condition2 = data2$condition,
      condition1_idx = idx1,
      condition2_idx = idx2,
      count = pmax(data1$count, data2$count),
      shape1 = data1$shape,
      shape2 = data2$shape
    )

    # Calculate changes
    arrows_data$sender_change <- arrows_data$xend - arrows_data$x
    arrows_data$receiver_change <- arrows_data$yend - arrows_data$y
    arrows_data$total_change <- sqrt(arrows_data$sender_change^2 +
                                       arrows_data$receiver_change^2)

    # Determine direction
    arrows_data$sender_direction <- ifelse(arrows_data$sender_change > 0,
                                           "increase", "decrease")
    arrows_data$receiver_direction <- ifelse(arrows_data$receiver_change > 0,
                                             "increase", "decrease")

    # Assign quadrants
    arrows_data$quadrant <- ifelse(
      arrows_data$sender_change > 0 & arrows_data$receiver_change > 0, "I",
      ifelse(arrows_data$sender_change < 0 & arrows_data$receiver_change > 0, "II",
             ifelse(arrows_data$sender_change < 0 & arrows_data$receiver_change < 0, "III", "IV")))

    # Add comparison identifier
    arrows_data$comparison <- paste(data2$condition[1], "vs", data1$condition[1])

    # Highlight specific cells if requested
    if (!is.null(highlight.cells)) {
      arrows_data$highlighted <- arrows_data$celltype %in% highlight.cells
    } else {
      arrows_data$highlighted <- FALSE
    }

    # Store the data
    all_plot_data[[p]] <- arrows_data

    # Skip creating individual plots if we're doing a combined plot
    if (plot.type == "combined" && length(comparison_pairs) > 1) {
      next
    }

    # Create individual plot
    p_plot <- createScatterPlot2D(
      data1, data2, arrows_data,
      colors = colors,
      group.colors = group.colors,
      group.shapes = group.shapes,
      use.group.colors = use.group.colors,
      show.group.legend = show.group.legend,
      arrow.size = arrow.size,
      arrow.alpha = arrow.alpha,
      arrow.color = arrow.color,
      arrow.type = arrow.type,
      label.cell = label.cell,
      label.size = label.size,
      label.color = label.color,
      label.repel = label.repel,
      label.box = label.box,
      label.min.segment.length = label.min.segment.length,
      label.max.overlaps = label.max.overlaps,
      add.quadrants = add.quadrants,
      quadrant.line.type = quadrant.line.type,
      quadrant.line.color = quadrant.line.color,
      quadrant.line.alpha = quadrant.line.alpha,
      quadrant.label.size = quadrant.label.size,
      quadrant.label.alpha = quadrant.label.alpha,
      point.size.range = point.size.range,
      point.alpha.values = point.alpha.values,
      highlight.color = highlight.color,
      log.scale = log.scale,
      convex.hull = convex.hull,
      hull.alpha = hull.alpha,
      show.group.labels = show.group.labels,
      group.label.size = group.label.size,
      title = paste0(title, " - ", arrows_data$comparison[1]),
      title.size = title.size,
      title.face = title.face,
      axis.text.size = axis.text.size,
      axis.title.size = axis.title.size,
      legend.position = legend.position,
      base.theme = base.theme,
      pathways.used = length(pathways)
    )

    # Apply common scale if requested
    if (common.scale) {
      p_plot <- p_plot +
        ggplot2::xlim(c(x_min, x_max)) +
        ggplot2::ylim(c(y_min, y_max))
    }

    # Store the plot
    all_plots[[p]] <- p_plot
  }

  # Combine all arrow data
  combined_arrows <- do.call(rbind, all_plot_data)

  # Create legend title based on comparison method
  if (comparison_method == "all_vs_ref") {
    legend_title <- paste0("Group vs ", condition_names[reference])
  } else {
    legend_title <- "Group"
  }

  # Create final plot based on plot type
  if (plot.type == "combined" && length(comparison_pairs) > 1) {
    # All trajectories on one plot
    p <- createScatterPlot2D(
      combined_data, NULL, combined_arrows,
      colors = colors,
      group.colors = group.colors,
      group.shapes = group.shapes,
      use.group.colors = use.group.colors,
      show.group.legend = show.group.legend,
      arrow.size = arrow.size,
      arrow.alpha = arrow.alpha,
      arrow.color = arrow.color,
      arrow.type = arrow.type,
      label.cell = label.cell,
      label.size = label.size,
      label.color = label.color,
      label.repel = label.repel,
      label.box = label.box,
      label.min.segment.length = label.min.segment.length,
      label.max.overlaps = label.max.overlaps,
      add.quadrants = add.quadrants,
      quadrant.line.type = quadrant.line.type,
      quadrant.line.color = quadrant.line.color,
      quadrant.line.alpha = quadrant.line.alpha,
      quadrant.label.size = quadrant.label.size,
      quadrant.label.alpha = quadrant.label.alpha,
      point.size.range = point.size.range,
      point.alpha.values = point.alpha.values,
      highlight.color = highlight.color,
      log.scale = log.scale,
      convex.hull = convex.hull,
      hull.alpha = hull.alpha,
      show.group.labels = show.group.labels,
      group.label.size = group.label.size,
      legend_title = legend_title,
      title = title,
      title.size = title.size,
      title.face = title.face,
      axis.text.size = axis.text.size,
      axis.title.size = axis.title.size,
      legend.position = legend.position,
      base.theme = base.theme,
      pathways.used = length(pathways)
    )

    # Apply common scale if requested
    if (common.scale) {
      p <- p +
        ggplot2::xlim(c(x_min, x_max)) +
        ggplot2::ylim(c(y_min, y_max))
    }

  } else if (plot.type == "facet" && length(comparison_pairs) > 1) {
    # Use facets to separate comparisons
    p <- createScatterPlot2D(
      combined_data, NULL, combined_arrows,
      colors = colors,
      group.colors = group.colors,
      group.shapes = group.shapes,
      use.group.colors = use.group.colors,
      show.group.legend = show.group.legend,
      arrow.size = arrow.size,
      arrow.alpha = arrow.alpha,
      arrow.color = arrow.color,
      arrow.type = arrow.type,
      label.cell = label.cell,
      label.size = label.size,
      label.color = label.color,
      label.repel = label.repel,
      label.box = label.box,
      label.min.segment.length = label.min.segment.length,
      label.max.overlaps = label.max.overlaps,
      add.quadrants = add.quadrants,
      quadrant.line.type = quadrant.line.type,
      quadrant.line.color = quadrant.line.color,
      quadrant.line.alpha = quadrant.line.alpha,
      quadrant.label.size = quadrant.label.size,
      quadrant.label.alpha = quadrant.label.alpha,
      point.size.range = point.size.range,
      point.alpha.values = point.alpha.values,
      highlight.color = highlight.color,
      log.scale = log.scale,
      convex.hull = convex.hull,
      hull.alpha = hull.alpha,
      show.group.labels = show.group.labels,
      group.label.size = group.label.size,
      legend_title = legend_title,
      title = title,
      title.size = title.size,
      title.face = title.face,
      axis.text.size = axis.text.size,
      axis.title.size = axis.title.size,
      legend.position = legend.position,
      base.theme = base.theme,
      pathways.used = length(pathways),
      use_facets = TRUE,
      facet.ncol = facet.ncol
    )

    # Apply common scale if requested
    if (common.scale) {
      p <- p +
        ggplot2::xlim(c(x_min, x_max)) +
        ggplot2::ylim(c(y_min, y_max))
    }

  } else if (plot.type == "grid" && length(comparison_pairs) > 1) {
    # Use the individual plots we created earlier
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      valid_plots <- all_plots[!sapply(all_plots, is.null)]

      if (length(valid_plots) > 0) {
        # Arrange in a grid
        p <- do.call(
          gridExtra::grid.arrange,
          c(valid_plots,
            list(
              ncol = min(length(valid_plots), facet.ncol),
              top = title
            )
          )
        )
      } else {
        warning("No valid plots to arrange in a grid")
        p <- NULL
      }
    } else {
      warning("gridExtra package is required for 'grid' plot type. Please install it.")
      p <- all_plots[[1]]  # Fall back to the first plot
    }
  } else {
    # Just use the first individual plot
    p <- all_plots[[1]]
  }

  # Save plot if requested
  if (save_plot && !is.null(p)) {
    if (requireNamespace("ggplot2", quietly = TRUE) &&
        ("ggplot" %in% class(p) || plot.type != "grid")) {
      ggplot2::ggsave(
        filename = save_name,
        plot = p,
        width = save_width,
        height = save_height
      )
      message(paste("Plot saved to", save_name))
    } else if (requireNamespace("grDevices", quietly = TRUE)) {
      # For grid-arranged plots
      grDevices::pdf(save_name, width = save_width, height = save_height)
      print(p)
      grDevices::dev.off()
      message(paste("Plot saved to", save_name))
    } else {
      warning("Either ggplot2 or grDevices is required for saving plots. Please install them.")
    }
  }

  # Return plot and data if requested
  if (return_data) {
    return(list(
      plot = p,
      data = combined_data,
      arrows = combined_arrows,
      all_plots = all_plots,
      comparison_pairs = comparison_pairs
    ))
  } else {
    return(p)
  }
}

createScatterPlot2D <- function(data1, data2, arrows_data,
                                colors,
                                group.colors,
                                group.shapes,
                                use.group.colors = FALSE,
                                show.group.legend = TRUE,
                                arrow.size, arrow.alpha, arrow.color, arrow.type,
                                label.cell, label.size, label.color, label.repel,
                                label.box, label.min.segment.length, label.max.overlaps,
                                add.quadrants, quadrant.line.type, quadrant.line.color,
                                quadrant.line.alpha, quadrant.label.size, quadrant.label.alpha,
                                point.size.range, point.alpha.values,
                                highlight.color, log.scale,
                                convex.hull = FALSE, hull.alpha = 0.2,
                                show.group.labels = FALSE, group.label.size = 3.5,
                                legend_title = "Group",
                                title, title.size, title.face,
                                axis.text.size, axis.title.size, legend.position, base.theme,
                                pathways.used, use_facets = FALSE, facet.ncol = 2) {

  # Create the plot
  p <- ggplot2::ggplot()

  # Add convex hulls if requested
  if (convex.hull && !is.null(data1)) {
    # For each condition, create a polygon of its hull
    hulls <- lapply(unique(arrows_data$condition2), function(cond) {
      cond_data <- arrows_data[arrows_data$condition2 == cond, ]
      if (nrow(cond_data) > 2) {
        hull_indices <- grDevices::chull(cond_data$xend, cond_data$yend)
        hull_data <- cond_data[hull_indices, ]
        hull_data$hull_group <- cond
        return(hull_data)
      }
      return(NULL)
    })

    # Remove NULL entries
    hulls <- hulls[!sapply(hulls, is.null)]

    # Combine hull data
    if (length(hulls) > 0) {
      hull_data <- do.call(rbind, hulls)

      # Add the convex hull polygons - we'll make them a light gray
      p <- p + ggplot2::geom_polygon(
        data = hull_data,
        ggplot2::aes(
          x = xend,
          y = yend,
          group = hull_group
        ),
        alpha = hull.alpha,
        fill = "gray80"  # Use a neutral color for all hulls
      )
    }
  }

  # Add points - always color by celltype, but use shape for condition
  if (!is.null(data1)) {
    p <- p + ggplot2::geom_point(
      data = data1,
      ggplot2::aes(
        x = sender,
        y = receiver,
        color = celltype,
        shape = condition,
        size = count
      ),
      alpha = point.alpha.values[1]
    )
  }

  if (!is.null(data2)) {
    p <- p + ggplot2::geom_point(
      data = data2,
      ggplot2::aes(
        x = sender,
        y = receiver,
        color = celltype,
        shape = condition,
        size = count
      ),
      alpha = point.alpha.values[2]
    )
  }

  # Add arrows - color by celltype
  p <- p + ggplot2::geom_segment(
    data = arrows_data,
    ggplot2::aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = celltype
    ),
    arrow = ggplot2::arrow(
      length = ggplot2::unit(arrow.size * 0.1, "inches"),
      type = arrow.type
    ),
    alpha = arrow.alpha
  )

  # For faceted plots
  if (use_facets) {
    p <- p + ggplot2::facet_wrap(~ comparison, ncol = facet.ncol)
  }

  # Add quadrant lines if requested
  if (add.quadrants) {
    # Find midpoints for reference lines
    x_mid <- mean(range(arrows_data$x, arrows_data$xend, na.rm = TRUE))
    y_mid <- mean(range(arrows_data$y, arrows_data$yend, na.rm = TRUE))

    p <- p +
      # Add reference lines
      ggplot2::geom_hline(
        yintercept = y_mid,
        linetype = quadrant.line.type,
        color = quadrant.line.color,
        alpha = quadrant.line.alpha
      ) +
      ggplot2::geom_vline(
        xintercept = x_mid,
        linetype = quadrant.line.type,
        color = quadrant.line.color,
        alpha = quadrant.line.alpha
      )

    # Add quadrant labels
    if (!use_facets) {
      # For non-faceted plots, add quadrant labels
      x_range <- range(c(arrows_data$x, arrows_data$xend), na.rm = TRUE)
      y_range <- range(c(arrows_data$y, arrows_data$yend), na.rm = TRUE)

      p <- p +
        ggplot2::annotate(
          "text",
          x = max(x_range) * 0.95,
          y = max(y_range) * 0.95,
          label = "High sender\nHigh receiver",
          size = quadrant.label.size,
          alpha = quadrant.label.alpha
        ) +
        ggplot2::annotate(
          "text",
          x = min(x_range) * 0.95,
          y = max(y_range) * 0.95,
          label = "Low sender\nHigh receiver",
          size = quadrant.label.size,
          alpha = quadrant.label.alpha
        ) +
        ggplot2::annotate(
          "text",
          x = min(x_range) * 0.95,
          y = min(y_range) * 0.95,
          label = "Low sender\nLow receiver",
          size = quadrant.label.size,
          alpha = quadrant.label.alpha
        ) +
        ggplot2::annotate(
          "text",
          x = max(x_range) * 0.95,
          y = min(y_range) * 0.95,
          label = "High sender\nLow receiver",
          size = quadrant.label.size,
          alpha = quadrant.label.alpha
        )
    }
  }

  # Add group labels if requested
  if (show.group.labels) {
    # Aggregate data by group to find centroid positions
    if (requireNamespace("dplyr", quietly = TRUE)) {
      if (use_facets) {
        # For faceted plots, calculate centroids within each facet
        group_centroids <- arrows_data %>%
          dplyr::group_by(comparison, condition2) %>%
          dplyr::summarize(
            x = mean(xend, na.rm = TRUE),
            y = mean(yend, na.rm = TRUE)
          )
      } else {
        # For non-faceted plots, calculate overall centroids
        group_centroids <- arrows_data %>%
          dplyr::group_by(condition2) %>%
          dplyr::summarize(
            x = mean(xend, na.rm = TRUE),
            y = mean(yend, na.rm = TRUE)
          )
      }
    } else {
      # Fallback if dplyr is not available
      if (use_facets) {
        # For faceted plots
        facets <- unique(arrows_data$comparison)
        group_centroids <- data.frame()
        for (facet in facets) {
          facet_data <- arrows_data[arrows_data$comparison == facet, ]
          groups <- unique(facet_data$condition2)
          for (grp in groups) {
            grp_data <- facet_data[facet_data$condition2 == grp, ]
            group_centroids <- rbind(group_centroids,
                                     data.frame(
                                       comparison = facet,
                                       condition2 = grp,
                                       x = mean(grp_data$xend, na.rm = TRUE),
                                       y = mean(grp_data$yend, na.rm = TRUE)
                                     ))
          }
        }
      } else {
        # For non-faceted plots
        groups <- unique(arrows_data$condition2)
        group_centroids <- data.frame()
        for (grp in groups) {
          grp_data <- arrows_data[arrows_data$condition2 == grp, ]
          group_centroids <- rbind(group_centroids,
                                   data.frame(
                                     condition2 = grp,
                                     x = mean(grp_data$xend, na.rm = TRUE),
                                     y = mean(grp_data$yend, na.rm = TRUE)
                                   ))
        }
      }
    }

    # Add group labels with black text
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_label_repel(
        data = group_centroids,
        ggplot2::aes(
          x = x,
          y = y,
          label = condition2
        ),
        size = group.label.size,
        fontface = "bold",
        alpha = 0.9,
        box.padding = 1,
        point.padding = 0.5,
        min.segment.length = 0,
        color = "black",
        fill = "white",
        segment.color = "gray50"
      )
    } else {
      p <- p + ggplot2::geom_label(
        data = group_centroids,
        ggplot2::aes(
          x = x,
          y = y,
          label = condition2
        ),
        size = group.label.size,
        fontface = "bold",
        alpha = 0.9,
        color = "black",
        fill = "white"
      )
    }
  }

  # Add text labels for cell types if requested
  if (label.cell) {
    # Group by celltype and get the entry with max total_change for each
    if (requireNamespace("dplyr", quietly = TRUE)) {
      label_data <- arrows_data %>%
        dplyr::group_by(celltype) %>%
        dplyr::arrange(dplyr::desc(total_change)) %>%
        dplyr::slice(1)

      if (use_facets) {
        # For faceted plots, ensure one label per celltype per facet
        label_data <- arrows_data %>%
          dplyr::group_by(comparison, celltype) %>%
          dplyr::arrange(dplyr::desc(total_change)) %>%
          dplyr::slice(1)
      }
    } else {
      # Fallback if dplyr is not available
      cells <- unique(arrows_data$celltype)
      label_data <- data.frame()

      if (use_facets) {
        # For faceted plots
        facets <- unique(arrows_data$comparison)
        for (facet in facets) {
          facet_data <- arrows_data[arrows_data$comparison == facet, ]
          for (cell in cells) {
            if (cell %in% facet_data$celltype) {
              cell_data <- facet_data[facet_data$celltype == cell, ]
              max_row <- which.max(cell_data$total_change)
              label_data <- rbind(label_data, cell_data[max_row, ])
            }
          }
        }
      } else {
        # For non-faceted plots
        for (cell in cells) {
          cell_data <- arrows_data[arrows_data$celltype == cell, ]
          max_row <- which.max(cell_data$total_change)
          label_data <- rbind(label_data, cell_data[max_row, ])
        }
      }
    }

    if (label.repel && requireNamespace("ggrepel", quietly = TRUE)) {
      # Use ggrepel for smart label placement
      if (label.box) {
        # Use label_repel with box
        p <- p + ggrepel::geom_label_repel(
          data = label_data,
          ggplot2::aes(
            x = xend,
            y = yend,
            label = celltype,
            color = celltype
          ),
          size = label.size,
          box.padding = ggplot2::unit(0.35, "lines"),
          point.padding = ggplot2::unit(0.3, "lines"),
          min.segment.length = label.min.segment.length,
          max.overlaps = label.max.overlaps,
          segment.color = "grey50",
          segment.alpha = 0.6,
          seed = 42  # For reproducibility
        )
      } else {
        # Use text_repel without box
        p <- p + ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(
            x = xend,
            y = yend,
            label = celltype,
            color = celltype
          ),
          size = label.size,
          box.padding = ggplot2::unit(0.35, "lines"),
          point.padding = ggplot2::unit(0.3, "lines"),
          min.segment.length = label.min.segment.length,
          max.overlaps = label.max.overlaps,
          segment.color = "grey50",
          segment.alpha = 0.6,
          seed = 42  # For reproducibility
        )
      }
    } else {
      # Use regular text labels
      p <- p + ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(
          x = xend,
          y = yend,
          label = celltype,
          color = celltype
        ),
        hjust = -0.2,
        vjust = 0.5,
        size = label.size
      )
    }
  }

  # Apply scales - color for cell types, shape for groups
  p <- p +
    ggplot2::scale_color_manual(values = colors, name = "Cell Type") +
    ggplot2::scale_shape_manual(values = group.shapes, name = legend_title)

  # Add size scale
  p <- p + ggplot2::scale_size_continuous(range = point.size.range)

  # Continue with the rest of the plot
  p <- p +
    # Labels
    ggplot2::labs(
      x = ifelse(log.scale, "Sender Role (log scale)", "Sender Role"),
      y = ifelse(log.scale, "Receiver Role (log scale)", "Receiver Role"),
      title = paste0(title, "\n(", pathways.used, " pathways)"),
      size = "LR Count"
    ) +
    base.theme +
    ggplot2::theme(
      legend.position = legend.position,
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = title.size,
        face = title.face
      ),
      axis.text = ggplot2::element_text(size = axis.text.size),
      axis.title = ggplot2::element_text(size = axis.title.size),
      strip.background = ggplot2::element_rect(fill = "lightgrey"),
      strip.text = ggplot2::element_text(size = 10, face = "bold")
    )

  return(p)
}
