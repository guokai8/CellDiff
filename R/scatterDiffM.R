#' Multi-Group Scatter Plot of Pathway-Level Signaling Role Changes
#'
#' @description
#' Creates a scatter plot showing how cell types change their signaling roles
#' across multiple conditions, using pathway-level information from the netP slot.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (NULL for all objects)
#' @param reference Integer or character, index or name of the reference object to compare against.
#'   If NULL (default), the first object in comparison will be used as reference.
#' @param pathways Vector of pathway names to include (NULL for all)
#' @param measure.type "sender", "receiver", or "influence" (default)
#' @param cell.type.strategy Character, strategy for aligning cell types: "shared" (default) or "union"
#' @param arrow.size Size of trajectory arrows
#' @param arrow.alpha Transparency of trajectory arrows
#' @param thresh P-value threshold for significant interactions (default: 0.05)
#' @param title Plot title
#' @param label.cell Whether to label cell types
#' @param label.size Size of cell type labels (default: 3)
#' @param color.use Vector of colors for cell types
#' @param group.colors Vector of colors for groups/conditions
#' @param point.size Base size for points (default: 3)
#' @param point.alpha Base alpha for points (default: 0.8)
#' @param show.legend Whether to show the legend (default: TRUE)
#' @param legend.position Position of the legend (default: "right")
#' @param plot.type Type of plot: "combined" (all trajectories on one plot) or
#'   "facet" (separate plot for each condition vs reference) or
#'   "group_colored" (all trajectories on one plot with group coloring)
#' @param facet.ncol Number of columns when using faceted plot
#' @param common.scale Whether to use a common scale for all plots
#' @param highlight.cells Vector of cell types to highlight
#' @param highlight.color Color for highlighted cells
#' @param comparison_method Method for comparing groups: "all_vs_ref" (each group vs reference)
#'   or "sequential" (each group vs next) or "custom_pairs" (specified pairs)
#' @param custom_comparisons List of custom comparison pairs
#' @param show_group_labels Whether to show group labels near points (only for group_colored)
#' @param group_label_size Size of group labels (default: 3)
#' @param group_label_alpha Transparency of group labels (default: 0.9)
#' @param group_point_shape Shape to use for points by group (NULL for default shapes)
#' @param convex_hull Whether to draw convex hulls around points of same group
#' @param hull_alpha Transparency of convex hulls (default: 0.2)
#' @param save_plot Whether to save the plot to file
#' @param save_name Filename for saving the plot
#' @param save_width Width of the saved plot in inches
#' @param save_height Height of the saved plot in inches
#' @param return_data Whether to return the data along with the plot
#'
#' @return A ggplot object or a list with plot and data
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(cellchatlist)
#'
#' # Basic multi-condition scatter plot (shared cell types)
#' plot <- scatterDiffM(cellchatlist)
#'
#' # Use union strategy to include all cell types
#' plot <- scatterDiffM(cellchatlist,
#'                      cell.type.strategy = "union")
#'
#' # Compare specific conditions with custom reference
#' plot <- scatterDiffM(cellchatlist,
#'                      comparison = c(1, 2, 3),
#'                      reference = 1,
#'                      measure.type = "sender",
#'                      cell.type.strategy = "union")
#'
#' # Faceted plot with union strategy
#' plot <- scatterDiffM(cellchatlist,
#'                      measure.type = "influence",
#'                      cell.type.strategy = "union",
#'                      plot.type = "facet",
#'                      facet.ncol = 2)
#'
#' # Group-colored plot with convex hulls
#' plot <- scatterDiffM(cellchatlist,
#'                      cell.type.strategy = "union",
#'                      plot.type = "group_colored",
#'                      convex_hull = TRUE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment theme_bw labs scale_color_manual
#'   scale_size scale_alpha guides facet_wrap theme element_text ggtitle
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr arrange group_by summarize
#' @export
scatterDiffM <- function(object.list, comparison = NULL, reference = NULL,
                         pathways = NULL, measure.type = "influence",
                         cell.type.strategy = c("shared", "union"),
                         arrow.size = 1, arrow.alpha = 0.8,
                         thresh = 0.05, title = "Changes in Signaling Roles",
                         label.cell = TRUE, label.size = 3,
                         color.use = NULL, group.colors = NULL,
                         point.size = 3, point.alpha = 0.8,
                         show.legend = TRUE, legend.position = "right",
                         plot.type = c("combined", "facet", "group_colored"),
                         facet.ncol = 2, common.scale = TRUE,
                         highlight.cells = NULL, highlight.color = "red",
                         comparison_method = c("all_vs_ref", "sequential", "custom_pairs"),
                         custom_comparisons = NULL,
                         show_group_labels = FALSE, group_label_size = 3,
                         group_label_alpha = 0.9, group_point_shape = NULL,
                         convex_hull = FALSE, hull_alpha = 0.2,
                         save_plot = FALSE, save_name = "scatterDiffM.pdf",
                         save_width = 10, save_height = 8,
                         return_data = FALSE) {

  # Match arguments
  plot.type <- match.arg(plot.type)
  comparison_method <- match.arg(comparison_method)
  cell.type.strategy <- match.arg(cell.type.strategy)

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

  # Align cell types using the specified strategy
  alignment <- alignCellTypes(object.list, indices = comparison, strategy = cell.type.strategy)
  common_cell_types <- alignment$cell_types

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
  if (is.null(color.use)) {
    if (exists("global_colors") && length(global_colors) >= length(common_cell_types)) {
      color.use <- global_colors[1:length(common_cell_types)]
    } else {
      color.use <- rainbow(length(common_cell_types))
    }
    names(color.use) <- common_cell_types
  } else if (length(color.use) < length(common_cell_types)) {
    # Recycle colors if needed
    color.use <- rep(color.use, length.out = length(common_cell_types))
    names(color.use) <- common_cell_types
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
  if (is.null(group_point_shape) && plot.type == "group_colored") {
    # Common distinctive shapes
    shape_options <- c(16, 17, 15, 18, 8, 7, 9, 10, 11, 12, 13, 14)
    group_point_shape <- shape_options[1:length(comparison)]
    names(group_point_shape) <- condition_names[comparison]
  }

  # Initialize lists to store results
  all_scores <- list()
  all_arrows <- list()

  # Get signaling scores for all conditions
  for (i in 1:length(comparison)) {
    idx <- comparison[i]
    all_scores[[i]] <- calculateSignaling(
      object.list[[idx]],
      common_cell_types,
      pathways,
      measure.type
    )
  }

  # Create a combined data frame for scores
  scores_df <- data.frame(
    celltype = rep(common_cell_types, length(comparison)),
    score = unlist(all_scores),
    condition = rep(condition_names[comparison], each = length(common_cell_types))
  )

  # Find global min/max for consistent scales if requested
  if (common.scale) {
    score_range <- range(scores_df$score)
    x_lim <- score_range
    y_lim <- score_range

    # Add some padding
    padding <- 0.1 * diff(score_range)
    x_lim <- c(max(0, x_lim[1] - padding), x_lim[2] + padding)
    y_lim <- c(max(0, y_lim[1] - padding), y_lim[2] + padding)
  }

  # Create data for each comparison pair
  all_plot_data <- list()
  for (p in 1:length(comparison_pairs)) {
    pair <- comparison_pairs[[p]]
    idx1 <- comparison[pair[1]]
    idx2 <- comparison[pair[2]]

    # Get scores
    scores1 <- all_scores[[pair[1]]]
    scores2 <- all_scores[[pair[2]]]

    # Get interaction counts
    counts1_all <- getInteractionCounts(
      object.list[[idx1]],
      signaling = pathways,
      slot.name = "net",
      thresh = thresh
    )
    counts2_all <- getInteractionCounts(
      object.list[[idx2]],
      signaling = pathways,
      slot.name = "net",
      thresh = thresh
    )

    # Extract counts for common cell types only (set to 0 if not present)
    counts1 <- sapply(common_cell_types, function(ct) {
      if (ct %in% names(counts1_all)) counts1_all[ct] else 0
    })
    counts2 <- sapply(common_cell_types, function(ct) {
      if (ct %in% names(counts2_all)) counts2_all[ct] else 0
    })

    # Use the maximum count for each cell type
    counts <- pmax(counts1, counts2)

    # Create data frame for plotting
    plot_data <- data.frame(
      celltype = common_cell_types,
      condition1 = scores1,
      condition2 = scores2,
      condition1_name = condition_names[idx1],
      condition2_name = condition_names[idx2],
      condition1_idx = idx1,
      condition2_idx = idx2,
      count = counts
    )

    # Calculate changes
    plot_data$change <- plot_data$condition2 - plot_data$condition1
    plot_data$perc_change <- (plot_data$condition2 - plot_data$condition1) /
      (plot_data$condition1 + 1e-10) * 100
    plot_data$abs_perc_change <- abs(plot_data$perc_change)

    # Create significance categories based on percentage change
    plot_data$change_category <- cut(
      plot_data$abs_perc_change,
      breaks = c(-Inf, 25, 50, 75, 100, Inf),
      labels = c("<25%", "25-50%", "50-75%", "75-100%", ">100%")
    )

    # Add comparison identifier
    plot_data$comparison <- paste(condition_names[idx2], "vs", condition_names[idx1])
    plot_data$pair_id <- p

    # Store the data
    all_plot_data[[p]] <- plot_data
  }

  # Combine all plot data
  combined_data <- do.call(rbind, all_plot_data)

  # Highlight specific cells if requested
  if (!is.null(highlight.cells)) {
    combined_data$highlighted <- combined_data$celltype %in% highlight.cells
  } else {
    combined_data$highlighted <- FALSE
  }

  # Create legend title based on comparison method
  if (comparison_method == "all_vs_ref") {
    legend_title <- paste0("Group vs ", condition_names[reference])
  } else {
    legend_title <- "Group"
  }

  # Create base plot
  if (plot.type == "combined") {
    # All trajectories on one plot
    p <- ggplot2::ggplot(combined_data,
                         ggplot2::aes(x = condition1, y = condition2, color = celltype)) +
      # Identity line
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
      # Points with size based on count
      ggplot2::geom_point(
        ggplot2::aes(size = count, alpha = change_category),
        shape = 16
      ) +
      # Arrows to show direction of change
      ggplot2::geom_segment(
        ggplot2::aes(
          x = condition1,
          y = condition1,  # Starting point on the identity line
          xend = condition1,
          yend = condition2,
          color = celltype
        ),
        arrow = ggplot2::arrow(
          length = ggplot2::unit(arrow.size * 0.1, "inches"),
          type = "closed"
        ),
        alpha = arrow.alpha
      ) +
      # Labels
      ggplot2::labs(
        x = paste(measure.type, "score"),
        y = paste(measure.type, "score"),
        title = title,
        size = "Count",
        alpha = "% Change"
      ) +
      ggplot2::scale_color_manual(values = color.use) +
      ggplot2::scale_alpha_discrete(range = c(0.4, 1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = if (show.legend) legend.position else "none",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
      )

    # Add text labels if requested
    if (label.cell) {
      # Group by celltype and get max value for positioning
      if (requireNamespace("dplyr", quietly = TRUE)) {
        label_data <- combined_data %>%
          dplyr::group_by(celltype) %>%
          dplyr::arrange(dplyr::desc(abs(change))) %>%
          dplyr::slice(1)
      } else {
        # Fallback if dplyr is not available
        cells <- unique(combined_data$celltype)
        label_data <- data.frame()
        for (cell in cells) {
          cell_data <- combined_data[combined_data$celltype == cell, ]
          max_row <- which.max(abs(cell_data$change))
          label_data <- rbind(label_data, cell_data[max_row, ])
        }
      }

      # Add labels with ggrepel
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(label = celltype),
          size = label.size,
          box.padding = 0.5,
          point.padding = 0.2,
          min.segment.length = 0,
          max.overlaps = 20,
          color = ifelse(label_data$highlighted, highlight.color, "black")
        )
      } else {
        p <- p + ggplot2::geom_text(
          data = label_data,
          ggplot2::aes(label = celltype),
          size = label.size,
          hjust = -0.2,
          vjust = 0.5,
          color = ifelse(label_data$highlighted, highlight.color, "black")
        )
      }
    }

    # Apply common scale if requested
    if (common.scale) {
      p <- p +
        ggplot2::xlim(x_lim) +
        ggplot2::ylim(y_lim)
    }

  } else if (plot.type == "group_colored") {
    # Group-colored trajectories on one plot
    p <- ggplot2::ggplot(combined_data,
                         ggplot2::aes(x = condition1, y = condition2)) +
      # Identity line
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +

      # Add convex hulls if requested
      {if(convex_hull) ggplot2::geom_polygon(
        ggplot2::aes(group = condition2_name, fill = condition2_name),
        alpha = hull_alpha,
        data = function(df) {
          hull_data <- do.call(rbind, lapply(split(df, df$condition2_name), function(group) {
            ch <- grDevices::chull(group$condition1, group$condition2)
            group[ch, ]
          }))
          return(hull_data)
        }
      )} +

      # Points with coloring and shapes based on condition
      ggplot2::geom_point(
        ggplot2::aes(
          size = count,
          color = condition2_name,
          shape = condition2_name
        ),
        alpha = point.alpha
      ) +
      ggplot2::scale_shape_manual(values = group_point_shape) +

      # Arrows to show direction of change
      ggplot2::geom_segment(
        ggplot2::aes(
          x = condition1,
          y = condition1,  # Starting point on the identity line
          xend = condition1,
          yend = condition2,
          color = condition2_name
        ),
        arrow = ggplot2::arrow(
          length = ggplot2::unit(arrow.size * 0.1, "inches"),
          type = "closed"
        ),
        alpha = arrow.alpha
      ) +
      # Labels
      ggplot2::labs(
        x = paste(measure.type, "score", "-", condition_names[reference]),
        y = paste(measure.type, "score"),
        title = title,
        size = "Count",
        color = legend_title,
        fill = legend_title,
        shape = legend_title
      ) +
      ggplot2::scale_color_manual(values = group.colors) +
      ggplot2::scale_fill_manual(values = group.colors) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = if (show.legend) legend.position else "none",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
      )

    # Add cell type labels if requested
    if (label.cell) {
      # Add labels with ggrepel
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          ggplot2::aes(label = celltype),
          size = label.size,
          box.padding = 0.5,
          point.padding = 0.2,
          min.segment.length = 0,
          max.overlaps = 20,
          color = "black"
        )
      } else {
        p <- p + ggplot2::geom_text(
          ggplot2::aes(label = celltype),
          size = label.size,
          hjust = -0.2,
          vjust = 0.5,
          color = "black"
        )
      }
    }

    # Add group labels if requested
    if (show_group_labels) {
      # Aggregate data by group to find centroid positions
      if (requireNamespace("dplyr", quietly = TRUE)) {
        group_centroids <- combined_data %>%
          dplyr::group_by(condition2_name) %>%
          dplyr::summarize(
            x = mean(condition1),
            y = mean(condition2)
          )
      } else {
        # Fallback if dplyr is not available
        groups <- unique(combined_data$condition2_name)
        group_centroids <- data.frame()
        for (grp in groups) {
          grp_data <- combined_data[combined_data$condition2_name == grp, ]
          group_centroids <- rbind(group_centroids,
                                   data.frame(
                                     condition2_name = grp,
                                     x = mean(grp_data$condition1),
                                     y = mean(grp_data$condition2)
                                   ))
        }
      }

      # Add group labels
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_label_repel(
          data = group_centroids,
          ggplot2::aes(x = x, y = y, label = condition2_name, color = condition2_name),
          size = group_label_size,
          fontface = "bold",
          alpha = group_label_alpha,
          box.padding = 1,
          point.padding = 0.5,
          min.segment.length = 0
        )
      } else {
        p <- p + ggplot2::geom_label(
          data = group_centroids,
          ggplot2::aes(x = x, y = y, label = condition2_name, color = condition2_name),
          size = group_label_size,
          fontface = "bold",
          alpha = group_label_alpha
        )
      }
    }

    # Apply common scale if requested
    if (common.scale) {
      p <- p +
        ggplot2::xlim(x_lim) +
        ggplot2::ylim(y_lim)
    }

  } else if (plot.type == "facet") {
    # Separate plot for each condition vs reference
    p <- ggplot2::ggplot(combined_data,
                         ggplot2::aes(x = condition1, y = condition2, color = celltype)) +
      # Identity line
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
      # Points with size based on count
      ggplot2::geom_point(
        ggplot2::aes(size = count, alpha = change_category),
        shape = 16
      ) +
      # Arrows to show direction of change
      ggplot2::geom_segment(
        ggplot2::aes(
          x = condition1,
          y = condition1,  # Starting point on the identity line
          xend = condition1,
          yend = condition2,
          color = celltype
        ),
        arrow = ggplot2::arrow(
          length = ggplot2::unit(arrow.size * 0.1, "inches"),
          type = "closed"
        ),
        alpha = arrow.alpha
      ) +
      # Use facets to separate comparisons
      ggplot2::facet_wrap(~ comparison, ncol = facet.ncol) +
      # Labels
      ggplot2::labs(
        x = paste(measure.type, "score"),
        y = paste(measure.type, "score"),
        title = title,
        size = "Count",
        alpha = "% Change"
      ) +
      ggplot2::scale_color_manual(values = color.use) +
      ggplot2::scale_alpha_discrete(range = c(0.4, 1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = if (show.legend) legend.position else "none",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        strip.background = ggplot2::element_rect(fill = "lightgrey"),
        strip.text = ggplot2::element_text(size = 10, face = "bold")
      )

    # Add text labels if requested
    if (label.cell) {
      # Add labels with ggrepel
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          ggplot2::aes(label = celltype),
          size = label.size,
          box.padding = 0.5,
          point.padding = 0.2,
          min.segment.length = 0,
          max.overlaps = 10,
          color = ifelse(combined_data$highlighted, highlight.color, "black")
        )
      } else {
        p <- p + ggplot2::geom_text(
          ggplot2::aes(label = celltype),
          size = label.size,
          hjust = -0.2,
          vjust = 0.5,
          color = ifelse(combined_data$highlighted, highlight.color, "black")
        )
      }
    }

    # Apply common scale if requested
    if (common.scale) {
      p <- p +
        ggplot2::xlim(x_lim) +
        ggplot2::ylim(y_lim)
    }
  }

  # Save plot if requested
  if (save_plot) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      ggplot2::ggsave(
        filename = save_name,
        plot = p,
        width = save_width,
        height = save_height
      )
      message(paste("Plot saved to", save_name))
    } else {
      warning("ggplot2 is required for saving plots. Please install it.")
    }
  }

  # Return plot and data if requested
  if (return_data) {
    return(list(
      plot = p,
      data = combined_data,
      comparison_pairs = comparison_pairs
    ))
  } else {
    return(p)
  }
}

