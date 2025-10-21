#' Compare signaling patterns across multiple conditions with a heatmap visualization
#'
#' This function creates a differential heatmap visualization comparing signaling patterns
#' across multiple conditions. It calculates normalized signaling strength for each cell type
#' and pathway, then visualizes the differences between conditions.
#'
#' @param object.list A list of CellChat objects to compare.
#' @param comparison A numeric vector specifying the indices of the objects to compare.
#'   If NULL (default), all objects will be compared.
#' @param reference Integer or character, index or name of the reference object to compare against.
#'   If NULL (default), the first object in comparison will be used as reference.
#' @param measure Character, specifies the signaling role to analyze: "sender" (outgoing signaling),
#'   "receiver" (incoming signaling), or "both" (overall signaling).
#' @param signaling Character vector of pathway names to include in the analysis. If NULL (default),
#'   all pathways present in any condition will be used.
#' @param color.heatmap Either a string specifying a predefined color scheme ("custom", "RdBu", "BuOr")
#'   or a vector of colors to create a custom gradient. Default is "custom" (blue-white-red).
#' @param color.use Vector of colors for cell type annotations. If NULL, colors will be generated automatically.
#' @param title Character, custom title for the heatmap. If NULL, a default title will be generated.
#' @param width Numeric, width of the heatmap in centimeters.
#' @param height Numeric, height of the heatmap in centimeters.
#' @param font.size Numeric, font size for row and column labels.
#' @param font.size.title Numeric, font size for the heatmap title.
#' @param cluster_rows Logical, whether to cluster rows (cell types in default orientation).
#' @param cluster_cols Logical, whether to cluster columns (pathways in default orientation).
#' @param show_row_dend Logical, whether to display the row dendrogram when clustering.
#' @param show_column_dend Logical, whether to display the column dendrogram when clustering.
#' @param show_rownames Logical, whether to display row names.
#' @param show_colnames Logical, whether to display column names.
#' @param show_values Logical, whether to display the difference values in cells.
#' @param value_digits Integer, number of decimal places to display for cell values.
#' @param signif_mark Logical, whether to mark significant differences.
#' @param border_color Character, color for cell borders.
#' @param show_heatmap_border Logical, whether to display cell borders.
#' @param border_width Numeric, width of cell borders.
#' @param use_log2fc Logical, whether to use log2 fold change instead of simple difference.
#' @param transpose Logical, whether to transpose the heatmap (pathways as rows, cell types as columns).
#' @param return_data Logical, whether to return a list with the heatmap and data matrices (TRUE)
#'   or just the heatmap object (FALSE).
#' @param save_plot Logical, whether to save the heatmap to a file.
#' @param save_name Character, filename for saving the plot (if save_plot=TRUE).
#' @param filter_zeros Logical, whether to filter out pathways with no signal.
#' @param filter_thresh Numeric, threshold for filtering pathways based on signal strength.
#' @param slot.name Character, name of the slot to extract data from (default: "netP").
#' @param cluster_conditions Logical, whether to cluster conditions (default: FALSE).
#' @param split_heatmap Logical, whether to split the heatmap by measure (useful when comparing multiple measures).
#' @param highlight_reference Logical, whether to highlight the reference condition (default: TRUE).
#' @param show_condition_names Logical, whether to display condition names (default: TRUE).
#' @param show_annotation Character, type of annotation to show: "both", "row", "column", or "none".
#' @param comparison_method Character, method for comparing conditions: "all_vs_ref" (all vs. reference),
#'   "all_vs_all" (pairwise comparisons between all conditions), or "custom_pairs" (specified pairs).
#' @param custom_comparisons List of custom comparison pairs (required when comparison_method = "custom_pairs").
#'   Each element should be a vector of length 2 containing either indices or names of conditions to compare.
#' @param big_heatmap Logical, whether to create a single big heatmap with side-by-side comparisons (default: FALSE).
#'
#' @return If return_data=TRUE, returns a list containing:
#'   \itemize{
#'     \item heatmap: The ComplexHeatmap object
#'     \item data_list: List of matrices for each condition
#'     \item normalized_list: List of normalized matrices for each condition
#'     \item diff_matrix: Difference matrix used for visualization (normalized)
#'     \item comparison: The comparison indices used
#'     \item reference: The reference index used
#'     \item condition_names: Names of the conditions compared
#'     \item measure: The measure type used
#'     \item pathways: Names of pathways included in the analysis
#'     \item cell_types: Names of cell types included in the analysis
#'   }
#'   If return_data=FALSE, returns only the ComplexHeatmap object.
#'
#' @examples
#' # Basic usage with default parameters
#' heatDiffM(cellchat.list, measure = "sender")
#'
#' # Advanced usage with customization
#' heatDiffM(cellchat.list,
#'          measure = "both",
#'          reference = "Normal",
#'          color.heatmap = c("purple", "white", "green"),
#'          use_log2fc = TRUE,
#'          cluster_rows = TRUE,
#'          cluster_cols = TRUE,
#'          show_values = TRUE,
#'          filter_thresh = 0.1)
#'
#' # Big side-by-side heatmap
#' heatDiffM(cellchat.list,
#'          measure = "sender",
#'          big_heatmap = TRUE)
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation anno_barplot
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#'
#' @export
heatDiffM <- function(object.list, comparison = NULL, reference = NULL,
                      measure = c("both", "sender", "receiver"), signaling = NULL,
                      color.heatmap = "custom", color.use = NULL,
                      title = NULL, width = 10, height = 8,
                      font.size = 8, font.size.title = 10,
                      cluster_rows = FALSE, cluster_cols = FALSE,
                      show_row_dend = TRUE, show_column_dend = TRUE,
                      show_rownames = TRUE, show_colnames = TRUE,
                      show_values = FALSE, value_digits = 2,
                      signif_mark = FALSE, border_color = "white",
                      show_heatmap_border = TRUE, border_width = 1,
                      use_log2fc = FALSE, transpose = TRUE,
                      return_data = FALSE, save_plot = FALSE,
                      save_name = "heatmap_diff.pdf",
                      filter_zeros = FALSE, filter_thresh = 0,
                      slot.name = "netP", cluster_conditions = FALSE,
                      split_heatmap = FALSE, highlight_reference = TRUE,
                      show_condition_names = TRUE,
                      show_annotation = c("both", "row", "column", "none"),
                      comparison_method = c("all_vs_ref", "all_vs_all", "custom_pairs"),
                      custom_comparisons = NULL,
                      big_heatmap = FALSE) {

  # Match arguments
  measure <- match.arg(measure)
  show_annotation <- match.arg(show_annotation)
  comparison_method <- match.arg(comparison_method)

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

  # Validate inputs
  if (length(comparison) < 2) {
    stop("At least two conditions are required for comparison")
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

  # Extract cell types (use only cell types present in all conditions)
  cell.types <- Reduce(
    intersect,
    lapply(comparison, function(idx) rownames(object.list[[idx]]@net$weight))
  )

  if (length(cell.types) == 0) {
    stop("No shared cell types found across all conditions")
  }

  # Extract ALL pathways from ALL conditions if not provided
  if (is.null(signaling)) {
    all_pathways <- unique(unlist(
      lapply(comparison, function(idx) unlist(slot(object.list[[idx]], slot.name)$pathways))
    ))

    signaling <- all_pathways
    message(paste("Using all", length(signaling), "pathways from all conditions"))
  }

  # Initialize data structures
  data_list <- list()
  data_norm_list <- list()

  # Extract data for each condition
  for (i in 1:length(comparison)) {
    idx <- comparison[i]
    obj <- object.list[[idx]]

    # Initialize matrix
    mat <- matrix(0, nrow = length(cell.types), ncol = length(signaling))
    rownames(mat) <- cell.types
    colnames(mat) <- signaling

    # Check if centrality is computed
    if (length(slot(obj, slot.name)$centr) == 0) {
      warning(paste0("Centrality metrics not computed for condition ", i,
                     ". Please run `netAnalysis_computeCentrality` first."))
      data_list[[i]] <- mat
      next
    }

    # Get centrality metrics
    centr <- slot(obj, slot.name)$centr

    # For each pathway, get centrality metrics
    for (p in 1:length(signaling)) {
      pathway <- signaling[p]

      if (pathway %in% names(centr)) {
        # Use centrality metrics directly
        for (cell in cell.types) {
          if (measure == "sender") {
            # Outgoing signaling
            if (cell %in% names(centr[[pathway]]$outdeg)) {
              mat[cell, p] <- centr[[pathway]]$outdeg[cell]
            }
          } else if (measure == "receiver") {
            # Incoming signaling
            if (cell %in% names(centr[[pathway]]$indeg)) {
              mat[cell, p] <- centr[[pathway]]$indeg[cell]
            }
          } else if (measure == "both") {
            # Both incoming and outgoing
            sum_val <- 0
            if (cell %in% names(centr[[pathway]]$outdeg)) {
              sum_val <- sum_val + centr[[pathway]]$outdeg[cell]
            }
            if (cell %in% names(centr[[pathway]]$indeg)) {
              sum_val <- sum_val + centr[[pathway]]$indeg[cell]
            }
            mat[cell, p] <- sum_val
          }
        }
      }
    }

    # Store the raw data matrix
    data_list[[i]] <- mat

    # Create normalized matrix for visualization
    mat_norm <- sweep(mat, 2L, apply(mat, 2, max), "/", check.margin = FALSE)
    mat_norm[is.na(mat_norm)] <- 0
    data_norm_list[[i]] <- mat_norm
  }

  # Filter out pathways with low signal if requested
  if (filter_zeros || filter_thresh > 0) {
    # Calculate maximum signal for each pathway across all conditions
    pathway_max_signal <- numeric(length(signaling))
    names(pathway_max_signal) <- signaling

    for (i in 1:length(data_list)) {
      for (p in 1:length(signaling)) {
        pathway <- signaling[p]
        pathway_max_signal[pathway] <- max(pathway_max_signal[pathway],
                                           sum(data_list[[i]][, pathway]),
                                           na.rm = TRUE)
      }
    }

    # Find pathways to keep
    pathways_to_keep <- which(pathway_max_signal > filter_thresh)

    if (length(pathways_to_keep) == 0) {
      message("No pathways pass the filtering threshold. Keeping all pathways.")
    } else if (length(pathways_to_keep) < length(signaling)) {
      message(paste("Filtering out", length(signaling) - length(pathways_to_keep),
                    "pathways with low signal."))
      # Filter matrices
      signaling <- signaling[pathways_to_keep]
      for (i in 1:length(data_list)) {
        data_list[[i]] <- data_list[[i]][, pathways_to_keep, drop = FALSE]
        data_norm_list[[i]] <- data_norm_list[[i]][, pathways_to_keep, drop = FALSE]
      }
    }
  }

  # Calculate difference matrices based on comparison method
  diff_matrices <- list()
  comparison_labels <- list()

  if (comparison_method == "all_vs_ref") {
    # Get reference matrices
    ref_idx <- which(comparison == reference)
    ref_mat <- data_list[[ref_idx]]
    ref_mat_norm <- data_norm_list[[ref_idx]]

    # Calculate differences (each condition vs. reference)
    # Create a counter for non-reference conditions
    diff_count <- 0

    for (i in 1:length(comparison)) {
      if (i != ref_idx) {  # Skip the reference condition itself
        diff_count <- diff_count + 1

        if (use_log2fc) {
          # Calculate log2 fold change with pseudocount
          pseudo_count <- 0.01
          diff_mat <- log2((data_norm_list[[i]] + pseudo_count) /
                             (ref_mat_norm + pseudo_count))
          # Cap extreme values
          cap_value <- 3
          diff_mat[diff_mat > cap_value] <- cap_value
          diff_mat[diff_mat < -cap_value] <- -cap_value
          diff_mat[is.infinite(diff_mat)] <- NA
        } else {
          # Simple difference
          diff_mat <- data_norm_list[[i]] - ref_mat_norm
        }

        # Store the difference matrix
        diff_matrices[[diff_count]] <- diff_mat

        # Create a label for this comparison
        comparison_labels[[diff_count]] <- paste(condition_names[comparison[i]],
                                                 "vs.",
                                                 condition_names[reference])
      }
    }
  } else if (comparison_method == "all_vs_all") {
    # Calculate all pairwise differences
    diff_count <- 0
    for (i in 1:(length(comparison)-1)) {
      for (j in (i+1):length(comparison)) {
        diff_count <- diff_count + 1

        if (use_log2fc) {
          # Calculate log2 fold change with pseudocount
          pseudo_count <- 0.01
          diff_mat <- log2((data_norm_list[[j]] + pseudo_count) /
                             (data_norm_list[[i]] + pseudo_count))
          # Cap extreme values
          cap_value <- 3
          diff_mat[diff_mat > cap_value] <- cap_value
          diff_mat[diff_mat < -cap_value] <- -cap_value
          diff_mat[is.infinite(diff_mat)] <- NA
        } else {
          # Simple difference
          diff_mat <- data_norm_list[[j]] - data_norm_list[[i]]
        }

        diff_matrices[[diff_count]] <- diff_mat
        comparison_labels[[diff_count]] <- paste(condition_names[comparison[j]],
                                                 "vs.",
                                                 condition_names[comparison[i]])
      }
    }
  } else if (comparison_method == "custom_pairs") {
    # Use custom comparison pairs
    if (is.null(custom_comparisons) || !is.list(custom_comparisons)) {
      stop("For custom_pairs comparison method, please provide a list of comparison pairs")
    }

    # Validate and process custom comparisons
    comparison_pairs <- list()
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

      # Convert object indices to positions within comparison vector
      pos1 <- which(comparison == pair[1])
      pos2 <- which(comparison == pair[2])
      comparison_pairs[[i]] <- c(pos1, pos2)
    }

    # Calculate differences for custom pairs
    for (p in 1:length(comparison_pairs)) {
      pair <- comparison_pairs[[p]]
      # pair now contains positions in the comparison vector
      idx1 <- pair[1]
      idx2 <- pair[2]

      if (use_log2fc) {
        # Calculate log2 fold change with pseudocount
        pseudo_count <- 0.01
        diff_mat <- log2((data_norm_list[[idx2]] + pseudo_count) /
                           (data_norm_list[[idx1]] + pseudo_count))
        # Cap extreme values
        cap_value <- 3
        diff_mat[diff_mat > cap_value] <- cap_value
        diff_mat[diff_mat < -cap_value] <- -cap_value
        diff_mat[is.infinite(diff_mat)] <- NA
      } else {
        # Simple difference
        diff_mat <- data_norm_list[[idx2]] - data_norm_list[[idx1]]
      }

      diff_matrices[[p]] <- diff_mat
      comparison_labels[[p]] <- paste(condition_names[comparison[idx2]],
                                      "vs.",
                                      condition_names[comparison[idx1]])
    }
  }

  # Check if any difference matrices were created
  if (length(diff_matrices) == 0) {
    stop("No viable comparisons found. Check your reference and comparison settings.")
  }

  # Set colors for the heatmap
  if (is.character(color.heatmap) && length(color.heatmap) == 1) {
    if (color.heatmap == "custom") {
      # Custom color palette: blue-white-red
      color.heatmap.use <- colorRampPalette(c("blue", "white", "red"))(100)
    } else if (color.heatmap == "BuOr") {
      # Blue-Orange palette
      color.heatmap.use <- colorRampPalette(c("blue", "white", "orange"))(100)
    } else if (color.heatmap %in% rownames(RColorBrewer::brewer.pal.info)) {
      # Default RColorBrewer palette
      color.heatmap.use <- (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9,
                                                                                  name = color.heatmap))))(100)
    } else {
      # Fallback to default blue-white-red
      color.heatmap.use <- colorRampPalette(c("blue", "white", "red"))(100)
    }
  } else if (is.character(color.heatmap) && length(color.heatmap) >= 2) {
    # User provided a vector of colors
    color.heatmap.use <- colorRampPalette(color.heatmap)(100)
  } else {
    # Fallback to default blue-white-red
    color.heatmap.use <- colorRampPalette(c("blue", "white", "red"))(100)
  }

  # Set cell type colors for annotation
  if (is.null(color.use)) {
    if (exists("global_colors")) {
      color.use <- global_colors[1:length(cell.types)]
    } else {
      color.use <- rainbow(length(cell.types))
    }
    names(color.use) <- cell.types
  } else {
    # Ensure color.use has names even if user provided it
    if (is.null(names(color.use))) {
      if (length(color.use) >= length(cell.types)) {
        names(color.use) <- cell.types
      } else {
        warning("color.use has fewer colors than cell types. Extending with rainbow colors.")
        color.use <- c(color.use, rainbow(length(cell.types) - length(color.use)))
        names(color.use) <- cell.types
      }
    }
  }

  # Create title with measure description
  measure_description <- switch(
    measure,
    "sender" = "Outgoing",
    "receiver" = "Incoming",
    "both" = "Overall"
  )

  if (is.null(title)) {
    if (comparison_method == "all_vs_ref") {
      title <- paste0("Differential ", measure_description, " signaling patterns vs. ",
                      condition_names[reference])
    } else {
      title <- paste0("Pairwise differential ", measure_description, " signaling patterns")
    }

    if (use_log2fc) {
      title <- paste0(title, " (log2FC)")
    }
  }

  # Create a big side-by-side heatmap if requested
  if (big_heatmap && length(diff_matrices) > 1) {
    # For big single heatmap with side-by-side comparisons

    # Create a heatmap list to store individual heatmaps
    heatmap_list <- list()

    # Calculate a global legend break for consistent colors across all heatmaps
    all_values <- unlist(lapply(diff_matrices, function(x) as.vector(x)))
    all_values <- all_values[!is.na(all_values)]  # Remove NA values

    # Check if we have valid values
    if (length(all_values) == 0) {
      warning("All difference values are NA. Using default legend breaks.")
      legend_break <- c(-1, 0, 1)
    } else {
      min_val <- min(all_values, na.rm = TRUE)
      max_val <- max(all_values, na.rm = TRUE)

      # Check for infinite values
      if (is.infinite(min_val) || is.infinite(max_val)) {
        warning("Infinite values detected. Using default legend breaks.")
        legend_break <- c(-1, 0, 1)
      } else if (min_val == max_val) {
        # All values are the same
        legend_break <- c(min_val - 0.5, min_val, min_val + 0.5)
      } else {
        # Create symmetric breaks around zero for divergent color scale
        if (min_val < 0 && max_val > 0) {
          abs_max <- max(abs(min_val), abs(max_val))
          legend_break <- c(-abs_max, 0, abs_max)
        } else {
          # If all values are positive or all negative
          legend_break <- c(min_val, (min_val + max_val)/2, max_val)
        }

        # Round to one decimal place
        legend_break <- round(legend_break, 1)

        # Ensure values are distinct after rounding
        if (length(unique(legend_break)) < 2) {
          legend_break <- c(min_val - 0.5, min_val, min_val + 0.5)
        }
      }
    }

    # Create the color mapping function (shared across all heatmaps)
    col_fun <- circlize::colorRamp2(legend_break, color.heatmap.use[c(1, 50, 100)])

    # Legend name (will only be used for the first heatmap)
    legend_name <- if(use_log2fc) "Log2FC" else "Difference"

    # Generate the heatmaps with specific named legends only for the first one
    for (i in 1:length(diff_matrices)) {
      if (transpose) {
        # For transposed version (pathways as rows, cell types as columns)
        t_diff_matrix <- t(diff_matrices[[i]])

        # Generate the heatmap
        # Critical fix: Use a space character name for non-first heatmaps
        if (i == 1) {
          # First heatmap gets a proper legend
          hm <- ComplexHeatmap::Heatmap(
            t_diff_matrix,
            col = col_fun,
            na_col = "white",
            name = legend_name,  # Use the legend name for first heatmap only
            cluster_rows = cluster_cols,
            cluster_columns = cluster_rows,
            show_row_dend = show_column_dend && cluster_cols,
            show_column_dend = show_row_dend && cluster_rows,
            show_row_names = show_colnames && i==1,  # Only show row names on the first heatmap
            show_column_names = show_rownames,
            row_names_side = "left",
            row_names_rot = 0,
            row_names_gp = grid::gpar(fontsize = font.size),
            column_names_gp = grid::gpar(fontsize = font.size),
            width = grid::unit(width/length(diff_matrices), "cm"),  # Adjust width
            height = grid::unit(height, "cm"),
            column_title = comparison_labels[[i]],  # Use comparison label as title
            column_title_gp = grid::gpar(fontsize = font.size),
            column_names_rot = 90,
            # Cell function to display values
            cell_fun = if(show_values) {
              function(j, i, x, y, width, height, fill) {
                value <- t_diff_matrix[i, j]
                if (!is.na(value) && abs(value) > 1e-10) {
                  display_val <- format(round(value, value_digits), nsmall = value_digits)
                  grid::grid.text(display_val, x, y,
                                  gp = grid::gpar(fontsize = font.size * 0.8))
                }
              }
            } else NULL,
            # Border settings
            border = if(show_heatmap_border) border_color else NA,
            border_gp = grid::gpar(lwd = border_width),
            rect_gp = grid::gpar(col = if(show_heatmap_border) border_color else NA,
                                 lwd = border_width),
            # Legend parameters
            heatmap_legend_param = list(
              title_gp = grid::gpar(fontsize = 8, fontface = "plain"),
              title_position = "leftcenter-rot",
              border = NA,
              at = legend_break,
              legend_height = grid::unit(20, "mm"),
              labels_gp = grid::gpar(fontsize = 8),
              grid_width = grid::unit(2, "mm")
            )
          )
        } else {
          # All other heatmaps get NO legend by using unique space-based names
          hm <- ComplexHeatmap::Heatmap(
            t_diff_matrix,
            col = col_fun,
            na_col = "white",
            name = paste0(rep(" ", i), collapse = ""),  # Unique name with multiple spaces
            show_heatmap_legend = FALSE,  # Explicitly don't show legend
            cluster_rows = cluster_cols,
            cluster_columns = cluster_rows,
            show_row_dend = show_column_dend && cluster_cols,
            show_column_dend = show_row_dend && cluster_rows,
            show_row_names = FALSE,  # Don't show row names on non-first heatmaps
            show_column_names = show_rownames,
            row_names_side = "left",
            row_names_rot = 0,
            row_names_gp = grid::gpar(fontsize = font.size),
            column_names_gp = grid::gpar(fontsize = font.size),
            width = grid::unit(width/length(diff_matrices), "cm"),  # Adjust width - FIXED!
            height = grid::unit(height, "cm"),  # FIXED!
            column_title = comparison_labels[[i]],  # Use comparison label as title
            column_title_gp = grid::gpar(fontsize = font.size),
            column_names_rot = 90,
            # Cell function to display values
            cell_fun = if(show_values) {
              function(j, i, x, y, width, height, fill) {
                value <- t_diff_matrix[i, j]
                if (!is.na(value) && abs(value) > 1e-10) {
                  display_val <- format(round(value, value_digits), nsmall = value_digits)
                  grid::grid.text(display_val, x, y,
                                  gp = grid::gpar(fontsize = font.size * 0.8))
                }
              }
            } else NULL,
            # Border settings
            border = if(show_heatmap_border) border_color else NA,
            border_gp = grid::gpar(lwd = border_width),
            rect_gp = grid::gpar(col = if(show_heatmap_border) border_color else NA,
                                 lwd = border_width)
          )
        }

        heatmap_list[[i]] <- hm

      } else {
        # For regular orientation (cell types as rows, pathways as columns)
        diff_matrix <- diff_matrices[[i]]

        # Generate the heatmap
        # Critical fix: Use a space character name for non-first heatmaps
        if (i == 1) {
          # First heatmap gets a proper legend
          hm <- ComplexHeatmap::Heatmap(
            diff_matrix,
            col = col_fun,
            na_col = "white",
            name = legend_name,  # Use the legend name for first heatmap only
            cluster_rows = cluster_rows,
            cluster_columns = cluster_cols,
            show_row_dend = show_row_dend && cluster_rows,
            show_column_dend = show_column_dend && cluster_cols,
            show_row_names = show_rownames && i==1,  # Only show row names on the first heatmap
            show_column_names = show_colnames,
            row_names_side = "left",
            row_names_rot = 0,
            row_names_gp = grid::gpar(fontsize = font.size),
            column_names_gp = grid::gpar(fontsize = font.size),
            width = grid::unit(width/length(diff_matrices), "cm"),  # Adjust width
            height = grid::unit(height, "cm"),
            column_title = comparison_labels[[i]],  # Use comparison label as title
            column_title_gp = grid::gpar(fontsize = font.size),
            column_names_rot = 90,
            # Cell function to display values
            cell_fun = if(show_values) {
              function(j, i, x, y, width, height, fill) {
                value <- diff_matrix[i, j]
                if (!is.na(value) && abs(value) > 1e-10) {
                  display_val <- format(round(value, value_digits), nsmall = value_digits)
                  grid::grid.text(display_val, x, y,
                                  gp = grid::gpar(fontsize = font.size * 0.8))
                }
              }
            } else NULL,
            # Border settings
            border = if(show_heatmap_border) border_color else NA,
            border_gp = grid::gpar(lwd = border_width),
            rect_gp = grid::gpar(col = if(show_heatmap_border) border_color else NA,
                                 lwd = border_width),
            # Legend parameters
            heatmap_legend_param = list(
              title_gp = grid::gpar(fontsize = 8, fontface = "plain"),
              title_position = "leftcenter-rot",
              border = NA,
              at = legend_break,
              legend_height = grid::unit(20, "mm"),
              labels_gp = grid::gpar(fontsize = 8),
              grid_width = grid::unit(2, "mm")
            )
          )
        } else {
          # All other heatmaps get NO legend by using unique space-based names
          hm <- ComplexHeatmap::Heatmap(
            diff_matrix,
            col = col_fun,
            na_col = "white",
            name = paste0(rep(" ", i), collapse = ""),  # Unique name with multiple spaces
            show_heatmap_legend = FALSE,  # Explicitly don't show legend
            cluster_rows = cluster_rows,
            cluster_columns = cluster_cols,
            show_row_dend = show_row_dend && cluster_rows,
            show_column_dend = show_column_dend && cluster_cols,
            show_row_names = FALSE,  # Don't show row names on non-first heatmaps
            show_column_names = show_colnames,
            row_names_side = "left",
            row_names_rot = 0,
            row_names_gp = grid::gpar(fontsize = font.size),
            column_names_gp = grid::gpar(fontsize = font.size),
            width = grid::unit(width/length(diff_matrices), "cm"),  # Adjust width
            height = grid::unit(height, "cm"),
            column_title = comparison_labels[[i]],  # Use comparison label as title
            column_title_gp = grid::gpar(fontsize = font.size),
            column_names_rot = 90,
            # Cell function to display values
            cell_fun = if(show_values) {
              function(j, i, x, y, width, height, fill) {
                value <- diff_matrix[i, j]
                if (!is.na(value) && abs(value) > 1e-10) {
                  display_val <- format(round(value, value_digits), nsmall = value_digits)
                  grid::grid.text(display_val, x, y,
                                  gp = grid::gpar(fontsize = font.size * 0.8))
                }
              }
            } else NULL,
            # Border settings
            border = if(show_heatmap_border) border_color else NA,
            border_gp = grid::gpar(lwd = border_width),
            rect_gp = grid::gpar(col = if(show_heatmap_border) border_color else NA,
                                 lwd = border_width)
          )
        }

        heatmap_list[[i]] <- hm
      }
    }

    # Create the correct HeatmapList using + operator
    if (length(heatmap_list) > 1) {
      # Combine heatmaps
      final_heatmap <- heatmap_list[[1]]
      for (i in 2:length(heatmap_list)) {
        final_heatmap <- final_heatmap + heatmap_list[[i]]
      }

      # Add a global title
      draw_params <- list(
        column_title = title,
        column_title_gp = grid::gpar(fontsize = font.size.title, fontface = "bold")
      )
    } else {
      final_heatmap <- heatmap_list[[1]]
      draw_params <- list()
    }

  } else {
    # Create individual heatmaps for standard display
    heatmap_list <- list()

    # Process each difference matrix
    for (i in 1:length(diff_matrices)) {
      diff_matrix <- diff_matrices[[i]]
      comp_name <- comparison_labels[[i]]

      # Create annotations
      if (show_annotation %in% c("both", "row") && !transpose) {
        row_anno <- ComplexHeatmap::rowAnnotation(
          Strength = ComplexHeatmap::anno_barplot(
            rowSums(abs(diff_matrix)),
            border = FALSE,
            gp = grid::gpar(fill = color.use[rownames(diff_matrix)])
          ),
          show_annotation_name = FALSE
        )
      } else {
        row_anno <- NULL
      }

      if (show_annotation %in% c("both", "column") && !transpose) {
        col_anno <- ComplexHeatmap::HeatmapAnnotation(
          Strength = ComplexHeatmap::anno_barplot(
            colSums(abs(diff_matrix)),
            border = FALSE
          ),
          show_annotation_name = FALSE
        )
      } else {
        col_anno <- NULL
      }

      # Calculate legend breaks
      min_val <- min(diff_matrix, na.rm = TRUE)
      max_val <- max(diff_matrix, na.rm = TRUE)

      # Check for infinite values (happens when all values are NA)
      if (is.infinite(min_val) || is.infinite(max_val)) {
        warning("All values in difference matrix are NA. Using default legend breaks.")
        legend.break <- c(-1, 0, 1)
      } else if (min_val == max_val) {
        # All values are the same
        legend.break <- c(min_val - 0.5, min_val, min_val + 0.5)
      } else {
        # Create symmetric breaks around zero for divergent color scale
        if (min_val < 0 && max_val > 0) {
          abs_max <- max(abs(min_val), abs(max_val))
          legend.break <- c(-abs_max, 0, abs_max)
        } else {
          # If all values are positive or all negative
          legend.break <- c(min_val, (min_val + max_val)/2, max_val)
        }

        # Round to one decimal place
        legend.break <- round(legend.break, 1)

        # Ensure values are distinct after rounding
        if (length(unique(legend.break)) < 2) {
          legend.break <- c(min_val - 0.5, min_val, min_val + 0.5)
        }
      }

      # Create the color mapping function
      col_fun <- circlize::colorRamp2(legend.break, color.heatmap.use[c(1, 50, 100)])

      # Determine heatmap subtitle
      hm_subtitle <- comp_name

      # Create the heatmap
      if (transpose) {
        # For transposed version (pathways as rows, cell types as columns)
        t_diff_matrix <- t(diff_matrix)

        # Create appropriate annotations for the transposed matrix
        if (show_annotation %in% c("both", "row")) {
          t_row_anno <- ComplexHeatmap::rowAnnotation(
            Strength = ComplexHeatmap::anno_barplot(
              rowSums(abs(t_diff_matrix)),
              border = FALSE
            ),
            show_annotation_name = FALSE
          )
        } else {
          t_row_anno <- NULL
        }

        if (show_annotation %in% c("both", "column")) {
          t_col_anno <- ComplexHeatmap::HeatmapAnnotation(
            Strength = ComplexHeatmap::anno_barplot(
              colSums(abs(t_diff_matrix)),
              border = FALSE,
              gp = grid::gpar(fill = color.use[colnames(t_diff_matrix)])
            ),
            show_annotation_name = FALSE
          )
        } else {
          t_col_anno <- NULL
        }

        heatmap_list[[i]] <- ComplexHeatmap::Heatmap(
          t_diff_matrix,
          col = col_fun,
          na_col = "white",
          name = if(use_log2fc) "Log2FC" else "Difference",
          right_annotation = t_row_anno,
          top_annotation = t_col_anno,
          cluster_rows = cluster_cols,  # Swap clustering settings
          cluster_columns = cluster_rows,
          show_row_dend = show_column_dend && cluster_cols,
          show_column_dend = show_row_dend && cluster_rows,
          show_row_names = show_colnames,
          show_column_names = show_rownames,
          row_names_side = "left",
          row_names_rot = 0,
          row_names_gp = grid::gpar(fontsize = font.size),
          column_names_gp = grid::gpar(fontsize = font.size),
          width = grid::unit(width, "cm"),
          height = grid::unit(height, "cm"),
          row_title = hm_subtitle,
          row_title_gp = grid::gpar(fontsize = font.size.title),
          column_names_rot = 90,
          row_names_max_width = grid::unit(10, "cm"),
          # Cell function to display values
          cell_fun = if(show_values) {
            function(j, i, x, y, width, height, fill) {
              value <- t_diff_matrix[i, j]
              if (!is.na(value) && abs(value) > 1e-10) {
                display_val <- format(round(value, value_digits), nsmall = value_digits)
                grid::grid.text(display_val, x, y,
                                gp = grid::gpar(fontsize = font.size * 0.8))
              }
            }
          } else NULL,
          # Border settings
          border = if(show_heatmap_border) border_color else NA,
          border_gp = grid::gpar(lwd = border_width),
          rect_gp = grid::gpar(col = if(show_heatmap_border) border_color else NA,
                               lwd = border_width),
          # Legend parameters
          heatmap_legend_param = list(
            title_gp = grid::gpar(fontsize = 8, fontface = "plain"),
            title_position = "leftcenter-rot",
            border = NA,
            at = legend.break,
            legend_height = grid::unit(20, "mm"),
            labels_gp = grid::gpar(fontsize = 8),
            grid_width = grid::unit(2, "mm")
          )
        )
      } else {
        # Original orientation (cell types as rows, pathways as columns)
        heatmap_list[[i]] <- ComplexHeatmap::Heatmap(
          diff_matrix,
          col = col_fun,
          na_col = "white",
          name = if(use_log2fc) "Log2FC" else "Difference",
          right_annotation = row_anno,
          top_annotation = col_anno,
          cluster_rows = cluster_rows,
          cluster_columns = cluster_cols,
          show_row_dend = show_row_dend && cluster_rows,
          show_column_dend = show_column_dend && cluster_cols,
          show_row_names = show_rownames,
          show_column_names = show_colnames,
          row_names_side = "left",
          row_names_rot = 0,
          row_names_gp = grid::gpar(fontsize = font.size),
          column_names_gp = grid::gpar(fontsize = font.size),
          width = grid::unit(width, "cm"),
          height = grid::unit(height, "cm"),
          column_title = hm_subtitle,
          column_title_gp = grid::gpar(fontsize = font.size.title),
          column_names_rot = 90,
          row_names_max_width = grid::unit(10, "cm"),
          # Cell function to display values
          cell_fun = if(show_values) {
            function(j, i, x, y, width, height, fill) {
              value <- diff_matrix[i, j]
              if (!is.na(value) && abs(value) > 1e-10) {
                display_val <- format(round(value, value_digits), nsmall = value_digits)
                grid::grid.text(display_val, x, y,
                                gp = grid::gpar(fontsize = font.size * 0.8))
              }
            }
          } else NULL,
          # Border settings
          border = if(show_heatmap_border) border_color else NA,
          border_gp = grid::gpar(lwd = border_width),
          rect_gp = grid::gpar(col = if(show_heatmap_border) border_color else NA,
                               lwd = border_width),
          # Legend parameters
          heatmap_legend_param = list(
            title_gp = grid::gpar(fontsize = 8, fontface = "plain"),
            title_position = "leftcenter-rot",
            border = NA,
            at = legend.break,
            legend_height = grid::unit(20, "mm"),
            labels_gp = grid::gpar(fontsize = 8),
            grid_width = grid::unit(2, "mm")
          )
        )
      }
    }

    # Combine heatmaps
    final_heatmap <- NULL
    draw_params <- list()

    if (length(heatmap_list) == 1) {
      final_heatmap <- heatmap_list[[1]]
    } else if (length(heatmap_list) > 1) {
      # Try to arrange in a grid with shared legends
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        # Convert heatmaps to grobs
        hm_grobs <- lapply(heatmap_list, function(hm) {
          grid::grid.grabExpr(ComplexHeatmap::draw(hm, newpage = FALSE))
        })

        # Arrange in a grid
        final_grid <- do.call(
          gridExtra::grid.arrange,
          c(hm_grobs,
            list(ncol = length(heatmap_list),
                 top = grid::textGrob(title, gp = grid::gpar(fontsize = font.size.title,
                                                             fontface = "bold"))))
        )

        # Convert back to a heatmap-like object
        final_heatmap <- grid::grid.grabExpr(print(final_grid))
      } else {
        # Fall back to using ComplexHeatmap's list approach
        if (length(heatmap_list) > 1) {
          final_heatmap <- heatmap_list[[1]]
          for (i in 2:length(heatmap_list)) {
            final_heatmap <- final_heatmap + heatmap_list[[i]]
          }

          # Prepare draw parameters with title
          draw_params <- list(
            column_title = title,
            column_title_gp = grid::gpar(fontsize = font.size.title, fontface = "bold")
          )
        } else {
          final_heatmap <- heatmap_list[[1]]
        }
      }
    }
  }

  # Save the plot if requested
  if (save_plot && !is.null(final_heatmap)) {
    if (!is.character(save_name)) {
      save_name <- "heatmap_diff.pdf"
    }

    tryCatch({
      pdf(save_name, width = width/2.54, height = height/2.54)
      if (length(draw_params) > 0) {
        ComplexHeatmap::draw(final_heatmap, column_title = draw_params$column_title,
                             column_title_gp = draw_params$column_title_gp)
      } else {
        ComplexHeatmap::draw(final_heatmap)
      }
      invisible(dev.off())
      message(paste("Heatmap saved to", save_name))
    }, error = function(e) {
      warning(paste("Failed to save heatmap:", e$message))
    })
  }

  # Return results based on user preference
  if (return_data) {
    # Return a list with all components
    result <- list(
      heatmap = final_heatmap,
      data_list = data_list,
      normalized_list = data_norm_list,
      diff_matrices = diff_matrices,
      comparison = comparison,
      reference = reference,
      condition_names = condition_names[comparison],
      measure = measure,
      pathways = signaling,
      cell_types = cell.types,
      comparison_method = comparison_method,
      comparison_labels = comparison_labels,
      draw_params = if(exists("draw_params") && length(draw_params) > 0) draw_params else NULL
    )

    return(result)
  } else {
    # Return just the heatmap - matches heatDiff pattern for correct RStudio display
    return(final_heatmap)
  }
}
