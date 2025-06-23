#' Compare signaling patterns between two conditions with a heatmap visualization
#'
#' This function creates a differential heatmap visualization comparing signaling patterns
#' between two conditions. It calculates normalized signaling strength for each cell type
#' and pathway, then visualizes the differences between conditions.
#'
#' @param object.list A list of CellChat objects to compare.
#' @param comparison A numeric vector of length 2 specifying the indices of the objects to compare.
#' @param measure Character, specifies the signaling role to analyze: "sender" (outgoing signaling),
#'   "receiver" (incoming signaling), or "both" (overall signaling).
#' @param signaling Character vector of pathway names to include in the analysis. If NULL (default),
#'   all pathways present in either condition will be used.
#' @param color.heatmap Either a string specifying a predefined color scheme ("custom", "RdBu", "BuOr")
#'   or a vector of colors to create a custom gradient. Default is "custom" (deepskyblue-white-darkorange).
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
#' @param signif_mark Logical, whether to mark significant differences (currently not used).
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
#'
#' @return If return_data=TRUE (default), returns a list containing:
#'   \itemize{
#'     \item heatmap: The ComplexHeatmap object
#'     \item mat1: Raw data matrix for condition 1
#'     \item mat2: Raw data matrix for condition 2
#'     \item mat1_norm: Normalized data matrix for condition 1
#'     \item mat2_norm: Normalized data matrix for condition 2
#'     \item diff_matrix: Difference matrix used for visualization (normalized)
#'     \item diff_matrix_raw: Raw difference matrix (unnormalized)
#'     \item transposed: Logical indicating if the heatmap is transposed
#'     \item pathways: Names of pathways included in the analysis
#'     \item cell_types: Names of cell types included in the analysis
#'     \item measure: The measure type used
#'     \item comparison: The comparison indices used
#'     \item condition_names: Names of the conditions compared
#'   }
#'   If return_data=FALSE, returns only the ComplexHeatmap object.
#'
#' @examples
#' # Basic usage with default parameters
#' result <- heatDiff_all(cellchat.list, comparison = c(1, 2), measure = "sender")
#'
#' # Advanced usage with customization
#' result <- heatDiff_all(cellchat.list, comparison = c(1, 2),
#'                       measure = "both",
#'                       color.heatmap = c("purple", "white", "green"),
#'                       use_log2fc = TRUE,
#'                       cluster_rows = TRUE,
#'                       cluster_cols = TRUE,
#'                       show_values = TRUE,
#'                       filter_thresh = 0.1)
#'
#' # Just get the heatmap for plotting
#' heatmap <- heatDiff_all(cellchat.list, comparison = c(1, 2),
#'                        measure = "receiver", return_data = FALSE)
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation anno_barplot
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#'
#' @export
heatDiff <- function(object.list, comparison = c(1, 2),
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
                         slot.name = "netP") {

  # Validate inputs
  measure <- match.arg(measure)
  if (length(comparison) != 2) {
    stop("Please provide exactly 2 indices for comparison")
  }

  # Get condition names
  condition_names <- names(object.list)
  if (is.null(condition_names)) {
    condition_names <- paste0("Condition_", 1:length(object.list))
  }

  # Check for object list validity
  if (!is.list(object.list) || length(object.list) < max(comparison)) {
    stop("object.list must be a list containing at least as many elements as the maximum index in comparison")
  }

  # Extract cell types (use only cell types present in both conditions)
  cell.types <- intersect(
    rownames(object.list[[comparison[1]]]@net$weight),
    rownames(object.list[[comparison[2]]]@net$weight)
  )

  if (length(cell.types) == 0) {
    stop("No shared cell types found between conditions")
  }

  # Extract ALL pathways from BOTH conditions if not provided
  if (is.null(signaling)) {
    pathways1 <- unlist(slot(object.list[[comparison[1]]], slot.name)$pathways)
    pathways2 <- unlist(slot(object.list[[comparison[2]]], slot.name)$pathways)

    # Use UNION to get all pathways
    signaling <- union(pathways1, pathways2)

    message(paste("Using all", length(signaling), "pathways from both conditions"))
  }

  # Initialize matrices for both conditions
  mat1 <- matrix(0, nrow = length(cell.types), ncol = length(signaling))
  mat2 <- matrix(0, nrow = length(cell.types), ncol = length(signaling))

  rownames(mat1) <- rownames(mat2) <- cell.types
  colnames(mat1) <- colnames(mat2) <- signaling

  # Fill matrices for each condition
  for (i in 1:2) {
    obj <- object.list[[comparison[i]]]
    mat <- if (i == 1) mat1 else mat2

    # Check if centrality is computed
    if (length(slot(obj, slot.name)$centr) == 0) {
      stop(paste0("Please run `netAnalysis_computeCentrality` for condition ", i, " to compute the network centrality scores!"))
    }

    # Get centrality metrics
    centr <- slot(obj, slot.name)$centr

    # For each pathway, get centrality metrics
    for (p in 1:length(signaling)) {
      pathway <- signaling[p]

      if (pathway %in% names(centr)) {
        # Use centrality metrics directly
        for (cell in rownames(mat)) {
          if (cell %in% names(centr[[pathway]]$outdeg) || cell %in% names(centr[[pathway]]$indeg)) {
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
    }

    # Update the appropriate matrix
    if (i == 1) {
      mat1 <- mat
    } else {
      mat2 <- mat
    }
  }

  # Filter out pathways with low signal if requested
  if (filter_zeros || filter_thresh > 0) {
    # Calculate a metric to identify pathways with meaningful signal
    rowSums1 <- rowSums(mat1)
    rowSums2 <- rowSums(mat2)
    colSums1 <- colSums(mat1)
    colSums2 <- colSums(mat2)

    # Find pathways that have signal above threshold in at least one condition
    pathways_to_keep <- which(colSums1 > filter_thresh | colSums2 > filter_thresh)

    if (length(pathways_to_keep) == 0) {
      message("No pathways pass the filtering threshold. Keeping all pathways.")
    } else if (length(pathways_to_keep) < ncol(mat1)) {
      message(paste("Filtering out", ncol(mat1) - length(pathways_to_keep), "pathways with low signal."))
      # Filter matrices
      mat1 <- mat1[, pathways_to_keep, drop = FALSE]
      mat2 <- mat2[, pathways_to_keep, drop = FALSE]
      signaling <- signaling[pathways_to_keep]
    }
  }

  # Calculate the difference matrix (mat2 - mat1)
  diff_matrix_raw <- mat2 - mat1

  # Create normalized matrices for visualization (normalize by column maximum)
  mat1_norm <- sweep(mat1, 2L, apply(mat1, 2, max), "/", check.margin = FALSE)
  mat2_norm <- sweep(mat2, 2L, apply(mat2, 2, max), "/", check.margin = FALSE)

  # Replace NAs with 0s
  mat1_norm[is.na(mat1_norm)] <- 0
  mat2_norm[is.na(mat2_norm)] <- 0

  # Calculate difference matrix using normalized matrices
  if (use_log2fc) {
    # Add small pseudo-count to avoid log(0)
    pseudo_count <- 0.01
    # Calculate log2 fold change
    diff_matrix <- log2((mat2_norm + pseudo_count) / (mat1_norm + pseudo_count))
    # Cap extreme values
    cap_value <- 3
    diff_matrix[diff_matrix > cap_value] <- cap_value
    diff_matrix[diff_matrix < -cap_value] <- -cap_value
    diff_matrix[is.infinite(diff_matrix)] <- NA  # Handle any remaining Inf values
  } else {
    # Use simple difference
    diff_matrix <- mat2_norm - mat1_norm
  }

  # Set colors for the heatmap
  if (is.character(color.heatmap) && length(color.heatmap) == 1) {
    if (color.heatmap == "custom") {
      # Custom color palette: blue-white-orange
      color.heatmap.use <- colorRampPalette(c("deepskyblue", "white", "darkorange"))(100)
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
    color.use <- scPalette(length(rownames(diff_matrix)))
    names(color.use) <- rownames(diff_matrix)
  }

  # Create title with measure description
  measure_description <- switch(
    measure,
    "sender" = "Outgoing",
    "receiver" = "Incoming",
    "both" = "Overall"
  )

  if (is.null(title)) {
    title <- paste0("Differential ", measure_description, " signaling patterns")
    if (use_log2fc) {
      title <- paste0(title, " (log2FC)")
    }
  } else {
    title <- paste0("Differential ", measure_description, " signaling patterns - ", title)
    if (use_log2fc) {
      title <- paste0(title, " (log2FC)")
    }
  }

  # Create annotations
  row_anno <- ComplexHeatmap::rowAnnotation(
    Strength = ComplexHeatmap::anno_barplot(rowSums(abs(diff_matrix)), border = FALSE),
    show_annotation_name = FALSE
  )

  col_anno <- ComplexHeatmap::HeatmapAnnotation(
    Strength = ComplexHeatmap::anno_barplot(colSums(abs(diff_matrix)), border = FALSE,
                                            gp = grid::gpar(fill = color.use[1:length(colSums(abs(diff_matrix)))])),
    show_annotation_name = FALSE
  )

  # Set legend breaks (follow reference style)
  min_val <- min(diff_matrix, na.rm = TRUE)
  max_val <- max(diff_matrix, na.rm = TRUE)

  # Check if all values are the same
  if (min_val == max_val) {
    # Force different values for color scale
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
  }

  # Round to one decimal place
  legend.break <- round(legend.break, 1)

  # Ensure values are distinct
  if (length(unique(legend.break)) < 2) {
    legend.break <- c(min_val - 1, min_val, min_val + 1)
  }

  # Create the heatmap (following reference style)
  legend_title <- if(use_log2fc) "Log2FC" else "Difference"

  # Use circlize::colorRamp2 for more control over the color mapping
  col_fun <- circlize::colorRamp2(legend.break, color.heatmap.use[c(1, 50, 100)])

  # Create the heatmap
  if (transpose) {
    # For transposed version, we need to create new annotations specifically for the transposed matrix
    t_diff_matrix <- t(diff_matrix)

    # Create appropriate annotations for the transposed matrix
    t_row_anno <- ComplexHeatmap::rowAnnotation(
      Strength = ComplexHeatmap::anno_barplot(rowSums(abs(t_diff_matrix)), border = FALSE),
      show_annotation_name = FALSE
    )

    t_col_anno <- ComplexHeatmap::HeatmapAnnotation(
      Strength = ComplexHeatmap::anno_barplot(colSums(abs(t_diff_matrix)), border = FALSE,
                                              gp = grid::gpar(fill = color.use[1:length(colSums(abs(t_diff_matrix)))])),
      show_annotation_name = FALSE
    )

    # Transposed version (pathways as rows, cell types as columns)
    ht1 <- ComplexHeatmap::Heatmap(
      t_diff_matrix,
      col = col_fun,
      na_col = "white",
      name = legend_title,
      right_annotation = t_row_anno,
      top_annotation = t_col_anno,
      cluster_rows = cluster_cols,  # Swap clustering settings
      cluster_columns = cluster_rows,
      show_row_dend = show_column_dend && cluster_cols,  # Show dendrogram only if clustering and show_dend is TRUE
      show_column_dend = show_row_dend && cluster_rows,
      show_row_names = show_colnames,  # Swap name visibility settings
      show_column_names = show_rownames,
      row_names_side = "left",
      row_names_rot = 0,
      row_names_gp = grid::gpar(fontsize = font.size),
      column_names_gp = grid::gpar(fontsize = font.size),
      width = grid::unit(height, "cm"),  # Swap width and height
      height = grid::unit(width, "cm"),
      row_title = title,  # Title now for rows
      row_title_gp = grid::gpar(fontsize = font.size.title),
      column_names_rot = 90,
      row_names_max_width = grid::unit(10, "cm"),
      # Add cell function to display values
      cell_fun = if(show_values) {
        function(j, i, x, y, width, height, fill) {
          # Format the value based on specified digits
          value <- t_diff_matrix[i, j]
          if (!is.na(value) && abs(value) > 1e-10) {  # Hide values that are effectively zero
            display_val <- format(round(value, value_digits), nsmall = value_digits)
            grid::grid.text(display_val, x, y, gp = grid::gpar(fontsize = font.size * 0.8))
          }
        }
      } else NULL,
      # Control border display
      border = if(show_heatmap_border) border_color else NA,
      border_gp = grid::gpar(lwd = border_width),
      rect_gp = grid::gpar(col = if(show_heatmap_border) border_color else NA, lwd = border_width),
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
    ht1 <- ComplexHeatmap::Heatmap(
      diff_matrix,
      col = col_fun,
      na_col = "white",
      name = legend_title,
      right_annotation = row_anno,
      top_annotation = col_anno,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_cols,
      show_row_dend = show_row_dend && cluster_rows,  # Show dendrogram only if clustering and show_dend is TRUE
      show_column_dend = show_column_dend && cluster_cols,
      show_row_names = show_rownames,
      show_column_names = show_colnames,
      row_names_side = "left",
      row_names_rot = 0,
      row_names_gp = grid::gpar(fontsize = font.size),
      column_names_gp = grid::gpar(fontsize = font.size),
      width = grid::unit(height, "cm"),
      height = grid::unit(width, "cm"),
      column_title = title,
      column_title_gp = grid::gpar(fontsize = font.size.title),
      column_names_rot = 90,
      row_names_max_width = grid::unit(10, "cm"),
      # Add cell function to display values
      cell_fun = if(show_values) {
        function(j, i, x, y, width, height, fill) {
          # Format the value based on specified digits
          value <- diff_matrix[i, j]
          if (!is.na(value) && abs(value) > 1e-10) {  # Hide values that are effectively zero
            display_val <- format(round(value, value_digits), nsmall = value_digits)
            grid::grid.text(display_val, x, y, gp = grid::gpar(fontsize = font.size * 0.8))
          }
        }
      } else NULL,
      # Control border display
      border = if(show_heatmap_border) border_color else NA,
      border_gp = grid::gpar(lwd = border_width),
      rect_gp = grid::gpar(col = if(show_heatmap_border) border_color else NA, lwd = border_width),
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

  # Save the plot if requested
  if (save_plot) {
    if (!is.character(save_name)) {
      save_name <- "heatmap_diff.pdf"
    }

    tryCatch({
      pdf(save_name, width = width/2.54, height = height/2.54)
      print(ht1)
      invisible(dev.off())
      message(paste("Heatmap saved to", save_name))
    }, error = function(e) {
      warning(paste("Failed to save heatmap:", e$message))
    })
  }

  # Return results based on user preference
  if (return_data) {
    # Return a list containing the heatmap and the matrices
    result <- list(
      heatmap = ht1,
      mat1 = mat1,
      mat2 = mat2,
      mat1_norm = mat1_norm,
      mat2_norm = mat2_norm,
      diff_matrix = diff_matrix,
      diff_matrix_raw = diff_matrix_raw,
      transposed = transpose,
      pathways = colnames(mat1),
      cell_types = rownames(mat1),
      measure = measure,
      comparison = comparison,
      condition_names = condition_names[comparison]
    )

    return(result)
  } else {
    # Return just the heatmap
    return(ht1)
  }
}
