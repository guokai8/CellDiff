#' Rank and Visualize Significant Signaling Pathways Across Multiple Groups
#'
#' This function ranks signaling pathways based on their differential activity
#' across multiple groups and visualizes the significant ones.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (default NULL, uses all objects)
#' @param reference Integer or character, index or name of the reference object to compare against.
#'   If NULL (default), the first object in comparison will be used as reference.
#' @param measure Measure type, either "weight" or "count"
#' @param slot.name Slot name containing the network information
#' @param color.use Vector of colors for the groups
#' @param color.pathway Vector of colors for the pathway (slope plot)
#' @param pThresh P-value threshold for significant differences
#' @param tol Tolerance for fold change significance
#' @param top.n Number of top pathways to show
#' @param show.pval Whether to show p-values on the plot
#' @param show.heatmap Whether to generate a heatmap
#' @param title Title for the plot
#' @param sources.use Optional vector of source cell types to filter
#' @param targets.use Optional vector of target cell types to filter
#' @param comparison_method Method for comparing groups: "all_vs_ref" (each group vs reference)
#'   or "sequential" (each group vs next) or "custom_pairs" (specified pairs)
#' @param custom_comparisons List of custom comparison pairs
#' @param show.dendrogram Whether to show dendrogram in the heatmap
#' @param show.barplot Whether to include barplots showing pathway strength
#' @param cluster_groups Whether to cluster groups in the heatmap
#' @param use_log2fc Whether to use log2 fold change instead of absolute differences
#' @param min_pathways Minimum number of pathways to include in visualization
#' @param filter_min_change Minimum fold change to consider a pathway significant
#' @param order_by Method to order pathways: "foldchange", "pvalue", or "pathway_name"
#' @param show_all_pairwise Whether to show all pairwise comparisons in a grid
#' @param save_plot Whether to save the plot to file
#' @param save_name Filename for saving the plot
#' @param save_width Width of the saved plot in inches
#' @param save_height Height of the saved plot in inches
#' @param return_top_paths Whether to return only top_n pathways or all significant pathways
#' @param show_comparison_barplot Show faceted barplot for all comparisons (default: TRUE)
#' @param show_comparison_heatmap Show heatmap for all comparisons (default: FALSE)
#' @param show_comparison_slope Show slope plot for all comparisons (default: FALSE)
#' @param comparison_plot_type Default plot type for comparison: "barplot", "heatmap", or "slope"
#' @param heatmap_colors Vector of 3 colors for heatmap gradient (low, mid, high)
#' @param text_size Text size for plot elements
#' @param label_angle Angle for x-axis labels
#' @param border_color Border color for tiles in heatmap
#' @param cell.type.strategy Character, strategy for aligning cell types across objects.
#'   "shared" (default) uses only cell types common to all objects (intersection).
#'   "union" uses all cell types from all objects (fills missing with 0).
#' @param show.all Logical, whether to show all pathways (TRUE) or only significant pathways
#'   (FALSE, default) in the comparison_barplot and comparison_heatmap. When TRUE,
#'   non-significant pathways are shown with faded colors and significance stars
#'   (*p<0.05, **p<0.01, ***p<0.001) indicate which pathways are significant in each comparison.
#'
#' @return A list containing:
#'   \item{plots}{List of individual barplots for each comparison}
#'   \item{heatmaps}{List of individual heatmaps for each comparison (if requested)}
#'   \item{combined_plot}{Combined plot of all pairwise comparisons (if requested)}
#'   \item{data}{Data frame with all pathway information and statistics}
#'   \item{all_significant_paths}{List of significant pathways for each comparison}
#'   \item{all_significant_paths_full}{List of all significant pathways without top.n limit}
#'   \item{params}{List of parameters used for the analysis}
#'   \item{comparison_barplot}{Faceted barplot comparing all conditions to reference (if requested)}
#'   \item{comparison_heatmap}{Heatmap comparing all conditions to reference (if requested)}
#'   \item{comparison_slope_plot}{Slope plot showing pathway trends across conditions (if requested)}
#'   \item{top_paths}{Vector of top pathways based on frequency across comparisons}
#'
#' @examples
#' # Basic usage with default parameters
#' results <- rankDiffM(
#'   object.list = cellchat.list,
#'   comparison_method = "all_vs_ref",
#'   reference = "WT",
#'   cell.type.strategy = "shared"
#' )
#'
#' # Using log2 fold change and custom visualization options
#' results <- rankDiffM(
#'   object.list = cellchat.list,
#'   comparison_method = "all_vs_ref",
#'   reference = "WT",
#'   use_log2fc = TRUE,
#'   comparison_plot_type = "heatmap",
#'   heatmap_colors = c("blue", "white", "red"),
#'   cell.type.strategy = "union"
#' )
#'
#' # Using sequential comparison with custom pathway ordering
#' results <- rankDiffM(
#'   object.list = cellchat.list,
#'   comparison_method = "sequential",
#'   order_by = "pathway_name",
#'   show_comparison_slope = TRUE,
#'   cell.type.strategy = "shared"
#' )
#'
#' # Show all pathways (not just significant) in comparison plots
#' results <- rankDiffM(
#'   object.list = cellchat.list,
#'   comparison_method = "all_vs_ref",
#'   reference = "WT",
#'   show_comparison_barplot = TRUE,
#'   show.all = TRUE
#' )
#'
#' @export
rankDiffM <- function(object.list, comparison = NULL, reference = NULL,
                      measure = "weight", slot.name = "netP", color.use = NULL,
                      color.pathway = NULL,
                      pThresh = 0.05, tol = 0.05, top.n = NULL,
                      show.pval = TRUE, show.heatmap = FALSE, title = NULL,
                      sources.use = NULL, targets.use = NULL,
                      comparison_method = c("all_vs_ref", "sequential", "custom_pairs"),
                      custom_comparisons = NULL, show.dendrogram = FALSE,
                      show.barplot = TRUE, cluster_groups = FALSE,
                      use_log2fc = FALSE, min_pathways = 3,
                      filter_min_change = 0.25, order_by = c("foldchange", "pvalue", "pathway_name"),
                      show_all_pairwise = FALSE, save_plot = FALSE,
                      save_name = "rankDiff_plot.pdf", save_width = 10, save_height = 8,
                      return_top_paths = TRUE,
                      # New parameters for multi-condition visualization
                      show_comparison_barplot = TRUE,    # Show faceted barplot for all comparisons
                      show_comparison_heatmap = FALSE,   # Show heatmap for all comparisons
                      show_comparison_slope = FALSE,     # Show slope plot for all comparisons
                      comparison_plot_type = c("barplot", "heatmap", "slope"), # Default plot type for comparison
                      heatmap_colors = c("deepskyblue", "white", "darkorange"), # Colors for heatmap
                      text_size = 10,                   # Text size for plot elements
                      label_angle = 45,                 # Angle for x-axis labels
                      border_color = "grey60",          # Border color for tiles in heatmap
                      cell.type.strategy = c("shared", "union"),  # Strategy for aligning cell types
                      show.all = FALSE) {  # Show all pathways in comparison plots

  # Load required packages
  require(plyr)
  require(reshape2)

  # Match arguments
  comparison_method <- match.arg(comparison_method)
  order_by <- match.arg(order_by)
  comparison_plot_type <- match.arg(comparison_plot_type)
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

  # Basic validation of CellChat objects
  if (!is.list(object.list) || length(object.list) < max(comparison)) {
    stop("object.list must be a list containing at least ", max(comparison), " CellChat objects")
  }

  # Get condition names
  condition_names <- names(object.list)
  if (is.null(condition_names) || any(condition_names == "")) {
    condition_names <- paste0("Condition_", 1:length(object.list))
    names(object.list) <- condition_names
  }

  # Set group colors if not provided
  if (is.null(color.use)) {
    if (exists("global_colors") && length(global_colors) >= length(comparison)) {
      color.use <- global_colors[1:length(comparison)]
    } else {
      color.use <- rainbow(length(comparison))
    }
    names(color.use) <- condition_names[comparison]
  } else if (length(color.use) < length(comparison)) {
    # Recycle colors if needed
    color.use <- rep(color.use, length.out = length(comparison))
    names(color.use) <- condition_names[comparison]
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

      # Convert object indices to positions within comparison vector
      pos1 <- which(comparison == pair[1])
      pos2 <- which(comparison == pair[2])
      comparison_pairs[[i]] <- c(pos1, pos2)
    }
  }

  if (length(comparison_pairs) == 0) {
    stop("No valid comparison pairs were generated")
  }

  # Align cell types across all objects in comparison
  alignment <- alignCellTypes(object.list, indices = comparison, strategy = cell.type.strategy)
  aligned_cells <- alignment$cell_types

  # Initialize lists to store results
  all_results <- list()
  all_plots <- list()
  all_heatmaps <- list()
  all_significant_paths <- list()
  all_significant_paths_full <- list()
  all_df_combined <- list()  # Store original contribution data for comparison plots
  # Process each comparison pair
  for (p in 1:length(comparison_pairs)) {
    pair <- comparison_pairs[[p]]
    idx1 <- comparison[pair[1]]
    idx2 <- comparison[pair[2]]

    # Skip if indices are invalid
    if (idx1 > length(object.list) || idx2 > length(object.list)) {
      warning(paste("Invalid comparison pair:", idx1, "vs", idx2, "- skipping"))
      next
    }

    # Process the data for comparison
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx.list <- list()
    pSum.original.all <- c()
    object.names.pair <- condition_names[c(idx1, idx2)]

    for (i in 1:2) {
      idx <- if (i == 1) idx1 else idx2
      object.data <- object.list[[idx]]@netP

      prob <- object.data$prob
      # Apply source/target filtering if provided
      # Apply source/target filtering if provided
      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (all(sources.use %in% rownames(prob))) {
            sources.idx <- match(sources.use, rownames(prob))
          } else {
            stop("The input `sources.use` should be cell group names or a numerical vector!")
          }
        } else {
          sources.idx <- sources.use
        }
        exclude.idx <- setdiff(1:nrow(prob), sources.idx)
        prob[exclude.idx, , ] <- 0
      }

      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (all(targets.use %in% rownames(prob))) {
            targets.idx <- match(targets.use, rownames(prob))
          } else {
            stop("The input `targets.use` should be cell group names or a numerical vector!")
          }
        } else {
          targets.idx <- targets.use
        }
        exclude.idx <- setdiff(1:nrow(prob), targets.idx)
        prob[, exclude.idx, ] <- 0
      }

      if (measure == "count") {
        prob <- 1 * (prob > 0)
      }

      prob.list[[i]] <- prob

      if (sum(prob) == 0) {
        warning(paste("No inferred communications for object", i, "in pair", p))
        next
      }

      pSum.original[[i]] <- apply(prob, 3, sum)

      if (measure == "weight") {
        pSum[[i]] <- -1/log(pSum.original[[i]] + 1e-10)
        pSum[[i]][is.na(pSum[[i]])] <- 0
        idx.list[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)

        if (length(idx.list[[i]]) > 0) {
          pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx.list[[i]]])
        }
      } else if (measure == "count") {
        pSum[[i]] <- pSum.original[[i]]
      }

      pair.name[[i]] <- names(pSum.original[[i]])
    }

    # Replace values with max pSum for infinity cases
    if (measure == "weight" && length(unlist(idx.list)) > 0) {
      values.assign <- seq(max(unlist(pSum), na.rm = TRUE) * 1.1,
                           max(unlist(pSum), na.rm = TRUE) * 1.5,
                           length.out = length(unlist(idx.list)))
      position <- sort(pSum.original.all, index.return = TRUE)$ix

      for (i in 1:2) {
        if (length(idx.list[[i]]) > 0) {
          if (i == 1) {
            pSum[[i]][idx.list[[i]]] <- values.assign[match(1:length(idx.list[[i]]), position)]
          } else {
            pSum[[i]][idx.list[[i]]] <- values.assign[match(length(unlist(idx.list[1:i-1])) + 1:length(idx.list[[i]]), position)]
          }
        }
      }
    }

    # Extract all pathway names
    pair.name.all <- as.character(unique(unlist(pair.name)))

    if (length(pair.name.all) == 0) {
      warning(paste("No pathway names found for pair", p))
      next
    }

    # Create data frames with pathway contributions
    df.list <- list()
    for (i in 1:2) {
      df.list[[i]] <- data.frame(
        name = pair.name.all,
        contribution = 0,
        contribution.scaled = 0,
        group = object.names.pair[i],
        stringsAsFactors = FALSE
      )

      matching_names <- intersect(pair.name[[i]], pair.name.all)

      if (length(matching_names) > 0) {
        df.list[[i]][match(matching_names, df.list[[i]]$name), "contribution.scaled"] <- pSum[[i]][matching_names]
        df.list[[i]][match(matching_names, df.list[[i]]$name), "contribution"] <- pSum.original[[i]][matching_names]
      }
    }

    # Combine data frames
    df_combined <- do.call(rbind, df.list)

    # Filter out pathways with zero contribution
    pathway_sums <- tapply(df_combined$contribution, df_combined$name, sum)
    valid_pathways <- names(pathway_sums[pathway_sums > 0])

    if (length(valid_pathways) == 0) {
      warning(paste("No pathways with non-zero contribution for pair", p))
      next
    }

    df_combined <- df_combined[df_combined$name %in% valid_pathways, ]

    # Set factor levels for proper ordering
    df_combined$group <- factor(df_combined$group, levels = object.names.pair)

    # Calculate aggregate values by group
    df_aggregated <- aggregate(
      cbind(contribution, contribution.scaled) ~ group + name,
      data = df_combined,
      FUN = mean
    )

    # Calculate relative changes
    wide_df <- reshape2::dcast(
      df_aggregated[, c("group", "name", "contribution")],
      name ~ group,
      value.var = "contribution"
    )

    # Calculate fold change or log2 fold change
    if (use_log2fc) {
      # Avoid division by zero with pseudocount
      pseudo_count <- 1e-10
      wide_df$fold_change <- log2((wide_df[[object.names.pair[2]]] + pseudo_count) /
                                    (wide_df[[object.names.pair[1]]] + pseudo_count))
    } else {
      # Simple fold change
      wide_df$fold_change <- wide_df[[object.names.pair[2]]] /
        (wide_df[[object.names.pair[1]]] + 1e-10)
    }
    # Calculate p-values with Wilcoxon test
    p_values <- numeric(length(valid_pathways))
    names(p_values) <- valid_pathways

    for (i in 1:length(valid_pathways)) {
      pathway <- valid_pathways[i]

      # Extract values from each condition separately to handle different dimensions
      prob_values_list <- list()

      for (j in 1:2) {
        if (pathway %in% pair.name[[j]]) {
          # Extract all values for this pathway and flatten to vector
          prob_values_list[[j]] <- as.vector(prob.list[[j]][, , pathway])
        } else {
          prob_values_list[[j]] <- numeric(0)
        }
      }

      # Find the maximum length to create the matrix
      max_len <- max(length(prob_values_list[[1]]), length(prob_values_list[[2]]))

      # Create matrix with appropriate size
      prob_values <- matrix(NA, nrow = max_len, ncol = 2)

      # Fill in the values
      if (length(prob_values_list[[1]]) > 0) {
        prob_values[1:length(prob_values_list[[1]]), 1] <- prob_values_list[[1]]
      }
      if (length(prob_values_list[[2]]) > 0) {
        prob_values[1:length(prob_values_list[[2]]), 2] <- prob_values_list[[2]]
      }

      # Remove rows with all zeros or all NAs
      prob_values <- prob_values[rowSums(prob_values, na.rm = TRUE) != 0, , drop = FALSE]

      # Check if paired test is applicable
      paired_test <- TRUE
      if (nrow(prob.list[[1]]) != nrow(prob.list[[2]]) ||
          ncol(prob.list[[1]]) != ncol(prob.list[[2]])) {
        paired_test <- FALSE
      }

      if (nrow(prob_values) > 3 && sum(is.na(prob_values)) == 0) {
        p_val <- suppressWarnings(
          wilcox.test(prob_values[, 1], prob_values[, 2], paired = paired_test)$p.value
        )
        p_values[pathway] <- p_val
      } else {
        p_values[pathway] <- 0
      }
    }

    # Replace NA p-values with 0
    p_values[is.na(p_values)] <- 0

    # Add p-values to the data frame
    wide_df$pvalues <- p_values[wide_df$name]

    # Apply multiple testing correction
    wide_df$padj <- p.adjust(wide_df$pvalues, method = "BH")

    # Identify significant pathways
    significant <- if (use_log2fc) {
      (abs(wide_df$fold_change) > filter_min_change) & (wide_df$padj < pThresh)
    } else {
      (wide_df$fold_change < (1 - tol) | wide_df$fold_change > (1 + tol)) &
        (wide_df$padj < pThresh)
    }

    # Add a flag to indicate which pathways are unique to each group
    wide_df$group1_only <- wide_df$name %in% pair.name[[1]] & !(wide_df$name %in% pair.name[[2]])
    wide_df$group2_only <- wide_df$name %in% pair.name[[2]] & !(wide_df$name %in% pair.name[[1]])

    # Include both significant paths and unique paths
    significant_paths <- unique(c(
      wide_df$name[significant],
      wide_df$name[wide_df$group1_only | wide_df$group2_only]
    ))

    if (length(significant_paths) == 0) {
      warning(paste("No significant pathways found for pair", p, "with the current thresholds"))
      significant_paths <- character(0)
    }

    # Order pathways by fold change, p-value, or name
    if (order_by == "foldchange") {
      if (use_log2fc) {
        wide_df <- wide_df[order(abs(wide_df$fold_change), decreasing = TRUE), ]
      } else {
        # For regular fold change, order by distance from 1 (no change)
        wide_df <- wide_df[order(abs(wide_df$fold_change - 1), decreasing = TRUE), ]
      }
    } else if (order_by == "pvalue") {
      wide_df <- wide_df[order(wide_df$padj), ]
    } else if (order_by == "pathway_name") {
      wide_df <- wide_df[order(wide_df$name), ]
    }

    # Store results
    all_results[[p]] <- wide_df
    all_significant_paths_full[[p]] <- significant_paths

    # Store df_combined for comparison plots (needed to access all pathways)
    all_df_combined[[p]] <- df_combined

    # Don't limit significant_paths for plotting
    plot_significant_paths <- significant_paths

    # For the actual results output, still limit if requested
    if (!is.null(top.n) && top.n > 0 && length(significant_paths) > top.n) {
      # Limit to top pathways based on significance and fold change
      path_order <- order(
        !significant[match(significant_paths, wide_df$name)],
        !(wide_df$group1_only | wide_df$group2_only)[match(significant_paths, wide_df$name)],
        -abs(wide_df$fold_change[match(significant_paths, wide_df$name)])
      )

      significant_paths <- significant_paths[path_order][1:top.n]
    } else if (length(significant_paths) < min_pathways && nrow(wide_df) >= min_pathways) {
      # If too few significant paths, use top min_pathways by fold change
      significant_paths <- wide_df$name[1:min(min_pathways, nrow(wide_df))]
      warning(paste("Using top", length(significant_paths),
                    "pathways by fold change for pair", p,
                    "as too few significant pathways were found"))
    }

    # Store for final output
    all_significant_paths[[p]] <- significant_paths

    # Filter data for plotting - use all significant pathways for visualization
    plot_data <- df_combined[df_combined$name %in% plot_significant_paths, ]

    if (nrow(plot_data) == 0) {
      warning(paste("No data available for plotting for pair", p))
      all_plots[[p]] <- NULL
      all_heatmaps[[p]] <- NULL
      next
    }

    # Calculate fold changes for ordering plot
    fc_values <- numeric(length(plot_significant_paths))
    names(fc_values) <- plot_significant_paths

    # Calculate fold change for each pathway based on the mean contribution.scaled values
    for (pathway in plot_significant_paths) {
      # Calculate means for each group
      group1_values <- df_combined$contribution.scaled[df_combined$name == pathway &
                                                         df_combined$group == object.names.pair[1]]
      group2_values <- df_combined$contribution.scaled[df_combined$name == pathway &
                                                         df_combined$group == object.names.pair[2]]

      group1_mean <- mean(group1_values, na.rm = TRUE)
      group2_mean <- mean(group2_values, na.rm = TRUE)

      # Calculate fold change
      if (length(group1_values) > 0 && length(group2_values) > 0 &&
          !is.na(group1_mean) && !is.na(group2_mean)) {
        # Add small epsilon to avoid division by zero
        fc_values[pathway] <- (group2_mean + 1e-10) / (group1_mean + 1e-10)
      } else {
        # If a pathway is missing in one group, set a default value
        fc_values[pathway] <- 1  # No change
      }
    }

    # Use log2 of fold changes for better comparison
    log2_fc <- log2(fc_values)

    # Order by absolute log2 fold change (most changed first)
    ordered_names <- names(log2_fc)[order(abs(log2_fc), decreasing = TRUE)]

    # Update plot_data with the custom ordering
    plot_data$name <- factor(as.character(plot_data$name), levels = ordered_names)

    # Create text colors based on significance and direction
    colors.text <- rep("black", length(ordered_names))
    names(colors.text) <- ordered_names

    for (i in 1:length(ordered_names)) {
      path <- ordered_names[i]
      row_idx <- which(wide_df$name == path)

      if (length(row_idx) > 0) {
        p_val <- wide_df$padj[row_idx]

        if (use_log2fc) {
          # Log2FC interpretation
          fc <- wide_df$fold_change[row_idx]

          if ((fc < -filter_min_change) && (p_val < pThresh)) {
            colors.text[i] <- color.use[object.names.pair[1]]  # First group is stronger
          } else if ((fc > filter_min_change) && (p_val < pThresh)) {
            colors.text[i] <- color.use[object.names.pair[2]]  # Second group is stronger
          }
        } else {
          # Regular fold change interpretation
          ratio <- wide_df$fold_change[row_idx]

          if ((ratio < (1 - tol)) && (p_val < pThresh)) {
            colors.text[i] <- color.use[object.names.pair[1]]  # First group is stronger
          } else if ((ratio > (1 + tol)) && (p_val < pThresh)) {
            colors.text[i] <- color.use[object.names.pair[2]]  # Second group is stronger
          }
        }
      }
    }

    # Create the barplot
    pair_title <- paste(object.names.pair[2], "vs", object.names.pair[1])
    if (!is.null(title)) {
      if (length(title) == 1) {
        pair_title <- paste0(title, " (", pair_title, ")")
      } else if (length(title) >= p) {
        pair_title <- title[p]
      }
    }

    gg <- ggplot2::ggplot(plot_data,
                          ggplot2::aes(x = name, y = contribution, fill = group)) +
      ggplot2::geom_bar(stat = "identity", position = "fill") +
      ggplot2::xlab("") +
      ggplot2::ylab("Relative information flow") +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                          color = "grey50", size = 0.5) +
      ggplot2::theme_classic() +
      ggplot2::coord_flip() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(
          colour = colors.text[levels(plot_data$name)],
          size = text_size
        ),
        legend.position = "bottom"
      ) +
      ggplot2::scale_fill_manual(values = color.use[object.names.pair]) +
      ggplot2::ggtitle(pair_title)

    # Add p-values to the plot if requested
    if (show.pval && length(plot_significant_paths) > 0) {
      # Prepare p-value annotations
      pval_data <- wide_df[wide_df$name %in% plot_significant_paths, ]

      # Convert p-values to stars
      pval_labels <- sapply(pval_data$padj, function(p) {
        if (is.na(p)) return("")
        if (p < 0.001) return("***")
        if (p < 0.01) return("**")
        if (p < 0.05) return("*")
        return("")
      })

      # Add text labels to the plot - only if we have at least one significant p-value
      if (any(pval_labels != "")) {
        # Create label data
        label_data <- data.frame(
          name = pval_data$name,
          contribution = 0.95,  # Position near the top
          label = pval_labels,
          group = object.names.pair[1]  # Doesn't matter which group, just for position
        )

        # Make sure name is a factor with same levels as plot_data
        label_data$name <- factor(label_data$name, levels = levels(plot_data$name))

        gg <- gg +
          ggplot2::geom_text(
            data = label_data,
            ggplot2::aes(label = label),
            position = ggplot2::position_dodge(width = 1),
            hjust = 1
          )
      }
    }

    # Store the plot
    all_plots[[p]] <- gg

    # Create heatmap if requested
    if (show.heatmap && length(plot_significant_paths) >= 2) {  # Need at least 2 for clustering
      # Prepare data for heatmap
      heatmap_data <- reshape2::dcast(
        plot_data,
        name ~ group,
        value.var = "contribution"
      )

      # Set rownames for the heatmap
      rownames(heatmap_data) <- heatmap_data$name
      heatmap_data$name <- NULL

      # Create annotation for groups
      annotation_col <- data.frame(
        Group = factor(colnames(heatmap_data))
      )
      rownames(annotation_col) <- colnames(heatmap_data)

      # Create heatmap - safely check if pheatmap is available
      heatmap <- NULL
      if (requireNamespace("pheatmap", quietly = TRUE)) {
        tryCatch({
          heatmap <- pheatmap::pheatmap(
            heatmap_data,
            cluster_cols = cluster_groups,
            cluster_rows = show.dendrogram,
            scale = "row",
            annotation_col = annotation_col,
            annotation_colors = list(Group = color.use[object.names.pair]),
            border_color = "white",
            color = colorRampPalette(heatmap_colors)(128),
            main = paste0("Heatmap: ", pair_title),
            fontsize_row = text_size,
            fontsize_col = text_size,
            silent = TRUE
          )
        }, error = function(e) {
          warning(paste("Error creating heatmap for pair", p, ":", e$message))
          heatmap <- NULL
        })
      }

      all_heatmaps[[p]] <- heatmap
    }
  }

  # Add names to lists for better readability
  comparison_names <- sapply(1:length(comparison_pairs), function(p) {
    pair <- comparison_pairs[[p]]
    idx1 <- comparison[pair[1]]
    idx2 <- comparison[pair[2]]
    paste(condition_names[idx2], "vs", condition_names[idx1])
  })
  names(all_results) <- comparison_names
  names(all_significant_paths_full) <- comparison_names
  names(all_significant_paths) <- comparison_names
  if (length(all_plots) > 0) names(all_plots) <- comparison_names
  if (length(all_heatmaps) > 0) names(all_heatmaps) <- comparison_names
  names(all_df_combined) <- comparison_names

  create_combined_comparison_plot <- function(all_significant_paths_full, all_results,
                                              all_df_combined,
                                              comparison_pairs, condition_names, comparison,
                                              reference, color.use,
                                              text_size = 10,
                                              label_angle = 45,
                                              heatmap_colors = c("deepskyblue", "white", "darkorange"),
                                              border_color = "grey60",
                                              pThresh = 0.05,
                                              comparison_method = "all_vs_ref",
                                              color.pathway = NULL,
                                              show.all = FALSE) {
    # Debug output
    if (getOption("CellDiff.debug", FALSE)) {
      cat(sprintf("\n[COMBINED PLOT] show.all = %s\n", show.all))
    }

    # Get pathways to include in the plot
    if (show.all) {
      # Get all pathways from all_results (includes non-significant pathways)
      all_paths <- unique(unlist(lapply(all_results, function(x) {
        if (!is.null(x)) x$name else NULL
      })))
      if (getOption("CellDiff.debug", FALSE)) {
        cat(sprintf("[COMBINED PLOT] show.all=TRUE: Using all pathways from all_results (%d pathways)\n", length(all_paths)))
      }
    } else {
      # Get only significant pathways
      all_paths <- unique(unlist(all_significant_paths_full))
      if (getOption("CellDiff.debug", FALSE)) {
        cat(sprintf("[COMBINED PLOT] show.all=FALSE: Using only significant pathways (%d pathways)\n", length(all_paths)))
      }
    }

    if (length(all_paths) == 0) {
      return(NULL)
    }

    # Create a data frame for the combined plot
    ref_idx <- which(comparison == reference)
    ref_name <- condition_names[comparison[ref_idx]]

    # Collect contribution data for all pathways in all comparisons
    plot_data <- data.frame()

    for (p in 1:length(comparison_pairs)) {
      if (is.null(all_df_combined[[p]])) next

      pair <- comparison_pairs[[p]]

      # For all_vs_ref: only include comparisons against the reference
      # For custom_pairs: include all comparison pairs
      if (comparison_method == "all_vs_ref" && pair[1] != ref_idx) next

      cond_idx1 <- pair[1]
      cond_idx2 <- pair[2]
      cond_name1 <- condition_names[comparison[cond_idx1]]
      cond_name2 <- condition_names[comparison[cond_idx2]]

      # Get the contribution data for this comparison pair
      df_comb <- all_df_combined[[p]]

      # Filter to include only the pathways we want to show
      # When show.all = FALSE, only include pathways significant in THIS comparison
      if (show.all) {
        # Show all pathways
        df_comb_filtered <- df_comb[df_comb$name %in% all_paths, ]
        paths_to_show <- all_paths
      } else {
        # Only show pathways significant in THIS specific comparison
        paths_to_show <- all_significant_paths_full[[p]]
        df_comb_filtered <- df_comb[df_comb$name %in% paths_to_show, ]
      }

      if (nrow(df_comb_filtered) == 0) next

      # Add comparison label and significance information
      for (path in paths_to_show) {
        path_data <- df_comb_filtered[df_comb_filtered$name == path, ]

        if (nrow(path_data) == 0) next

        # Get p-value for this pathway from the CURRENT comparison
        p_val <- NA
        is_in_results <- FALSE
        if (!is.null(all_results[[p]]) && path %in% all_results[[p]]$name) {
          p_val <- all_results[[p]]$padj[all_results[[p]]$name == path]
          is_in_results <- TRUE
        }

        # Add significance flag and rename columns for consistency
        # A pathway is significant in this comparison ONLY if:
        # 1. It exists in all_results for this comparison, AND
        # 2. Its adjusted p-value is below the threshold
        path_data$significant <- is_in_results && !is.na(p_val) && p_val < pThresh
        path_data$comparison <- paste(cond_name2, "vs", cond_name1)
        path_data$condition <- path_data$group
        path_data$pathway <- path_data$name  # Rename for consistency

        # Calculate proportions
        total_contrib <- sum(path_data$contribution)
        if (total_contrib > 0) {
          path_data$proportion <- path_data$contribution / total_contrib
        } else {
          path_data$proportion <- 0.5  # Equal if both are zero
        }

        # Ensure pathway column is character, not factor
        path_data$pathway <- as.character(path_data$pathway)
        path_data$comparison <- as.character(path_data$comparison)

        # Use rbind.fill to handle potentially missing columns
        if (nrow(plot_data) == 0) {
          plot_data <- path_data
        } else {
          plot_data <- plyr::rbind.fill(plot_data, path_data)
        }
      }
    }

    if (nrow(plot_data) == 0) {
      return(NULL)
    }

    # Ensure all required columns exist and have consistent types
    if (!"pathway" %in% colnames(plot_data)) {
      stop("Internal error: pathway column missing from plot_data")
    }
    if (!"proportion" %in% colnames(plot_data)) {
      stop("Internal error: proportion column missing from plot_data")
    }

    # Remove any rows with NA pathway or proportion
    plot_data <- plot_data[!is.na(plot_data$pathway) & !is.na(plot_data$proportion), ]

    if (nrow(plot_data) == 0) {
      return(NULL)
    }

    # Debug: check data structure
    if (length(unique(sapply(plot_data, length))) > 1) {
      # Columns have different lengths - this is the problem!
      col_lengths <- sapply(plot_data, length)
      bad_cols <- names(col_lengths[col_lengths != nrow(plot_data)])
      warning(paste("Columns with inconsistent length:", paste(bad_cols, collapse=", ")))
      # Fix by removing problematic columns or padding
      for (col in bad_cols) {
        if (length(plot_data[[col]]) < nrow(plot_data)) {
          plot_data[[col]] <- rep(plot_data[[col]][1], nrow(plot_data))
        }
      }
    }

    # Convert to data.frame if needed (plyr::rbind.fill returns a data.frame)
    plot_data <- as.data.frame(plot_data, stringsAsFactors = FALSE)

    # Order pathways by the magnitude of change
    # Calculate average difference from 0.5 (equal proportions)
    pathway_diff <- aggregate(
      abs(proportion - 0.5) ~ pathway,
      data = plot_data,
      FUN = mean
    )

    # Sort pathways by this difference
    ordered_paths <- pathway_diff$pathway[order(pathway_diff$`abs(proportion - 0.5)`, decreasing = TRUE)]

    # Set factor levels for proper ordering
    plot_data$pathway <- factor(plot_data$pathway, levels = ordered_paths)
    plot_data$condition <- factor(plot_data$condition)

    # Get all unique groups for proper factor ordering
    all_groups <- unique(plot_data$group)
    plot_data$group <- factor(plot_data$group, levels = all_groups)

    # Create a faceted barplot showing all comparisons
    # Determine the title based on comparison method
    barplot_title <- if (comparison_method == "all_vs_ref") {
      paste("Pathway Activity: All Conditions vs", ref_name)
    } else {
      "Pathway Activity: Custom Comparisons"
    }

    # Add alpha channel to show non-significant pathways with transparency
    if (show.all) {
      # Calculate significance per pathway-comparison combination
      plot_data$is_sig_in_comparison <- FALSE
      for (i in 1:nrow(plot_data)) {
        # Check if this pathway is significant in this specific comparison
        same_pathway_comp <- plot_data$pathway == plot_data$pathway[i] &
                            plot_data$comparison == plot_data$comparison[i]
        plot_data$is_sig_in_comparison[i] <- any(plot_data$significant[same_pathway_comp], na.rm = TRUE)
      }
      plot_data$bar_alpha <- ifelse(plot_data$is_sig_in_comparison, 1.0, 0.3)
    } else {
      plot_data$bar_alpha <- 1.0
    }

    barplot <- ggplot2::ggplot(plot_data,
                               ggplot2::aes(x = pathway, y = proportion, fill = group, alpha = bar_alpha)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_alpha_identity() +
      ggplot2::facet_grid(. ~ comparison) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = color.use) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = text_size),
        strip.text = ggplot2::element_text(size = text_size),
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = label_angle, hjust = 1)
      ) +
      ggplot2::labs(
        title = barplot_title,
        subtitle = if (show.all) "Faded bars: non-significant in that comparison. Stars: *p<0.05, **p<0.01, ***p<0.001" else NULL,
        x = "",
        y = "Relative Contribution"
      )

    # Add significance stars if showing all pathways
    if (show.all) {
      # Create significance labels for each pathway-comparison combination
      sig_labels_data <- data.frame()

      for (path in levels(plot_data$pathway)) {
        for (comp in unique(plot_data$comparison)) {
          # Get significance for this pathway in this comparison
          sig_rows <- plot_data[plot_data$pathway == path &
                                plot_data$comparison == comp, ]

          if (nrow(sig_rows) > 0) {
            is_sig <- any(sig_rows$significant)

            # Get the adjusted p-value from all_results
            # Need to find the correct comparison pair that matches this comparison string
            # Extract condition names from comparison string
            comp_parts <- strsplit(as.character(comp), " vs ")[[1]]
            cond_name2 <- comp_parts[1]
            cond_name1 <- comp_parts[2]

            # Find the comparison pair index that matches these conditions
            p_val <- NA
            matched_p_idx <- NA
            for (p_idx in 1:length(comparison_pairs)) {
              pair <- comparison_pairs[[p_idx]]
              pair_cond1 <- condition_names[comparison[pair[1]]]
              pair_cond2 <- condition_names[comparison[pair[2]]]

              # Check if this pair matches the comparison we're looking for
              if (pair_cond1 == cond_name1 && pair_cond2 == cond_name2) {
                # Found the right comparison pair
                matched_p_idx <- p_idx

                # IMPORTANT: Check if pathway is in all_significant_paths_full for this comparison
                # This ensures we only show stars for pathways that passed BOTH p-value AND fold-change thresholds
                is_in_significant_list <- path %in% all_significant_paths_full[[p_idx]]

                res <- all_results[[p_idx]]
                if (!is.null(res) && path %in% res$name && is_in_significant_list) {
                  # Only extract p-value if pathway is truly significant (in the significant list)
                  p_val <- res$padj[res$name == path]

                  # Debug output (only show if pathway will get a star or has low p-value but not in sig list)
                  if (getOption("CellDiff.debug", FALSE) && !is.na(p_val) && p_val < pThresh) {
                    cat(sprintf("  STAR DEBUG [%s in %s]: p_idx=%d, in_sig_list=TRUE, p_val=%.4f, will add star=%s\n",
                               path, comp, p_idx, p_val, ifelse(!is.na(p_val) && p_val < pThresh, "YES", "NO")))
                  }
                } else {
                  # Debug: pathway has low p-value but not in significant list (potential issue)
                  if (getOption("CellDiff.debug", FALSE)) {
                    in_res <- !is.null(res) && path %in% res$name
                    # Only show if pathway is in all_results (has p-value) but not in sig list
                    if (in_res) {
                      p_check <- res$padj[res$name == path]
                      if (!is.na(p_check) && p_check < pThresh) {
                        cat(sprintf("  STAR DEBUG [%s in %s]: p_idx=%d, in_all_results=TRUE, in_sig_list=FALSE, no star (failed fold-change threshold)\n",
                                   path, comp, p_idx))
                      }
                    }
                  }
                }
                break
              }
            }

            # Convert p-value to stars
            sig_label <- ""
            if (!is.na(p_val)) {
              if (p_val < 0.001) {
                sig_label <- "***"
              } else if (p_val < 0.01) {
                sig_label <- "**"
              } else if (p_val < 0.05) {
                sig_label <- "*"
              }
            }

            if (sig_label != "") {
              sig_labels_data <- rbind(sig_labels_data,
                                      data.frame(
                                        pathway = path,
                                        comparison = comp,
                                        proportion = 1.02,  # Position slightly above the bar
                                        label = sig_label,
                                        stringsAsFactors = FALSE
                                      ))
            }
          }
        }
      }

      # Add significance stars to the plot
      if (nrow(sig_labels_data) > 0) {
        sig_labels_data$pathway <- factor(sig_labels_data$pathway, levels = levels(plot_data$pathway))

        # Position stars at the end of the stacked bar (y=1.0) plus a small offset
        barplot <- barplot +
          ggplot2::geom_text(
            data = sig_labels_data,
            ggplot2::aes(x = pathway, y = 1.02, label = label),
            inherit.aes = FALSE,
            hjust = 0,
            vjust = 0.5,
            size = text_size / 3,
            fontface = "bold"
          )
      }
    }

    # Create heatmap
    # Create a heatmap showing pathway changes across conditions
    heatmap_data <- data.frame()

    if (comparison_method == "all_vs_ref") {
      # For all_vs_ref: show differences from reference
      all_conditions <- unique(plot_data$condition[plot_data$condition != ref_name])

      for (path in levels(plot_data$pathway)) {
        for (cond in all_conditions) {
          # Get reference proportion for this specific comparison
          comp_name <- paste(cond, "vs", ref_name)

          # Check if this comparison exists in the data
          if (comp_name %in% unique(plot_data$comparison)) {
            # Get reference proportion
            ref_props <- plot_data$proportion[plot_data$pathway == path &
                                                plot_data$group == ref_name &
                                                plot_data$comparison == comp_name]

            # Get condition proportion
            cond_props <- plot_data$proportion[plot_data$pathway == path &
                                                 plot_data$group == cond &
                                                 plot_data$comparison == comp_name]

            # Only proceed if we have both values
            if (length(ref_props) > 0 && length(cond_props) > 0) {
              # Calculate difference from reference
              diff_from_ref <- cond_props - ref_props

              # Get significance
              is_sig <- plot_data$significant[plot_data$pathway == path &
                                                plot_data$group == cond &
                                                plot_data$comparison == comp_name]

              # Add to heatmap data
              heatmap_data <- rbind(heatmap_data,
                                    data.frame(
                                      pathway = path,
                                      condition = cond,
                                      difference = diff_from_ref,
                                      significant = is_sig,
                                      stringsAsFactors = FALSE
                                    ))
            }
          }
        }
      }
    } else {
      # For custom_pairs: show differences for each comparison
      all_comparisons <- unique(plot_data$comparison)

      for (path in levels(plot_data$pathway)) {
        for (comp in all_comparisons) {
          # Get the two conditions for this comparison
          comp_data <- plot_data[plot_data$comparison == comp & plot_data$pathway == path, ]

          if (nrow(comp_data) == 2) {
            # Calculate difference (second condition - first condition)
            diff_value <- comp_data$proportion[2] - comp_data$proportion[1]

            # Get significance
            is_sig <- any(comp_data$significant)

            # Add to heatmap data
            heatmap_data <- rbind(heatmap_data,
                                  data.frame(
                                    pathway = path,
                                    condition = comp,
                                    difference = diff_value,
                                    significant = is_sig,
                                    stringsAsFactors = FALSE
                                  ))
          }
        }
      }
    }

    # Only create heatmap if we have data
    if (nrow(heatmap_data) > 0) {
      heatmap_title <- if (comparison_method == "all_vs_ref") {
        paste("Pathway Changes Relative to", ref_name)
      } else {
        "Pathway Changes: Custom Comparisons"
      }

      heatmap <- ggplot2::ggplot(heatmap_data,
                                 ggplot2::aes(x = condition, y = pathway, fill = difference)) +
        ggplot2::geom_tile(color = border_color) +
        ggplot2::scale_fill_gradient2(
          low = heatmap_colors[1], mid = heatmap_colors[2], high = heatmap_colors[3],
          midpoint = 0,
          name = "Difference"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = text_size),
          axis.text.x = ggplot2::element_text(angle = label_angle, hjust = 1),
          panel.grid.major = ggplot2::element_blank()
        ) +
        ggplot2::labs(
          title = heatmap_title,
          x = "",
          y = ""
        )

      # Add significance stars to heatmap if showing all pathways
      if (show.all) {
        # Get p-values for each pathway-condition combination and create labels
        heatmap_data$sig_label <- ""

        for (i in 1:nrow(heatmap_data)) {
          path <- as.character(heatmap_data$pathway[i])
          cond <- as.character(heatmap_data$condition[i])

          # Find the corresponding p-value from the correct comparison
          # The heatmap shows differences, so we need to find which comparison this belongs to
          p_val <- NA
          for (p_idx in 1:length(comparison_pairs)) {
            pair <- comparison_pairs[[p_idx]]
            pair_cond1 <- condition_names[comparison[pair[1]]]
            pair_cond2 <- condition_names[comparison[pair[2]]]

            # For all_vs_ref, condition is the non-reference condition
            # For custom_pairs, condition is the comparison label
            if (comparison_method == "all_vs_ref") {
              # Check if this is the right comparison (non-ref condition matches)
              if (cond == pair_cond2 || cond == pair_cond1) {
                # IMPORTANT: Check if pathway is in all_significant_paths_full for this comparison
                is_in_significant_list <- path %in% all_significant_paths_full[[p_idx]]
                res <- all_results[[p_idx]]
                if (!is.null(res) && path %in% res$name && is_in_significant_list) {
                  p_val <- res$padj[res$name == path]
                  break
                }
              }
            } else {
              # For custom pairs, use the comparison string match
              comp_str <- paste(pair_cond2, "vs", pair_cond1)
              if (cond == comp_str) {
                # IMPORTANT: Check if pathway is in all_significant_paths_full for this comparison
                is_in_significant_list <- path %in% all_significant_paths_full[[p_idx]]
                res <- all_results[[p_idx]]
                if (!is.null(res) && path %in% res$name && is_in_significant_list) {
                  p_val <- res$padj[res$name == path]
                  break
                }
              }
            }
          }

          # Convert p-value to stars
          if (!is.na(p_val)) {
            if (p_val < 0.001) {
              heatmap_data$sig_label[i] <- "***"
            } else if (p_val < 0.01) {
              heatmap_data$sig_label[i] <- "**"
            } else if (p_val < 0.05) {
              heatmap_data$sig_label[i] <- "*"
            }
          }
        }

        # Add text labels to heatmap
        heatmap <- heatmap +
          ggplot2::geom_text(
            ggplot2::aes(label = sig_label),
            color = "black",
            size = text_size / 3.5,
            fontface = "bold"
          )
      }
    } else {
      heatmap <- NULL
    }

    # Create a slope graph to visualize trends across conditions
    # Note: Slope plot works best for all_vs_ref method
    # For custom_pairs, we'll create a simpler version
    slope_plot <- NULL

    if (comparison_method == "all_vs_ref") {
      # Get all significant pathways from any comparison
      significant_pathways <- unique(plot_data$pathway[plot_data$significant == TRUE])

      # Get data for non-reference conditions for significant pathways
      slope_data <- plot_data[plot_data$group != ref_name &
                                plot_data$pathway %in% significant_pathways, ]

      # Add position column for each condition
      conditions <- unique(slope_data$condition)
      for (i in 1:length(conditions)) {
        slope_data$position[slope_data$condition == conditions[i]] <- i
      }

if (nrow(slope_data) > 0) {
  # Handle pathway colors - create a color palette with enough colors
  if(is.null(color.pathway)){
    # Calculate how many colors we need
    n_needed <- length(significant_pathways)

    # Create a large diverse color palette
    if(exists("global_colors")) {
      # Use global_colors as base palette
      base_palette <- global_colors
    } else {
      # Create a colorful base palette with 20 distinct colors
      base_palette <- c(
        # Colorblind-friendly palette
        "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
        "#D55E00", "#CC79A7", "#000000", "#0000FF", "#FF0000",
        # Additional distinct colors
        "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
        "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD"
      )
    }

    # If we need more colors than in the base palette, use color interpolation
    if(n_needed > length(base_palette)) {
      # Initialize full color palette
      full_palette <- base_palette

      # Use colorRampPalette to generate additional colors
      additional_colors <- grDevices::colorRampPalette(base_palette)(n_needed - length(base_palette))

      # Combine base and additional colors
      full_palette <- c(base_palette, additional_colors)

      # Ensure the colors are maximally distinct by reordering
      # This spreads similar colors further apart
      if(n_needed > 30) {
        reordered_indices <- seq(1, n_needed, by = 2)
        reordered_indices <- c(reordered_indices, seq(2, n_needed, by = 2))
        reordered_indices <- reordered_indices[1:n_needed]
        full_palette <- full_palette[reordered_indices]
      }

      color.pathway <- full_palette[1:n_needed]
    } else {
      # If we have enough colors in the base palette, just use them
      color.pathway <- base_palette[1:n_needed]
    }

    # Assign pathway names to colors
    names(color.pathway) <- significant_pathways
  }

  # Create slope plot with custom colors
  slope_plot <- ggplot2::ggplot(slope_data,
                               ggplot2::aes(x = position, y = proportion,
                                           group = pathway, color = pathway)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = 1:length(conditions), labels = conditions) +
    ggplot2::scale_color_manual(values = color.pathway) +  # Use custom color palette
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = label_angle, hjust = 1),
      text = ggplot2::element_text(size = text_size)
    ) +
    ggplot2::labs(
      title = paste("Pathway Trends Across Conditions (vs", ref_name, ")"),
      x = "Condition",
      y = "Proportion"
    )
} else {
  slope_plot <- NULL
}
    }  # End of if (comparison_method == "all_vs_ref")

    # Return all plots
    return(list(
      barplot = barplot,
      heatmap = heatmap,
      slope_plot = slope_plot
    ))
  }

  # Combine all significant pathways
  all_paths <- unique(unlist(all_significant_paths_full))

  # Prepare combined data for all pairs
  # Initialize combined_results with NULL instead of empty data frame
  combined_results <- NULL

  for (p in 1:length(comparison_pairs)) {
    if (is.null(all_results[[p]])) {
      next
    }

    pair <- comparison_pairs[[p]]
    idx1 <- comparison[pair[1]]
    idx2 <- comparison[pair[2]]

    # Add comparison info to results
    pair_result <- all_results[[p]]
    pair_result$comparison <- paste(condition_names[idx2], "vs", condition_names[idx1])

    # Combine results using rbind.fill which handles different column structures
    if(is.null(combined_results)) {
      combined_results <- pair_result
    } else {
      combined_results <- plyr::rbind.fill(combined_results, pair_result)
    }
  }

  # If still empty, create empty data frame
  if(is.null(combined_results)) {
    combined_results <- data.frame()
  }

  # Create a combined plot if requested and we have multiple pairs
  combined_plot <- NULL
  if (show_all_pairwise && length(all_plots) > 1) {
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      valid_plots <- all_plots[!sapply(all_plots, is.null)]

      if (length(valid_plots) > 0) {
        # Create a combined plot
        combined_plot <- do.call(
          gridExtra::grid.arrange,
          c(valid_plots, list(ncol = min(length(valid_plots), 2)))
        )
      }
    }
  }

  # Create combined visualization for all vs reference comparisons OR custom_pairs
  combined_comparison_plots <- NULL
  if (comparison_method %in% c("all_vs_ref", "custom_pairs")) {
    # Only create the combined plots if the specific visualization type is requested
    if (show_comparison_barplot || show_comparison_heatmap || show_comparison_slope ||
        comparison_plot_type %in% c("barplot", "heatmap", "slope")) {
      combined_comparison_plots <- create_combined_comparison_plot(
        all_significant_paths_full = all_significant_paths_full,
        all_results = all_results,
        all_df_combined = all_df_combined,
        comparison_pairs = comparison_pairs,
        condition_names = condition_names,
        comparison = comparison,
        reference = reference,
        color.use = color.use,
        text_size = text_size,
        label_angle = label_angle,
        heatmap_colors = heatmap_colors,
        border_color = border_color,
        pThresh = pThresh,
        comparison_method = comparison_method,
        color.pathway = color.pathway,
        show.all = show.all
      )
    }
  }

  # Save plots if requested
  if (save_plot) {
    # Determine which plot to save based on comparison_plot_type
    plot_to_save <- NULL

    if (comparison_method == "all_vs_ref" && !is.null(combined_comparison_plots)) {
      if (comparison_plot_type == "barplot" && !is.null(combined_comparison_plots$barplot)) {
        plot_to_save <- combined_comparison_plots$barplot
      } else if (comparison_plot_type == "heatmap" && !is.null(combined_comparison_plots$heatmap)) {
        plot_to_save <- combined_comparison_plots$heatmap
      } else if (comparison_plot_type == "slope" && !is.null(combined_comparison_plots$slope_plot)) {
        plot_to_save <- combined_comparison_plots$slope_plot
      }
    }

    # If no combined plot, use individual plots
    if (is.null(plot_to_save)) {
      if (is.null(combined_plot) && length(all_plots) == 1) {
        # Save single plot
        if (!is.null(all_plots[[1]])) {
          tryCatch({
            ggplot2::ggsave(
              filename = save_name,
              plot = all_plots[[1]],
              width = save_width,
              height = save_height
            )
            message(paste("Plot saved to", save_name))
          }, error = function(e) {
            warning(paste("Failed to save plot:", e$message))
          })
        }
      } else if (!is.null(combined_plot)) {
        # Save combined plot
        tryCatch({
          pdf(save_name, width = save_width, height = save_height)
          print(combined_plot)
          dev.off()
          message(paste("Combined plot saved to", save_name))
        }, error = function(e) {
          warning(paste("Failed to save combined plot:", e$message))
        })
      }
    } else {
      # Save the selected comparison plot
      tryCatch({
        ggplot2::ggsave(
          filename = save_name,
          plot = plot_to_save,
          width = save_width,
          height = save_height
        )
        message(paste("Plot saved to", save_name))
      }, error = function(e) {
        warning(paste("Failed to save plot:", e$message))
      })
    }
  }

  # Prepare output
  results <- list(
    plots = all_plots,
    heatmaps = all_heatmaps,
    combined_plot = combined_plot,
    data = combined_results,
    all_significant_paths = all_significant_paths,
    all_significant_paths_full = all_significant_paths_full,
    params = list(
      measure = measure,
      pThresh = pThresh,
      tol = tol
    )
  )

  # Add combined comparison plots if available
  if (!is.null(combined_comparison_plots)) {
    # Only add the plots that were requested or the default one
    if (show_comparison_barplot || comparison_plot_type == "barplot") {
      results$comparison_barplot <- combined_comparison_plots$barplot
    }
    if (show_comparison_heatmap || comparison_plot_type == "heatmap") {
      results$comparison_heatmap <- combined_comparison_plots$heatmap
    }
    if (show_comparison_slope || comparison_plot_type == "slope") {
      results$comparison_slope_plot <- combined_comparison_plots$slope_plot
    }
  }

  # Return top paths if requested
  if (return_top_paths && !is.null(top.n) && top.n > 0 && length(all_paths) > 0) {
    # Count pathway occurrences across all comparisons
    path_counts <- table(unlist(all_significant_paths_full))
    # Sort by frequency
    path_counts <- sort(path_counts, decreasing = TRUE)
    # Take top n
    top_paths <- names(path_counts)[1:min(top.n, length(path_counts))]
    results$top_paths <- top_paths
  } else {
    results$top_paths <- all_paths
  }

  return(results)
}
