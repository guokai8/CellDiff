#' Rank and Visualize Significant Signaling Pathways
#'
#' This function ranks signaling pathways based on their differential activity
#' between two groups and visualizes only the significant ones.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (length 2)
#' @param measure Measure type, either "weight" or "count"
#' @param slot.name Slot name containing the network information
#' @param color.use Vector of colors for the groups
#' @param pThresh P-value threshold for significant differences
#' @param tol Tolerance for fold change significance
#' @param top.n Number of top pathways to show
#' @param show.pval Whether to show p-values on the plot
#' @param show.heatmap Whether to generate a heatmap
#' @param title Title for the plot
#' @param sources.use Optional vector of source cell types to filter
#' @param targets.use Optional vector of target cell types to filter
#' @param return.data Logical, whether to return full results list (TRUE) or just the plot (FALSE, default)
#'
#' @return By default, returns the ggplot object. If return.data=TRUE, returns a list containing the ggplot object, heatmap, data table, and significant pathways
#' @export
rankDiff <- function(object.list, comparison = c(1, 2), measure = "weight",
                     slot.name = "netP", color.use = NULL,
                     pThresh = 0.05, tol = 0.05, top.n = NULL,
                     show.pval = TRUE, show.heatmap = TRUE, title = NULL,
                     sources.use = NULL, targets.use = NULL,
                     return.data = FALSE) {

  if (length(comparison) != 2) {
    stop("Please provide exactly 2 indices for comparison")
  }

  # Basic validation of CellChat objects
  if (!is.list(object.list) || length(object.list) < max(comparison)) {
    stop("object.list must be a list containing at least ", max(comparison), " CellChat objects")
  }

  # Process the data for comparison
  prob.list <- list()
  pSum <- list()
  pSum.original <- list()
  pair.name <- list()
  idx <- list()
  pSum.original.all <- c()
  object.names.comparison <- names(object.list)[comparison]

  for (i in 1:length(comparison)) {
    object.data <- object.list[[comparison[i]]]@netP

    prob <- object.data$prob
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
      warning("No inferred communications for object ", i)
      return(NULL)
    }

    pSum.original[[i]] <- apply(prob, 3, sum)

    if (measure == "weight") {
      pSum[[i]] <- -1/log(pSum.original[[i]] + 1e-10)
      pSum[[i]][is.na(pSum[[i]])] <- 0
      idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)

      if (length(idx[[i]]) > 0) {
        pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
      }
    } else if (measure == "count") {
      pSum[[i]] <- pSum.original[[i]]
    }

    pair.name[[i]] <- names(pSum.original[[i]])
  }

  # Replace values with max pSum for infinity cases - maintaining original logic
  if (measure == "weight" && length(unlist(idx)) > 0) {
    values.assign <- seq(max(unlist(pSum), na.rm = TRUE) * 1.1,
                         max(unlist(pSum), na.rm = TRUE) * 1.5,
                         length.out = length(unlist(idx)))
    position <- sort(pSum.original.all, index.return = TRUE)$ix

    for (i in 1:length(comparison)) {
      if (length(idx[[i]]) > 0) {
        if (i == 1) {
          pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), position)]
        } else {
          pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i-1])) + 1:length(idx[[i]]), position)]
        }
      }
    }
  }

  # Extract all pathway names
  pair.name.all <- as.character(unique(unlist(pair.name)))

  if (length(pair.name.all) == 0) {
    warning("No pathway names found")
    return(NULL)
  }

  # Create data frames with pathway contributions
  df <- list()
  for (i in 1:length(comparison)) {
    df[[i]] <- data.frame(
      name = pair.name.all,
      contribution = 0,
      contribution.scaled = 0,
      group = object.names.comparison[i],
      stringsAsFactors = FALSE
    )

    matching_names <- intersect(pair.name[[i]], pair.name.all)

    if (length(matching_names) > 0) {
      df[[i]][match(matching_names, df[[i]]$name), "contribution.scaled"] <- pSum[[i]][matching_names]
      df[[i]][match(matching_names, df[[i]]$name), "contribution"] <- pSum.original[[i]][matching_names]
    }
  }

  # Combine data frames
  df_combined <- do.call(rbind, df)

  # Filter out pathways with zero contribution
  pathway_sums <- tapply(df_combined$contribution, df_combined$name, sum)
  valid_pathways <- names(pathway_sums[pathway_sums > 0])

  if (length(valid_pathways) == 0) {
    warning("No pathways with non-zero contribution")
    return(NULL)
  }

  df_combined <- df_combined[df_combined$name %in% valid_pathways, ]

  # Set factor levels for proper ordering
  df_combined$group <- factor(df_combined$group, levels = object.names.comparison)

  # Calculate aggregate values by group
  df_aggregated <- aggregate(
    cbind(contribution, contribution.scaled) ~ group + name,
    data = df_combined,
    FUN = mean
  )

  # Calculate relative changes
  wide_df <- reshape(
    df_aggregated[, c("group", "name", "contribution")],
    timevar = "group",
    idvar = "name",
    direction = "wide"
  )

  # Rename columns to match conditions
  colnames(wide_df) <- gsub("^contribution\\.", "", colnames(wide_df))

  # Calculate relative contribution
  wide_df$contribution.relative <- wide_df[[object.names.comparison[2]]] /
    (wide_df[[object.names.comparison[1]]] + 1e-10)

  # Calculate p-values with Wilcoxon test - using the exact same approach as in the original code
  p_values <- numeric(length(valid_pathways))
  names(p_values) <- valid_pathways

  for (i in 1:length(valid_pathways)) {
    pathway <- valid_pathways[i]

    # Extract values from each condition separately to handle different dimensions
    prob_values_list <- list()

    for (j in 1:length(comparison)) {
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
    prob_values <- matrix(NA, nrow = max_len, ncol = length(comparison))

    # Fill in the values
    for (j in 1:length(comparison)) {
      if (length(prob_values_list[[j]]) > 0) {
        prob_values[1:length(prob_values_list[[j]]), j] <- prob_values_list[[j]]
      }
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
  wide_df$padj <- p_values[wide_df$name]  # No adjustment

  # Sort by contribution (keeping original ordering for significance testing)
  wide_df <- wide_df[order(wide_df[[object.names.comparison[1]]], decreasing = TRUE), ]

  # Identify significant pathways
  is_significant <- (wide_df$contribution.relative < (1 - tol) |
                       wide_df$contribution.relative > (1 + tol)) &
    wide_df$padj < pThresh

  significant_paths <- wide_df$name[is_significant]

  if (length(significant_paths) == 0) {
    warning("No significant pathways found with the current thresholds")
    return(list(
      plot = NULL,
      heatmap = NULL,
      data = wide_df,
      significant_paths = character(0)
    ))
  }

  # Limit to top pathways if specified
  if (!is.null(top.n) && top.n < length(significant_paths)) {
    # Keep the top.n significant pathways in order of contribution
    significant_paths <- significant_paths[1:top.n]
  }

  # Filter data for plotting
  plot_data <- df_combined[df_combined$name %in% significant_paths, ]

  # Set colors for groups
  if (is.null(color.use)) {
    color.use <- global_colors[1:2]
  }
  names(color.use) <- object.names.comparison

  # Create text colors based on significance and direction
  colors.text <- rep("black", length(significant_paths))
  names(colors.text) <- significant_paths

  for (i in 1:length(significant_paths)) {
    path <- significant_paths[i]
    row_idx <- which(wide_df$name == path)

    if (length(row_idx) > 0) {
      ratio <- wide_df$contribution.relative[row_idx]
      p_val <- wide_df$padj[row_idx]

      if ((ratio < (1 - tol)) && (p_val < pThresh)) {
        colors.text[i] <- color.use[1]  # First group is stronger
      } else if ((ratio > (1 + tol)) && (p_val < pThresh)) {
        colors.text[i] <- color.use[2]  # Second group is stronger
      }
    }
  }
  # Get unique pathway names
  pathway_names <- unique(plot_data$name)

  # Create empty vectors to store the fold changes
  fc_values <- numeric(length(pathway_names))
  names(fc_values) <- pathway_names
  group_names <- unique(as.character(plot_data$group))

  # Calculate fold changes manually
  for (pathway in pathway_names) {
    # Get values for each group
    nl_value <- plot_data$contribution.scaled[plot_data$name == pathway & plot_data$group == group_names[1]]
    ls_value <- plot_data$contribution.scaled[plot_data$name == pathway & plot_data$group == group_names[2]]

    # Calculate fold change if both values exist and are non-zero
    if (length(nl_value) > 0 && length(ls_value) > 0) {
      # Add small epsilon to avoid division by zero
      fc_values[pathway] <- (ls_value + 1e-10) / (nl_value + 1e-10)
    } else {
      # If a pathway is missing in one group, set a default value
      fc_values[pathway] <- 1  # No change
    }
  }

  # Convert fold changes to log2 scale for better comparison
  log2_fc <- fc_values

  # Order by absolute log2 fold change (most changed first)
  ordered_names <- names(log2_fc)[order(abs(log2_fc), decreasing = TRUE)]

  # Convert name to a factor with custom order
  plot_data$name <- factor(plot_data$name, levels = ordered_names)

  min_nonzero <- min(plot_data$contribution[plot_data$contribution > 0])
  # Compute scale factor so smallest non-zero becomes >= 1e-6
  scale_factor <- 1e-6 / min_nonzero
  # Create a scaled value ONLY for plotting
  plot_data$contribution <- plot_data$contribution * scale_factor
  # Replace exact zeros with tiny value to prevent rendering artifacts
  plot_data$contribution[plot_data$contribution == 0] <- 1e-20

  # Now create the plot with your ordered factor
  gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x = name, y = contribution, fill = group)) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::xlab("") +
    ggplot2::ylab("Relative information flow") +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::coord_flip() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(colour = colors.text[levels(plot_data$name)], size = 10),
      legend.position = "bottom"
    ) +
    ggplot2::scale_fill_manual(values = color.use)
  # Convert name to a factor with your custom order

  # Add title if provided
  if (!is.null(title)) {
    gg <- gg + ggplot2::ggtitle(as.character(title)[1])
  }

  # Add p-values to the plot if requested
  if (show.pval && length(significant_paths) > 0) {
    # Prepare p-value annotations
    pval_data <- wide_df[wide_df$name %in% significant_paths, ]

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
        group = object.names.comparison[1]  # Doesn't matter which group, just for position
      )

      gg <- gg +
        ggplot2::geom_text(
          data = label_data,
          ggplot2::aes(label = label),
          position = ggplot2::position_dodge(width = 1),
          hjust = 1
        )
    }
  }

  # Create heatmap if requested
  heatmap <- NULL
  if (show.heatmap && length(significant_paths) >= 2) {  # Need at least 2 for clustering
    # Prepare data for heatmap
    heatmap_data <- reshape(
      df_combined[df_combined$name %in% significant_paths, c("group", "name", "contribution")],
      timevar = "group",
      idvar = "name",
      direction = "wide"
    )

    # Set rownames for the heatmap
    rownames(heatmap_data) <- heatmap_data$name
    heatmap_data$name <- NULL

    # Rename columns
    colnames(heatmap_data) <- gsub("^contribution\\.", "", colnames(heatmap_data))

    # Create annotation for groups
    annotation_col <- data.frame(
      Group = factor(colnames(heatmap_data))
    )
    rownames(annotation_col) <- colnames(heatmap_data)

    # Create heatmap - safely check if pheatmap is available
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      tryCatch({
        heatmap <- pheatmap::pheatmap(
          heatmap_data,
          cluster_cols = FALSE,
          scale = "row",
          annotation_col = annotation_col,
          annotation_colors = list(Group = color.use),
          border_color = "white",
          color = colorRampPalette(c("deepskyblue", "white", "darkorange"))(128),
          main = "Significant Pathway Differences",
          fontsize_row = 10,
          fontsize_col = 10,
          silent = TRUE
        )
      }, error = function(e) {
        warning("Error creating heatmap: ", e$message)
        heatmap <- NULL
      })
    }
  }

  # Prepare output
  if (return.data) {
    # Return full results list
    results <- list(
      plot = gg,
      heatmap = heatmap,
      data = wide_df,
      significant_paths = significant_paths,
      plot_data=plot_data
    )
    return(results)
  } else {
    # Return just the plot by default
    return(gg)
  }
}
