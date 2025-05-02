#' Compare Ligand-Receptor Pair Contributions Across Multiple Groups
#'
#' This function compares the contributions of ligand-receptor pairs to specific signaling pathways
#' across multiple groups. It allows visualization of the relative contributions and their differences
#' in an integrated view.
#'
#' @param object.list A list containing CellChat objects to compare
#' @param signaling Character string, the name of the signaling pathway to be analyzed
#' @param comparison Optional vector of indices in object.list to include in the comparison.
#'        If NULL (default), all objects will be compared.
#' @param signaling.name Character string, custom name for the signaling pathway (if NULL, uses the value of signaling)
#' @param sources.use Character vector or numeric indices, source cell groups to include in analysis
#'        (if NULL, all source cell groups will be included)
#' @param targets.use Character vector or numeric indices, target cell groups to include in analysis
#'        (if NULL, all target cell groups will be included)
#' @param width Numeric, width of the bars (default: 0.1)
#' @param thresh Numeric, p-value threshold for including interactions (default: 0.05)
#' @param group.names Character vector, custom names for the groups being compared
#'        (default: NULL, which will use the names of the list elements in object.list if available)
#' @param return.data Logical, whether to return the data along with the plot (default: FALSE)
#' @param x.rotation Numeric, rotation angle for x-axis labels in degrees (default: 0)
#' @param title Character string, title of the plot (default: "Comparison of L-R pair contributions")
#' @param font.size Numeric, font size for axis text (default: 10)
#' @param font.size.title Numeric, font size for plot title (default: 10)
#' @param color.use Character vector, custom colors for the groups (default: uses global_colors)
#' @param stack.method Character string, visualization method, either "side-by-side" or "stacked"
#'        (default: "side-by-side")
#' @param reference Integer or character, reference group for comparison (default: 1)
#' @param show.reference Logical, whether to show the reference group in visualizations (default: TRUE)
#' @param normalize Logical, whether to normalize contributions within each group (default: TRUE)
#' @param show.heatmap Logical, whether to show a heatmap visualization in addition to the bar plot (default: FALSE)
#' @param filter.min.contrib Numeric, minimum contribution to include an L-R pair in the visualization (default: 0)
#' @param top.n Integer, limit visualization to top N L-R pairs by contribution (default: NULL, shows all)
#'
#' @return If return.data = FALSE, returns a list of ggplot objects. If return.data = TRUE, returns a list with:
#' \itemize{
#'   \item LR.contribution: A data frame containing the contribution values for each L-R pair in each group
#'   \item plots: A list of ggplot objects (bar plot and heatmap if requested)
#' }
#'
#' @details
#' The function analyzes the contribution of each ligand-receptor pair to a specific signaling pathway
#' and compares these contributions across multiple groups. When reference is specified, differences
#' are calculated relative to the reference group.
#'
#' The stack.method parameter controls how the data is visualized:
#' \itemize{
#'   \item "side-by-side": Creates a bar plot where the groups are displayed side by side for direct comparison
#'   \item "stacked": Creates a stacked bar plot showing the total contribution across groups
#' }
#'
#' @examples
#' # Basic usage with named list
#' cellchat.list <- list(Normal = cellchat1, Mild = cellchat2, Severe = cellchat3)
#' ContriDiffM(
#'   object.list = cellchat.list,
#'   signaling = "TGFb"
#' )
#'
#' # Compare with custom group names and colors
#' ContriDiffM(
#'   object.list = cellchat.list,
#'   signaling = "TGFb",
#'   group.names = c("Healthy", "Early Stage", "Late Stage"),
#'   color.use = c("blue", "green", "red")
#' )
#'
#' # Specify a reference group and show heatmap
#' ContriDiffM(
#'   object.list = cellchat.list,
#'   signaling = "WNT",
#'   reference = "Normal",
#'   show.heatmap = TRUE
#' )
#'
#' # Return data for further analysis
#' result <- ContriDiffM(
#'   object.list = cellchat.list,
#'   signaling = "VEGF",
#'   return.data = TRUE
#' )
#' lr_data <- result$LR.contribution
#' top_contributors <- head(lr_data[order(lr_data$contribution, decreasing = TRUE), ], 5)
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme_classic theme element_text xlab ylab coord_flip ggtitle
#' @importFrom cowplot ggdraw draw_label plot_grid
#' @importFrom reshape2 melt dcast
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#'
#' @export
ContriDiffM <- function(object.list, signaling, comparison = NULL, signaling.name = NULL,
                        sources.use = NULL, targets.use = NULL, width = 0.1,
                        thresh = 0.05, group.names = NULL, return.data = FALSE,
                        x.rotation = 0, title = "Comparison of L-R pair contributions",
                        font.size = 10, font.size.title = 10,
                        color.use = NULL, stack.method = "side-by-side",
                        reference = 1, show.reference = TRUE,
                        normalize = TRUE, show.heatmap = FALSE,
                        filter.min.contrib = 0, top.n = NULL) {

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

  # Check if object.list is valid
  if (!is.list(object.list) || length(object.list) < max(comparison)) {
    stop("object.list must be a list containing at least as many elements as the maximum index in comparison")
  }

  # Process the reference parameter
  if (is.character(reference) && length(reference) == 1) {
    if (reference %in% names(object.list)) {
      reference <- match(reference, names(object.list))
    } else {
      stop("Reference name not found in object.list")
    }
  }

  if (!reference %in% comparison) {
    warning("Reference index not in comparison indices. Using first index in comparison.")
    reference <- comparison[1]
  }

  # If group.names is NULL, try to use the names from the list
  if (is.null(group.names)) {
    if (!is.null(names(object.list)) && all(names(object.list)[comparison] != "")) {
      group.names <- names(object.list)[comparison]
    } else {
      group.names <- paste0("Group", comparison)
    }
  } else if (length(group.names) != length(comparison)) {
    warning("group.names length doesn't match comparison length. Using default names.")
    group.names <- paste0("Group", comparison)
  }

  # Create a list to store data for each group
  df.list <- list()
  lr.pairs.all <- c()

  # Process data for each group
  for (i in 1:length(comparison)) {
    idx <- comparison[i]
    object <- object.list[[idx]]

    # Search for ligand-receptor pairs
    pairLR <- CellChat::searchPair(signaling = signaling, pairLR.use = object@LR$LRsig,
                                   key = "pathway_name", matching.exact = TRUE, pair.only = TRUE)

    # Skip this group if no L-R pairs found
    if (nrow(pairLR) == 0) {
      warning(paste("No ligand-receptor pairs found for", signaling, "pathway in", group.names[i]))
      df.list[[i]] <- data.frame(name = character(0), contribution = numeric(0), group = character(0))
      next
    }

    pair.name.use = object@DB$interaction[rownames(pairLR), "interaction_name_2"]

    if (is.null(signaling.name)) {
      signaling.name <- signaling
    }

    net <- object@net
    pairLR.use.name <- dimnames(net$prob)[[3]]
    pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
    pairLR <- pairLR[pairLR.name, , drop = FALSE]

    # Skip if no valid L-R pairs found
    if (nrow(pairLR) == 0) {
      warning(paste("No valid ligand-receptor pairs found for", signaling, "pathway in", group.names[i]))
      df.list[[i]] <- data.frame(name = character(0), contribution = numeric(0), group = character(0))
      next
    }

    prob <- net$prob
    pval <- net$pval
    prob[pval > thresh] <- 0

    # Process sources.use parameter
    if (!is.null(sources.use)) {
      if (is.character(sources.use)) {
        if (all(sources.use %in% dimnames(prob)[[1]])) {
          sources.use <- match(sources.use, dimnames(prob)[[1]])
        }
        else {
          stop("The input `sources.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), sources.use)
      prob[idx.t, , ] <- 0
    }

    # Process targets.use parameter
    if (!is.null(targets.use)) {
      if (is.character(targets.use)) {
        if (all(targets.use %in% dimnames(prob)[[2]])) {
          targets.use <- match(targets.use, dimnames(prob)[[2]])
        }
        else {
          stop("The input `targets.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), targets.use)
      prob[, idx.t, ] <- 0
    }

    # Filter valid L-R pairs
    if (length(pairLR.name) > 1) {
      pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name, drop = FALSE],
                                           3, sum) != 0]
    }
    else {
      pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name, drop = FALSE]) != 0]
    }

    # Skip if no significant L-R pairs found
    if (length(pairLR.name.use) == 0) {
      warning(paste("No significant communication of", signaling.name, "pathway in", group.names[i]))
      df.list[[i]] <- data.frame(name = character(0), contribution = numeric(0), group = character(0))
      next
    }

    # Filter valid L-R pairs and normalize probabilities
    pairLR <- pairLR[pairLR.name.use, , drop = FALSE]
    prob <- prob[, , pairLR.name.use, drop = FALSE]

    # Handle single L-R pair case
    if (length(dim(prob)) == 2) {
      prob <- replicate(1, prob, simplify = "array")
      dimnames(prob)[[3]] <- pairLR.name.use
    }

    # Normalize probabilities
    if (normalize && max(prob) > min(prob)) {
      prob <- (prob - min(prob))/(max(prob) - min(prob))
    }

    # Calculate contribution for each L-R pair
    pSum <- apply(prob, 3, sum)
    pSum.max <- sum(prob)
    if (pSum.max > 0) {
      pSum <- pSum/pSum.max
    }
    pSum[is.na(pSum)] <- 0

    # Extract L-R pair names
    pair.name <- unlist(dimnames(prob)[3])
    pair.name <- factor(pair.name, levels = unique(pair.name))

    if (!is.null(pairLR.name.use)) {
      pair.name_df <- data.frame(row.names = pairLR.name.use)
      pair.name_df$interaction_name_2 <- object@DB$interaction[pairLR.name.use, "interaction_name_2"]
      pair.name <- pair.name_df[as.character(pair.name), "interaction_name_2", drop = TRUE]
      pair.name <- factor(pair.name, levels = unique(pair.name))
    }

    # Create dataframe
    df <- data.frame(name = pair.name, contribution = pSum, group = group.names[i])
    df$name <- as.character(df$name)

    # Add to list
    df.list[[i]] <- df

    # Collect all L-R pairs
    lr.pairs.all <- c(lr.pairs.all, as.character(df$name))
  }

  # Merge data
  lr.pairs.all <- unique(lr.pairs.all)

  # If no valid data in any groups, return empty plot
  if (length(lr.pairs.all) == 0) {
    warning("No valid ligand-receptor pairs found in any group")
    empty_plot <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No valid ligand-receptor pairs detected", size = 5)

    if (return.data) {
      return(list(
        LR.contribution = data.frame(name = character(0), contribution = numeric(0), group = character(0)),
        plots = list(barplot = empty_plot)
      ))
    } else {
      return(list(barplot = empty_plot))
    }
  }

  # Ensure dataframes contain all L-R pairs
  for (i in 1:length(comparison)) {
    missing_pairs <- setdiff(lr.pairs.all, df.list[[i]]$name)
    if (length(missing_pairs) > 0) {
      df_missing <- data.frame(
        name = missing_pairs,
        contribution = rep(0, length(missing_pairs)),
        group = rep(group.names[i], length(missing_pairs))
      )
      df.list[[i]] <- rbind(df.list[[i]], df_missing)
    }
  }

  # Combine data from all groups
  df.combined <- do.call(rbind, df.list)

  # Filter out LR pairs with contributions below threshold
  if (filter.min.contrib > 0) {
    filtered_pairs <- names(which(tapply(df.combined$contribution, df.combined$name, max) >= filter.min.contrib))
    df.combined <- df.combined[df.combined$name %in% filtered_pairs, ]

    if (nrow(df.combined) == 0) {
      warning("No ligand-receptor pairs meet the minimum contribution threshold")
      empty_plot <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0, y = 0,
                          label = "No ligand-receptor pairs meet the minimum contribution threshold",
                          size = 5)

      if (return.data) {
        return(list(
          LR.contribution = data.frame(name = character(0), contribution = numeric(0), group = character(0)),
          plots = list(barplot = empty_plot)
        ))
      } else {
        return(list(barplot = empty_plot))
      }
    }
  }

  # Limit to top N L-R pairs if specified
  if (!is.null(top.n) && top.n > 0) {
    pair_sums <- tapply(df.combined$contribution, df.combined$name, sum)
    top_pairs <- names(sort(pair_sums, decreasing = TRUE)[1:min(top.n, length(pair_sums))])
    df.combined <- df.combined[df.combined$name %in% top_pairs, ]
  }

  # Calculate differences from reference
  ref_idx <- which(comparison == reference)
  if (length(ref_idx) > 0 && show.reference) {
    # Get reference group data
    ref_data <- df.list[[ref_idx]]

    # Calculate differences for each other group
    diff.list <- list()
    for (i in 1:length(comparison)) {
      if (i != ref_idx) {
        curr_data <- df.list[[i]]
        diff_data <- data.frame(
          name = lr.pairs.all,
          contribution = 0,
          group = paste0(group.names[i], " vs ", group.names[ref_idx]),
          stringsAsFactors = FALSE
        )

        # Calculate differences
        for (pair in lr.pairs.all) {
          ref_val <- ref_data$contribution[ref_data$name == pair]
          curr_val <- curr_data$contribution[curr_data$name == pair]

          if (length(ref_val) > 0 && length(curr_val) > 0) {
            diff_data$contribution[diff_data$name == pair] <- curr_val - ref_val
          }
        }

        diff.list[[i]] <- diff_data
      }
    }

    # Combine difference data
    if (length(diff.list) > 0) {
      diff.combined <- do.call(rbind, Filter(function(x) !is.null(x), diff.list))
      # Add to main data if we want to include differences in visualization
      # df.combined <- rbind(df.combined, diff.combined)
    }
  }

  # Sort L-R pairs by contribution
  lr.order <- tapply(df.combined$contribution, df.combined$name, sum)
  lr.order <- names(sort(lr.order, decreasing = TRUE))
  df.combined$name <- factor(df.combined$name, levels = lr.order)

  # Set colors for groups
  if (is.null(color.use)) {
    if (exists("global_colors") && length(global_colors) >= length(group.names)) {
      color.use <- global_colors[1:length(group.names)]
    } else {
      color.use <- scales::hue_pal()(length(group.names))
    }
  } else if (length(color.use) < length(group.names)) {
    # Extend colors if needed
    color.use <- rep(color.use, length.out = length(group.names))
  }
  names(color.use) <- group.names

  # Create the main comparison plot
  if (stack.method == "side-by-side") {
    # Horizontal bar plot, side-by-side
    gg <- ggplot2::ggplot(df.combined, ggplot2::aes(x = name, y = contribution, fill = group)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
      ggplot2::scale_fill_manual(values = color.use) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, size = font.size, colour = "black"),
        axis.text = ggplot2::element_text(size = font.size),
        axis.title.y = ggplot2::element_text(size = font.size),
        axis.text.x = ggplot2::element_text(angle = x.rotation, hjust = 1, size = font.size),
        legend.position = "top",
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("Relative contribution") +
      ggplot2::coord_flip()
  } else if (stack.method == "stacked") {
    # Horizontal stacked bar plot
    gg <- ggplot2::ggplot(df.combined, ggplot2::aes(x = name, y = contribution, fill = group)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.6) +
      ggplot2::scale_fill_manual(values = color.use) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, size = font.size, colour = "black"),
        axis.text = ggplot2::element_text(size = font.size),
        axis.title.y = ggplot2::element_text(size = font.size),
        axis.text.x = ggplot2::element_text(angle = x.rotation, hjust = 1, size = font.size),
        legend.position = "top",
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("Relative contribution") +
      ggplot2::coord_flip()
  }

  # Add title
  gg <- gg + ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = font.size.title))

  # Create heatmap visualization if requested
  heatmap <- NULL
  if (show.heatmap) {
    # Convert to wide format for heatmap
    heatmap_data <- reshape2::dcast(df.combined, name ~ group, value.var = "contribution")
    rownames(heatmap_data) <- heatmap_data$name
    heatmap_data$name <- NULL

    # Define color function
    col_fun <- circlize::colorRamp2(
      c(0, max(heatmap_data, na.rm = TRUE)/2, max(heatmap_data, na.rm = TRUE)),
      c("white", "orange", "red")
    )

    # Create heatmap
    heatmap <- ComplexHeatmap::Heatmap(
      as.matrix(heatmap_data),
      name = "Contribution",
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_gp = grid::gpar(fontsize = font.size),
      column_names_gp = grid::gpar(fontsize = font.size),
      column_title = paste(title, "- Heatmap"),
      column_title_gp = grid::gpar(fontsize = font.size.title)
    )
  }

  # Prepare output
  plots <- list(barplot = gg)
  if (show.heatmap) {
    plots$heatmap <- heatmap
  }

  # Return data or plot
  if (return.data) {
    return(list(
      LR.contribution = df.combined,
      plots = plots
    ))
  } else {
    return(plots)
  }
}
