#' Compare Ligand-Receptor Pair Contributions Between Two Groups
#'
#' This function compares the contributions of ligand-receptor pairs to a specific signaling pathway
#' between two groups. It allows visualization of the relative contributions and their differences.
#'
#' @param object.list A list containing exactly two CellChat objects to compare
#' @param signaling Character string, the name of the signaling pathway to be analyzed
#' @param signaling.name Character string, custom name for the signaling pathway (if NULL, uses the value of signaling)
#' @param sources.use Character vector or numeric indices, source cell groups to include in analysis
#'        (if NULL, all source cell groups will be included)
#' @param targets.use Character vector or numeric indices, target cell groups to include in analysis
#'        (if NULL, all target cell groups will be included)
#' @param width Numeric, width of the bars (default: 0.1)
#' @param vertex.receiver Numeric vector, indices of receiver cells for hierarchical analysis
#'        (if NULL, a standard plot without hierarchy is generated)
#' @param thresh Numeric, p-value threshold for including interactions (default: 0.05)
#' @param group.names Character vector of length 2, custom names for the two groups being compared
#'        (default: NULL, which will use the names of the list elements in object.list if available,
#'        otherwise defaults to c("Group1", "Group2"))
#' @param return.data Logical, whether to return the data along with the plot (default: FALSE)
#' @param x.rotation Numeric, rotation angle for x-axis labels in degrees (default: 0)
#' @param title Character string, title of the plot (default: "Comparison of L-R pair contributions")
#' @param font.size Numeric, font size for axis text (default: 10)
#' @param font.size.title Numeric, font size for plot title (default: 10)
#' @param show.difference Logical, whether to show difference between the two groups (default: TRUE)
#' @param color.use Character vector of length 2, custom colors for the two groups
#'        (default: c("#4682B4", "#B4464B"))
#' @param stack.method Character string, visualization method, either "side-by-side" or "stacked"
#'        (default: "side-by-side")
#'
#' @return If return.data = FALSE, returns a ggplot object. If return.data = TRUE, returns a list with two elements:
#' \itemize{
#'   \item LR.contribution: A data frame containing the contribution values for each L-R pair in each group
#'   \item gg.obj: The ggplot object
#' }
#'
#' @details
#' The function analyzes the contribution of each ligand-receptor pair to a specific signaling pathway
#' and compares these contributions between two groups. When show.difference = TRUE, it also calculates
#' the difference in contribution between Group2 and Group1 for each L-R pair.
#'
#' The stack.method parameter controls how the data is visualized:
#' \itemize{
#'   \item "side-by-side": Creates a bar plot where the groups are displayed side by side for direct comparison
#'   \item "stacked": Creates a stacked bar plot showing the total contribution, and a separate plot showing the differences
#' }
#'
#' @examples
#' # Basic usage with named list
#' cellchat.list <- list(Normal = cellchat1, Tumor = cellchat2)
#' ContriDiff(
#'   object.list = cellchat.list,
#'   signaling = "TGFb"
#'   # group.names will automatically use "Normal" and "Tumor" from the list names
#' )
#'
#' # Or with explicit group names
#' ContriDiff(
#'   object.list = cellchat.list,
#'   signaling = "TGFb",
#'   group.names = c("Healthy", "Disease")  # Override list names
#' )
#'
#' # Show differences and customize
#' ContriDiff(
#'   object.list = cellchat.list,
#'   signaling = "WNT",
#'   show.difference = TRUE,
#'   stack.method = "stacked",
#'   color.use = c("#2E86C1", "#C0392B")
#' )
#'
#' # Return data for further analysis
#' result <- ContriDiff(
#'   object.list = cellchat.list,
#'   signaling = "VEGF",
#'   return.data = TRUE
#' )
#' lr_data <- result$LR.contribution
#' top_contributors <- head(lr_data[order(lr_data$contribution, decreasing = TRUE), ], 5)
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme_classic theme element_text xlab ylab coord_flip ggtitle
#' @importFrom cowplot ggdraw draw_label plot_grid
#'
#' @export
ContriDiff <- function(object.list, signaling, signaling.name = NULL, sources.use = NULL,
                       targets.use = NULL, width = 0.1, vertex.receiver = NULL,
                       thresh = 0.05, group.names = NULL, return.data = FALSE,
                       x.rotation = 0, title = "Comparison of L-R pair contributions",
                       font.size = 10, font.size.title = 10, show.difference = TRUE,
                       color.use = c("deepskyblue", "#B4464B"), stack.method = "side-by-side")
{
  # Check if object.list is a list containing two CellChat objects
  if (!is.list(object.list) || length(object.list) != 2) {
    stop("object.list must be a list containing two CellChat objects")
  }

  # If group.names is NULL, try to use the names from the list
  if (is.null(group.names)) {
    if (!is.null(names(object.list)) && all(names(object.list) != "")) {
      group.names <- names(object.list)
    } else {
      group.names <- c("Group1", "Group2")
    }
  }

  # Ensure group.names has exactly 2 elements
  if (length(group.names) != 2) {
    warning("group.names must have exactly 2 elements. Using default names.")
    group.names <- c("Group1", "Group2")
  }

  # Create a list to store data for each group
  df.list <- list()
  lr.pairs.all <- c()

  # Process data for each group
  for (i in 1:2) {
    object <- object.list[[i]]

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
    if (max(prob) > min(prob)) {
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

  # If no valid data in both groups, return empty plot
  if (length(lr.pairs.all) == 0) {
    warning("No valid ligand-receptor pairs found in both groups")
    empty_plot <- ggplot() +
      theme_void() +
      annotate("text", x = 0, y = 0, label = "No valid ligand-receptor pairs detected", size = 5)

    if (return.data) {
      return(list(
        LR.contribution = data.frame(name = character(0), contribution = numeric(0), group = character(0)),
        gg.obj = empty_plot
      ))
    } else {
      return(empty_plot)
    }
  }

  # Ensure dataframes contain all L-R pairs
  for (i in 1:2) {
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

  # Combine data from both groups
  df.combined <- do.call(rbind, df.list)

  # Calculate differences
  if (show.difference && length(unique(df.combined$group)) == 2) {
    df.diff <- data.frame()
    for (pair in unique(df.combined$name)) {
      contrib1 <- df.list[[1]][df.list[[1]]$name == pair, "contribution"]
      contrib2 <- df.list[[2]][df.list[[2]]$name == pair, "contribution"]

      if (length(contrib1) > 0 && length(contrib2) > 0) {
        diff_val <- contrib2 - contrib1

        diff_row <- data.frame(
          name = pair,
          contribution = diff_val,
          group = "Difference"
        )

        df.diff <- rbind(df.diff, diff_row)
      }
    }

    # Add difference data if available
    if (nrow(df.diff) > 0) {
      df.combined <- rbind(df.combined, df.diff)
    }
  }

  # Check if dataframe is not empty
  if (nrow(df.combined) == 0) {
    warning("No valid ligand-receptor contribution data")
    empty_plot <- ggplot() +
      theme_void() +
      annotate("text", x = 0, y = 0, label = "No valid ligand-receptor contribution data", size = 5)

    if (return.data) {
      return(list(
        LR.contribution = data.frame(name = character(0), contribution = numeric(0), group = character(0)),
        gg.obj = empty_plot
      ))
    } else {
      return(empty_plot)
    }
  }

  # Sort L-R pairs by contribution
  lr.order <- tapply(df.combined$contribution, df.combined$name, sum)
  lr.order <- names(sort(lr.order, decreasing = TRUE))
  df.combined$name <- factor(df.combined$name, levels = lr.order)

  # Create comparison plot
  if (stack.method == "side-by-side") {
    # Horizontal bar plot, side-by-side
    gg <- ggplot(df.combined, aes(x = name, y = contribution, fill = group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
      scale_fill_manual(values = c(color.use, "grey50")[1:length(unique(df.combined$group))]) +
      theme_classic() +
      theme(
        axis.text.y = element_text(angle = 0, hjust = 1, size = font.size, colour = "black"),
        axis.text = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.x = element_text(angle = x.rotation, hjust = 1, size = font.size),
        legend.position = "top",
        legend.title = element_blank()
      ) +
      xlab("") +
      ylab("Relative contribution") +
      coord_flip()
  } else if (stack.method == "stacked") {
    # Remove difference data for stacked plot
    df.stacked <- df.combined[df.combined$group != "Difference", ]

    # Horizontal stacked bar plot
    gg <- ggplot(df.stacked, aes(x = name, y = contribution, fill = group)) +
      geom_bar(stat = "identity", position = "stack", width = 0.6) +
      scale_fill_manual(values = color.use[1:length(unique(df.stacked$group))]) +
      theme_classic() +
      theme(
        axis.text.y = element_text(angle = 0, hjust = 1, size = font.size, colour = "black"),
        axis.text = element_text(size = font.size),
        axis.title.y = element_text(size = font.size),
        axis.text.x = element_text(angle = x.rotation, hjust = 1, size = font.size),
        legend.position = "top",
        legend.title = element_blank()
      ) +
      xlab("") +
      ylab("Relative contribution") +
      coord_flip()

    # If difference data exists, create difference plot
    if (show.difference && "Difference" %in% df.combined$group) {
      df.diff <- df.combined[df.combined$group == "Difference", ]

      # Create difference plot
      gg.diff <- ggplot(df.diff, aes(x = name, y = contribution, fill = contribution > 0)) +
        geom_bar(stat = "identity", width = 0.6) +
        scale_fill_manual(values = c("deepskyblue", "darkorange"),
                          labels = c("Decreased", "Increased"),
                          name = "Change") +
        theme_classic() +
        theme(
          axis.text.y = element_text(angle = 0, hjust = 1, size = font.size, colour = "black"),
          axis.text = element_text(size = font.size),
          axis.title.y = element_text(size = font.size),
          axis.text.x = element_text(angle = x.rotation, hjust = 1, size = font.size),
          legend.position = "top"
        ) +
        xlab("") +
        ylab("Contribution difference") +
        coord_flip()

      # Combine plots
      title <- cowplot::ggdraw() + cowplot::draw_label(title, fontface = "bold", size = font.size.title)
      gg.combined <- cowplot::plot_grid(gg, gg.diff, nrow = 1, labels = c("A", "B"))
      gg <- cowplot::plot_grid(title, gg.combined, ncol = 1, rel_heights = c(0.1, 1))
    }
  }

  # Add title if not already added for stacked difference plots
  if (!is.null(title) && !(stack.method == "stacked" && show.difference && "Difference" %in% df.combined$group)) {
    gg <- gg + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = font.size.title))
  }

  # Return data or plot
  if (return.data) {
    df.filtered <- subset(df.combined, contribution > 0)
    return(list(LR.contribution = df.filtered, gg.obj = gg))
  }
  else {
    return(gg)
  }
}
