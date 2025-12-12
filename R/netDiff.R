#' Visualize Cell-Cell Communication Differences with Chord Diagram
#'
#' This function creates a chord diagram showing the differences in cell-cell communication
#' between two conditions, highlighting changes in signaling strength. It shows all cell types
#' and includes pathways that are only present in one condition.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (length 2)
#' @param signaling Signaling pathway name
#' @param thresh Significance threshold for p-values (interactions with p-value > thresh will be set to 0)
#' @param color.use Vector of colors for increased and decreased interactions
#' @param cell_color Named vector of colors for cell types (default: NULL, uses global_colors or default palette)
#' @param title.name Title for the plot
#' @param sources.use Optional vector of source cell types to include
#' @param targets.use Optional vector of target cell types to include
#' @param remove.isolate Whether to remove isolated nodes (default: FALSE)
#' @param lab.cex Size of labels
#' @param small.gap Gap between segments
#' @param big.gap Gap between groups
#' @param directional Whether to show directional information
#' @param link.visible Whether to show links
#' @param link.border Border color for links
#' @param transparency Transparency of links
#' @param min_weight Minimum absolute weight to display (to filter out trivial differences)
#' @param return_data Logical, whether to return both the plot and data (TRUE) or just the plot (FALSE)
#' @param save_plot Logical, whether to save the plot to a file
#' @param save_name Character, filename for saving the plot
#' @param save_width Width of the saved plot in inches
#' @param save_height Height of the saved plot in inches
#' @param slot.name Character, name of the slot to extract data from (default: "net")
#' @param show_all_cell_types Logical, whether to show all cell types, even if not present in the signaling pathway
#' @param use_normalized Whether to normalize networks before comparison (helps when conditions have different scales)
#' @param normalization_method Method for normalizing data: "max" (divide by max), "sum" (divide by sum), "none" (no normalization)
#' @param order_sectors Character vector specifying the order of sectors
#' @param link_style Character, style of links: "default", "bezier", "straight", or "arc"
#' @param grid_border Logical, whether to draw borders around grid sectors
#' @param grid_border_color Color for grid borders if grid_border=TRUE
#' @param grid_border_width Width of grid borders
#' @param highlight_sources Optional vector of source cell types to highlight
#' @param highlight_targets Optional vector of target cell types to highlight
#' @param highlight_color Color for highlighted cells/links
#' @param text_color Color for labels
#'
#' @return If return_data=TRUE, returns a list containing the plot and data matrices.
#'   Otherwise, returns just the plot.
#' @importFrom circlize circos.clear circos.track chordDiagram
#' @importFrom dplyr filter
#' @importFrom reshape2 melt
#' @importFrom grDevices colorRampPalette pdf dev.off recordPlot
#' @export
netDiff <- function(object.list, comparison = c(1, 2), signaling = NULL,
                    thresh = 0.05, color.use = c("darkorange", "deepskyblue"),
                    cell_color = NULL, title.name = NULL, sources.use = NULL, targets.use = NULL,
                    remove.isolate = FALSE, lab.cex = 0.8, small.gap = 1,
                    big.gap = 10, directional = 0, link.visible = TRUE,
                    link.border = NA, transparency = 0.4, min_weight = 0,
                    return_data = FALSE, save_plot = FALSE,
                    save_name = "netDiff_plot.pdf", save_width = 8, save_height = 8,
                    slot.name = "net", show_all_cell_types = TRUE,
                    use_normalized = FALSE, normalization_method = c("max", "sum", "none"),
                    order_sectors = "default", link_style = "default",
                    grid_border = FALSE, grid_border_color = "grey", grid_border_width = 0.5,
                    highlight_sources = NULL, highlight_targets = NULL,
                    highlight_color = "red", text_color = "black") {

  # Match normalization method argument
  normalization_method <- match.arg(normalization_method)

  # Check link style
  valid_link_styles <- c("default", "bezier", "straight", "arc")
  if (!link_style %in% valid_link_styles) {
    warning(paste("Invalid link_style. Using default. Valid options are:", paste(valid_link_styles, collapse=", ")))
    link_style <- "default"
  }

  # Input validation
  if (length(comparison) != 2) {
    stop("Please provide exactly 2 indices for comparison")
  }

  if (is.null(signaling)) {
    stop("Please provide a signaling pathway name")
  }

  # Get condition names
  condition_names <- names(object.list)
  if (is.null(condition_names)) {
    condition_names <- paste0("Condition_", 1:length(object.list))
  }

  # Get all cell types from both conditions
  all_cellTypes <- c()
  for (i in comparison) {
    all_cellTypes <- union(all_cellTypes, rownames(object.list[[i]]@net$weight))
  }

  # First check if the signaling pathway exists in at least one condition
  signaling_exists <- FALSE
  # Check for existing pathways using LR database
  for (i in 1:2) {
    idx <- comparison[i]
    if (is.list(object.list[[idx]]@LR)) {
      # Search for the pathway in LRsig
      pairLR <- searchPair(signaling = signaling,
                           pairLR.use = object.list[[idx]]@LR$LRsig,
                           key = "pathway_name",
                           matching.exact = TRUE,
                           pair.only = TRUE)
      if (!is.null(pairLR) && nrow(pairLR) > 0) {
        signaling_exists <- TRUE
        break
      }
    }
  }

  if (!signaling_exists) {
    stop(paste0("Signaling pathway '", signaling, "' not found in either condition"))
  }

  # Extract networks from both objects directly from net slot
  extractNetworkFromNet <- function(obj, pathway, idx) {
    # Get network data from net slot
    net <- obj@net

    # Get L-R pairs related to the pathway
    pairLR <- searchPair(signaling = pathway,
                         pairLR.use = obj@LR$LRsig,
                         key = "pathway_name",
                         matching.exact = TRUE,
                         pair.only = TRUE)

    if (is.null(pairLR) || nrow(pairLR) == 0) {
      message(paste("Pathway not found in condition", idx, "- using empty matrix"))
      mat <- matrix(0, length(all_cellTypes), length(all_cellTypes))
      rownames(mat) <- colnames(mat) <- all_cellTypes
      return(mat)
    }

    # Get L-R pair names in the network data
    pairLR.use.name <- dimnames(net$prob)[[3]]

    # Find intersection of pathway-specific L-R pairs and those in the network
    pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)

    if (length(pairLR.name) == 0) {
      message(paste("No matching L-R pairs found for pathway in condition", idx))
      mat <- matrix(0, length(all_cellTypes), length(all_cellTypes))
      rownames(mat) <- colnames(mat) <- all_cellTypes
      return(mat)
    }

    # Filter L-R pairs to those that are in both the pathway and the network
    pairLR <- pairLR[pairLR.name, ]

    # Extract probability and p-value matrices
    prob <- net$prob
    pval <- net$pval

    # Apply p-value threshold (set probabilities to 0 where p-value > threshold)
    prob[pval > thresh] <- 0

    # Check which L-R pairs have non-zero communications after thresholding
    if (length(pairLR.name) > 1) {
      pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 3, sum) != 0]
    } else {
      pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 0]
    }

    if (length(pairLR.name.use) == 0) {
      message(paste("No significant communication for pathway in condition", idx))
      mat <- matrix(0, length(all_cellTypes), length(all_cellTypes))
      rownames(mat) <- colnames(mat) <- all_cellTypes
      return(mat)
    }

    # Extract the relevant probability matrices for the significant L-R pairs
    prob_subset <- prob[, , pairLR.name.use, drop = FALSE]

    # If there's only one L-R pair, ensure it's still a 3D array
    if (length(dim(prob_subset)) == 2) {
      prob_subset <- array(prob_subset, dim = c(dim(prob_subset), 1),
                           dimnames = list(rownames(prob_subset),
                                           colnames(prob_subset),
                                           pairLR.name.use))
    }

    # Aggregate probabilities across all L-R pairs for this pathway
    agg_prob <- apply(prob_subset, c(1, 2), sum)

    message(paste("Successfully extracted network for condition", idx,
                  "with", length(pairLR.name.use), "significant L-R pairs"))

    return(agg_prob)
  }

  # Extract networks
  net1 <- extractNetworkFromNet(object.list[[comparison[1]]], signaling, 1)
  net2 <- extractNetworkFromNet(object.list[[comparison[2]]], signaling, 2)

  # Ensure matrices have the same dimensions and cell types
  all_cells <- union(rownames(net1), rownames(net2))

  # Include all cell types if requested
  if (show_all_cell_types) {
    all_cells <- union(all_cells, all_cellTypes)
  }

  # Create standardized matrices
  std_net1 <- matrix(0, length(all_cells), length(all_cells))
  std_net2 <- matrix(0, length(all_cells), length(all_cells))
  rownames(std_net1) <- colnames(std_net1) <- all_cells
  rownames(std_net2) <- colnames(std_net2) <- all_cells

  # Fill standardized matrices
  for (r in rownames(net1)) {
    for (c in colnames(net1)) {
      if (r %in% all_cells && c %in% all_cells) {
        std_net1[r, c] <- net1[r, c]
      }
    }
  }

  for (r in rownames(net2)) {
    for (c in colnames(net2)) {
      if (r %in% all_cells && c %in% all_cells) {
        std_net2[r, c] <- net2[r, c]
      }
    }
  }

  # Normalize matrices if requested
  if (use_normalized) {
    if (normalization_method == "max") {
      # Normalize by maximum value
      max_val1 <- max(std_net1, na.rm = TRUE)
      max_val2 <- max(std_net2, na.rm = TRUE)
      if (max_val1 > 0) std_net1 <- std_net1 / max_val1
      if (max_val2 > 0) std_net2 <- std_net2 / max_val2
      message("Networks normalized by maximum values")
    } else if (normalization_method == "sum") {
      # Normalize by sum
      sum_val1 <- sum(std_net1, na.rm = TRUE)
      sum_val2 <- sum(std_net2, na.rm = TRUE)
      if (sum_val1 > 0) std_net1 <- std_net1 / sum_val1
      if (sum_val2 > 0) std_net2 <- std_net2 / sum_val2
      message("Networks normalized by sum of values")
    }
  }

  # Calculate network difference
  net_diff <- std_net2 - std_net1

  # Convert to data frame for chord diagram
  net_df <- reshape2::melt(net_diff, value.name = "weight")
  colnames(net_df)[1:2] <- c("source", "target")

  # Filter sources and targets if specified
  if (!is.null(sources.use)) {
    net_df <- net_df %>% dplyr::filter(source %in% sources.use)
  }

  if (!is.null(targets.use)) {
    net_df <- net_df %>% dplyr::filter(target %in% targets.use)
  }

  # Remove cells with zero interaction or below threshold
  net_df <- subset(net_df, abs(weight) > min_weight)

  if (nrow(net_df) == 0) {
    message("No differential interactions found above the threshold.")
    if (return_data) {
      return(list(
        plot = NULL,
        net1 = std_net1,
        net2 = std_net2,
        diff = net_diff
      ))
    } else {
      return(NULL)
    }
  }

  # Get cell types to include - either all from data or just those with interactions
  if (show_all_cell_types && !remove.isolate) {
    cells.use <- all_cells
  } else {
    cells.use <- union(net_df$source, net_df$target)
    # Remove isolated nodes if requested
    if (remove.isolate) {
      cells.use <- union(unique(net_df$source), unique(net_df$target))
    }
  }

  # Determine sector order
  if (is.character(order_sectors) && length(order_sectors) == 1) {
    if (order_sectors == "default") {
      order.sector <- cells.use
    } else if (order_sectors == "strength") {
      # Order by interaction strength
      cell_strength <- sapply(cells.use, function(cell) {
        sum(abs(net_diff[cell, ])) + sum(abs(net_diff[, cell]))
      })
      order.sector <- names(sort(cell_strength, decreasing = TRUE))
    } else if (order_sectors == "alphabetical") {
      order.sector <- sort(cells.use)
    } else if (order_sectors == "custom" && !is.null(sources.use) && !is.null(targets.use)) {
      order.sector <- c(sources.use, setdiff(targets.use, sources.use))
    } else {
      order.sector <- cells.use
    }
  } else if (is.character(order_sectors) && length(order_sectors) > 1) {
    # Use the provided order, but ensure all cells are included
    order.sector <- c(order_sectors, setdiff(cells.use, order_sectors))
  } else {
    order.sector <- cells.use
  }

  # Define grid colors for cell types
  if (!is.null(cell_color)) {
    # Use user-provided colors
    # Ensure all cells have colors assigned
    missing_cells <- setdiff(cells.use, names(cell_color))
    if (length(missing_cells) > 0) {
      # Assign default colors to missing cells
      extra_colors <- colorRampPalette(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                                         "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                                         "#9C755F", "#BAB0AB"))(length(missing_cells))
      cell_color <- c(cell_color, setNames(extra_colors, missing_cells))
    }
    # Only keep colors for cells.use
    grid.col <- cell_color[cells.use]
  } else {
    # Use default palette with better contrast
    if (exists("global_colors") && length(global_colors) > 0) {
      # Use global_colors if available
      default_colors <- colorRampPalette(c(global_colors, "#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                                           "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                                           "#9C755F", "#BAB0AB"))(length(cells.use))
    } else {
      # Fallback color palette
      default_colors <- colorRampPalette(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                                           "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                                           "#9C755F", "#BAB0AB"))(length(cells.use))
    }
    grid.col <- setNames(default_colors, cells.use)
  }

  # Define edge colors for positive and negative changes
  if (length(color.use) != 2) {
    color.use <- c("darkorange", "deepskyblue")  # Orange for increase, Blue for decrease
  }

  # Set default title if not provided
  if (is.null(title.name)) {
    title.name <- paste0(signaling, " Signaling Pathway Differences")
  }

  # Prepare link style parameters
  link_params <- list()
  if (link_style == "straight") {
    link_params$link.shape <- "straight"
  } else if (link_style == "bezier") {
    link_params$link.shape <- "bezier"
  } else if (link_style == "arc") {
    link_params$link.shape <- "arc"
  }

  # Create the chord diagram with enhanced options
  par(mar = c(2, 1, 2, 1), xpd = TRUE)

  # Highlight specific links if requested
  highlight_links <- rep(FALSE, nrow(net_df))
  if (!is.null(highlight_sources) || !is.null(highlight_targets)) {
    if (!is.null(highlight_sources) && !is.null(highlight_targets)) {
      # Highlight links between specified sources and targets
      highlight_links <- (net_df$source %in% highlight_sources) & (net_df$target %in% highlight_targets)
    } else if (!is.null(highlight_sources)) {
      # Highlight all links from specified sources
      highlight_links <- net_df$source %in% highlight_sources
    } else if (!is.null(highlight_targets)) {
      # Highlight all links to specified targets
      highlight_links <- net_df$target %in% highlight_targets
    }
  }

  # Create color vector for links
  link_colors <- ifelse(highlight_links, highlight_color,
                        ifelse(net_df$weight > 0, color.use[1], color.use[2]))

  # Add additional parameters to chordDiagram call
  chord_params <- c(
    list(
      x = net_df,
      order = order.sector,
      col = link_colors,
      grid.col = grid.col,
      transparency = transparency,
      link.border = link.border,
      directional = directional,
      annotationTrack = "grid",
      annotationTrackHeight = c(0.03),
      small.gap = small.gap,
      big.gap = big.gap,
      link.visible = link.visible,
      diffHeight = 0.2,  # Add height difference for clearer direction
      # Add grid border if requested
      grid.border = if (grid_border) grid_border_color else NA
    ),
    link_params
  )

  # Handle grid.border.lwd separately
  if (grid_border) {
    # Set border width globally instead of in chordDiagram parameters
    grid.border.lwd <- grid_border_width
  }

  # Call chordDiagram with all parameters
  do.call(chordDiagram, chord_params)
  library(circlize)
  # Add labels with custom text color
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")

    # Determine if this is a highlighted cell
    is_highlighted <- FALSE
    if (!is.null(highlight_sources) && sector.name %in% highlight_sources) {
      is_highlighted <- TRUE
    }
    if (!is.null(highlight_targets) && sector.name %in% highlight_targets) {
      is_highlighted <- TRUE
    }

    # Set label color based on highlighting
    label_color <- if (is_highlighted) highlight_color else text_color

    circos.text(
      mean(xlim),
      ylim[1],
      sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = lab.cex,
      col = label_color
    )
  }, bg.border = NA)

  # Add title and legend
  if (!is.null(title.name)) {
    text(0, 1.1, title.name, cex = 1.2, font = 2)

    # Add legend for colors
    legend(
      "bottomright",
      legend = c(
        paste("Stronger in", condition_names[comparison[2]]),
        paste("Stronger in", condition_names[comparison[1]])
      ),
      fill = color.use,
      border = NA,
      bty = "n",
      cex = 0.8
    )
  }

  # Store the plot
  plot <- recordPlot()

  # Save the plot if requested
  if (save_plot) {
    if (!is.character(save_name)) {
      save_name <- "netDiff_plot.pdf"
    }

    tryCatch({
      pdf(save_name, width = save_width, height = save_height)
      replayPlot(plot)
      invisible(dev.off())
      message(paste("Plot saved to", save_name))
    }, error = function(e) {
      warning(paste("Failed to save plot:", e$message))
    })
  }

  # Clear circos plot
  circos.clear()

  # Return results based on user preference
  if (return_data) {
    result <- list(
      plot = plot,
      net1 = std_net1,
      net2 = std_net2,
      diff = net_diff,
      signaling = signaling,
      comparison = comparison,
      condition_names = condition_names[comparison],
      cell_types = all_cells,
      params = list(
        thresh = thresh,
        normalization = list(
          used = use_normalized,
          method = normalization_method
        ),
        visual = list(
          colors = color.use,
          order = order_sectors,
          link_style = link_style,
          highlighted_sources = highlight_sources,
          highlighted_targets = highlight_targets
        )
      )
    )

    return(result)
  } else {
    # Return just the plot
    return(plot)
  }
}
