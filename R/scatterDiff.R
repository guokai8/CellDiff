#' Scatter Plot of Pathway-Level Signaling Role Changes
#'
#' @description
#' Creates a scatter plot showing how cell types change their signaling roles
#' between two conditions, using pathway-level information from the netP slot.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (length 2)
#' @param pathways Vector of pathway names to include (NULL for all)
#' @param measure.type "sender", "receiver", or "influence" (default)
#' @param arrow.size Size of trajectory arrows
#' @param arrow.alpha Transparency of trajectory arrows
#' @param thresh P-value threshold for significant interactions (default: 0.05)
#' @param title Plot title
#' @param label.cell Whether to label cell types
#' @param label.size Size of cell type labels (default: 3)
#' @return A ggplot object
#' @export
scatterDiff <- function(object.list, comparison = c(1, 2),
                        pathways = NULL, measure.type = "influence",
                         arrow.size = 1, arrow.alpha = 0.8,
                        thresh = 0.05,
                        title = "Changes in Signaling Roles",
                        label.cell = TRUE,label_size=3) {

  # Global colors for cell types

  if (length(comparison) != 2) {
    stop("Please provide exactly 2 indices for comparison")
  }

  # Extract cell types
  cell.types1 <- rownames(object.list[[comparison[1]]]@netP$prob)
  cell.types2 <- rownames(object.list[[comparison[2]]]@netP$prob)
  cell.types <- intersect(cell.types1, cell.types2)

  if (length(cell.types) == 0) {
    stop("No common cell types found between the two objects")
  }
  # Get pathways if not specified
  pathways1 <- dimnames(object.list[[comparison[1]]]@netP$prob)[[3]]
  pathways2 <- dimnames(object.list[[comparison[2]]]@netP$prob)[[3]]
  if(!is.null(pathways)){
    pathways <- intersect(pathways,intersect(pathways1,pathways2))
    pathways1 <- pathways2<- pathways
  }else{
    pathways <- intersect(pathways1, pathways2)
  }

  if (length(pathways) == 0) {
    stop("No common pathways found between the two objects")
  }

  # Get signaling scores for both conditions
  scores1 <- calculateSignaling(object.list[[comparison[1]]], cell.types1, pathways1, measure.type)
  scores2 <- calculateSignaling(object.list[[comparison[2]]], cell.types2, pathways2, measure.type)

  # Get counts for both conditions and average them
  counts1 <- getInteractionCounts(object.list[[comparison[1]]], signaling = pathways1, slot.name = "net", thresh = thresh)
  counts2 <- getInteractionCounts(object.list[[comparison[2]]], signaling =  pathways2, slot.name = "net", thresh = thresh)
  counts <- abs(counts2 - counts1)

  # Scale counts for dot size (add small value to ensure all dots are visible)
  #scaled_counts <- sqrt(counts) + 1

  # Create data frame for plotting
  plot_data <- data.frame(
    celltype = cell.types,
    condition1 = scores1,
    condition2 = scores2,
    count = counts  )

  # Calculate changes
  plot_data <- plot_data %>%
    dplyr::mutate(
      change = condition2 - condition1,
      perc_change = (condition2 - condition1) / (condition1 + 1e-10) * 100,
      abs_perc_change = abs(perc_change)
    )

  # Create significance categories based on percentage change
  plot_data$change_category <- cut(plot_data$abs_perc_change,
                                   breaks = c(-Inf, 25, 50, 75, 100, Inf),
                                   labels = c("<25%", "25-50%", "50-75%", "75-100%", ">100%"))

  # Set colors for cell types
  n_cells <- length(cell.types)
  idx <- 1:n_cells
  idx[idx > length(global_colors)] <- ((idx[idx > length(global_colors)]-1) %% length(global_colors)) + 1
  color.use <- global_colors[idx]
  names(color.use) <- cell.types

  # Create the plot
  measure_label <- switch(measure.type,
                          "sender" = "Outgoing Signaling Strength",
                          "receiver" = "Incoming Signaling Strength",
                          "influence" = "Signaling Influence")

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = condition1, y = condition2, color = celltype)) +
    # Identity line
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
    # Points with size based on LR count
    ggplot2::geom_point(ggplot2::aes(size = count, alpha = change_category)) +
    # Labels
    ggplot2::labs(
      x = paste(measure_label, "-", names(object.list)[comparison[1]]),
      y = paste(measure_label, "-", names(object.list)[comparison[2]]),
      title = paste0(title, "\n(", length(pathways), " pathways)"),
      alpha = "% Change",
      size = "Count Diff"
    ) +
    ggplot2::scale_color_manual(values = color.use) +
    ggplot2::scale_alpha_discrete(range = c(0.4, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    # Remove color legend for cell types
    ggplot2::guides(color = "none")

  # Add arrows to show direction of change
  p <- p +
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
    )

  if (label.cell) {
    # Add text labels for cell types
    p <- p +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = celltype),
        size = label.size,
        max.overlaps = 100
      )
  }

  return(p)
}
