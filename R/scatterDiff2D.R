#' 2D Scatter Plot of Sender and Receiver Pathway-Level Changes
#'
#' @description
#' Creates a 2D scatter plot showing how cell types change in both sender and receiver
#' dimensions between two conditions, using pathway-level information from the netP slot.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (length 2)
#' @param pathways Vector of pathway names to include (NULL for all)
#' @param cell.type.strategy Character, strategy for aligning cell types: "shared" (uses intersection) or "union" (default, uses all cell types)
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
#' @param point.shapes Vector of shapes for different conditions (default: c(16, 17))
#' @param point.alpha.values Vector of alpha values for conditions (default: c(0.7, 1))
#' @param margin.factor.x Horizontal expansion factor for plot margins (default: 0.3)
#' @param margin.factor.y Vertical expansion factor for plot margins (default: 0.3)
#' @param colors Vector of colors for cell types (default: NULL, uses global_colors)
#' @param title Plot title (default: "Changes in Sender-Receiver Roles")
#' @param title.size Size of plot title (default: 14)
#' @param title.face Font face of plot title (default: "bold")
#' @param axis.text.size Size of axis text (default: 10)
#' @param axis.title.size Size of axis titles (default: 12)
#' @param legend.position Position of the legend (default: "right")
#' @param base.theme Base theme for the plot (default: ggplot2::theme_bw())
#' @param return.data Whether to return data instead of plot (default: FALSE)
#'
#' @return A ggplot object or a list of data frames if return.data is TRUE
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(celllist)
#'
#' # Basic 2D scatter plot (union strategy by default for this function)
#' plot <- scatterDiff2D(celllist, comparison = c(1, 2))
#'
#' # Use shared strategy instead
#' plot <- scatterDiff2D(celllist, comparison = c(1, 2),
#'                       cell.type.strategy = "shared")
#'
#' # Customize visualization
#' plot <- scatterDiff2D(celllist, comparison = c(1, 2),
#'                       cell.type.strategy = "union",
#'                       log.scale = TRUE,
#'                       label.cell = TRUE,
#'                       add.quadrants = TRUE,
#'                       arrow.size = 1.2,
#'                       title = "Sender-Receiver Role Changes")
#'
#' # Analyze specific pathways
#' plot <- scatterDiff2D(celllist, comparison = c(1, 2),
#'                       pathways = c("CXCL", "CCL"),
#'                       cell.type.strategy = "union")
#' }
#'
#' @export
scatterDiff2D <- function(object.list, comparison = c(1, 2),
                          pathways = NULL,
                          cell.type.strategy = c("union", "shared"),
                          thresh = 0.05,
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
                          point.shapes = c(16, 17),
                          point.alpha.values = c(0.7, 1),
                          margin.factor.x = 0.3,
                          margin.factor.y = 0.3,
                          colors = NULL,
                          title = "Changes in Sender-Receiver Roles",
                          title.size = 14, title.face = "bold",
                          axis.text.size = 10, axis.title.size = 12,
                          legend.position = "right",
                          base.theme = ggplot2::theme_bw(),
                          return.data = FALSE) {

  # Validate inputs
  cell.type.strategy <- match.arg(cell.type.strategy)

  if (length(comparison) != 2) {
    stop("Please provide exactly 2 indices for comparison")
  }

  # Align cell types using the specified strategy
  alignment <- alignCellTypes(object.list, indices = comparison, strategy = cell.type.strategy)
  cell.types <- alignment$cell_types

  # Get cell types available in each object
  cell.types1 <- rownames(object.list[[comparison[1]]]@netP$prob)
  cell.types2 <- rownames(object.list[[comparison[2]]]@netP$prob)

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


  # Get scores for both conditions (use aligned cell types)
  scores1 <- calculateScores(object.list[[comparison[1]]], cell.types, pathways1, thresh)
  scores2 <- calculateScores(object.list[[comparison[2]]], cell.types, pathways2, thresh)

  # Get counts for both conditions and average them
  counts1_all <- getInteractionCounts(object.list[[comparison[1]]], signaling = pathways1, slot.name = "net", thresh = thresh)
  counts2_all <- getInteractionCounts(object.list[[comparison[2]]], signaling =  pathways2, slot.name = "net", thresh = thresh)

  # Align counts to cell.types (some may be missing in union strategy)
  counts1 <- sapply(cell.types, function(ct) {
    if (ct %in% names(counts1_all)) counts1_all[ct] else 0
  })
  counts2 <- sapply(cell.types, function(ct) {
    if (ct %in% names(counts2_all)) counts2_all[ct] else 0
  })
  counts <- abs(counts2 - counts1)

  # Scale counts for dot size (add small value to ensure all dots are visible)
  scaled_counts <- sqrt(counts) + 1

  # Create data frames for plotting
  data1 <- data.frame(
    celltype = cell.types,
    sender = scores1$sender,
    receiver = scores1$receiver,
    condition = rep(names(object.list)[comparison[1]], length(cell.types)),
    count = counts1,
    lr_count = counts
  )

  data2 <- data.frame(
    celltype = cell.types,
    sender = scores2$sender,
    receiver = scores2$receiver,
    condition = rep(names(object.list)[comparison[2]], length(cell.types)),
    count = counts2,
    lr_count = counts
  )

  # Apply log transformation if requested
  if (log.scale) {
    data1$sender <- log1p(data1$sender)
    data1$receiver <- log1p(data1$receiver)
    data2$sender <- log1p(data2$sender)
    data2$receiver <- log1p(data2$receiver)
  }

  # Combine data frames
  plot_data <- rbind(data1, data2)
  x_range <- range(plot_data$sender, na.rm = TRUE)
  y_range <- range(plot_data$receiver, na.rm = TRUE)

  # Expand the range by user-defined factors
  x_span <- diff(x_range)
  y_span <- diff(y_range)

  x_min <- x_range[1] - margin.factor.x * x_span
  x_max <- x_range[2] + margin.factor.x * x_span
  y_min <- y_range[1] - margin.factor.y * y_span
  y_max <- y_range[2] + margin.factor.y * y_span

  # Calculate changes for each cell type
  changes <- data.frame(
    celltype = cell.types,
    sender_change = data2$sender - data1$sender,
    receiver_change = data2$receiver - data1$receiver,
    total_change = sqrt((data2$sender - data1$sender)^2 + (data2$receiver - data1$receiver)^2),
    sender_direction = ifelse(data2$sender > data1$sender, "increase", "decrease"),
    receiver_direction = ifelse(data2$receiver > data1$receiver, "increase", "decrease"),
    lr_count = counts
  )

  # Assign quadrants
  changes$quadrant <- ifelse(changes$sender_change > 0 & changes$receiver_change > 0, "I",
                             ifelse(changes$sender_change < 0 & changes$receiver_change > 0, "II",
                                    ifelse(changes$sender_change < 0 & changes$receiver_change < 0, "III", "IV")))

  # Set colors for cell types
  n_cells <- length(cell.types)

  # Use provided colors or default to global_colors
  if (is.null(colors)) {
    idx <- 1:n_cells
    idx[idx > length(global_colors)] <- ((idx[idx > length(global_colors)]-1) %% length(global_colors)) + 1
    color.use <- global_colors[idx]
  } else {
    # If custom colors provided but not enough, recycle them
    if (length(colors) < n_cells) {
      colors <- rep(colors, length.out = n_cells)
    }
    color.use <- colors[1:n_cells]
  }
  names(color.use) <- cell.types

  # Create arrows data
  arrows_data <- data.frame(
    celltype = cell.types,
    x = data1$sender,
    y = data1$receiver,
    xend = data2$sender,
    yend = data2$receiver,  # Fixed: correctly use receiver here
    lr_count = changes$lr_count
  )

  # Return data if requested
  if (return.data) {
    return(list(
      data1 = data1,
      data2 = data2,
      changes = changes,
      arrows = arrows_data
    ))
  }

  # Create the plot
  p <- ggplot2::ggplot() +
    # Add points for both conditions
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = sender, y = receiver, color = celltype, size = count,
                   shape = condition, alpha = condition)
    )

  # Add arrows showing direction of change - handle color differently depending on arrow.color
  if (is.null(arrow.color)) {
    # Use celltype-based colors from aes mapping
    p <- p + ggplot2::geom_segment(
      data = arrows_data,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = celltype),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(arrow.size * 0.1, "inches"),
        type = arrow.type
      ),
      alpha = arrow.alpha
    )
  } else {
    # Use specified arrow color
    p <- p + ggplot2::geom_segment(
      data = arrows_data,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(arrow.size * 0.1, "inches"),
        type = arrow.type
      ),
      alpha = arrow.alpha,
      color = arrow.color
    )
  }

  # Continue with the rest of the plot
  p <- p +
    # Set colors, shapes, etc.
    ggplot2::scale_color_manual(values = color.use) +
    ggplot2::scale_shape_manual(values = point.shapes) +
    ggplot2::scale_alpha_manual(values = point.alpha.values) +
    ggplot2::scale_size_continuous(range = point.size.range) +
    # Labels
    ggplot2::labs(
      x = ifelse(log.scale, "Sender Role (log scale)", "Sender Role"),
      y = ifelse(log.scale, "Receiver Role (log scale)", "Receiver Role"),
      title = paste0(title, "\n(", length(pathways), " pathways)"),
      shape = "Condition",
      alpha = "Condition",
      size = "LR Count"
    ) +
    base.theme +
    ggplot2::theme(
      legend.position = legend.position,
      plot.title = ggplot2::element_text(hjust = 0.5, size = title.size, face = title.face),
      axis.text = ggplot2::element_text(size = axis.text.size),
      axis.title = ggplot2::element_text(size = axis.title.size)
    ) +
    # Remove color legend for cell types
    ggplot2::guides(color = "none")

  # Add quadrant lines if requested
  if (add.quadrants) {
    # Find midpoints for reference lines
    x_mid <- mean(c(min(plot_data$sender), max(plot_data$sender)))
    y_mid <- mean(c(min(plot_data$receiver), max(plot_data$receiver)))

    p <- p +
      # Add reference lines
      ggplot2::geom_hline(yintercept = y_mid,
                          linetype = quadrant.line.type,
                          color = quadrant.line.color,
                          alpha = quadrant.line.alpha) +
      ggplot2::geom_vline(xintercept = x_mid,
                          linetype = quadrant.line.type,
                          color = quadrant.line.color,
                          alpha = quadrant.line.alpha) +
      # Add quadrant labels
      ggplot2::annotate("text",
                        x = max(plot_data$sender) * 1.1,
                        y = max(plot_data$receiver) * 1.2,
                        label = "High sender\nHigh receiver",
                        size = quadrant.label.size,
                        alpha = quadrant.label.alpha) +
      ggplot2::annotate("text",
                        x = x_min * 0.5,
                        y = max(plot_data$receiver) * 1.2,
                        label = "Low sender\nHigh receiver",
                        size = quadrant.label.size,
                        alpha = quadrant.label.alpha) +
      ggplot2::annotate("text",
                        x = x_min * 0.5,
                        y = y_min * 0.5,
                        label = "Low sender\nLow receiver",
                        size = quadrant.label.size,
                        alpha = quadrant.label.alpha) +
      ggplot2::annotate("text",
                        x = max(plot_data$sender) * 1.1,
                        y = y_min * 0.5,
                        label = "High sender\nLow receiver",
                        size = quadrant.label.size,
                        alpha = quadrant.label.alpha)
  }

  # Add text labels for cell types if requested
  if (label.cell) {
    # Choose which condition's points to label (usually the second one)
    label_data <- data2

    if (label.repel) {
      # Use ggrepel for smart label placement
      p <- p +
        ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(x = sender, y = receiver, label = celltype),
          size = label.size,
          color = label.color,
          box.padding = ggplot2::unit(0.35, "lines"),
          point.padding = ggplot2::unit(0.3, "lines"),
          min.segment.length = label.min.segment.length,
          max.overlaps = label.max.overlaps,
          segment.color = "grey50",
          segment.alpha = 0.6,
          seed = 42  # For reproducibility
        )
    } else {
      # Use regular text labels
      p <- p +
        ggplot2::geom_text(
          data = label_data,
          ggplot2::aes(x = sender, y = receiver, label = celltype),
          hjust = -0.2,
          vjust = 1.5,
          size = label.size,
          color = label.color,
          max.overlaps = 100
        )
    }
  }

  # Set limits
  p <- p + ggplot2::xlim(c(x_min, x_max)) + ggplot2::ylim(c(y_min, y_max))

  return(p)
}
