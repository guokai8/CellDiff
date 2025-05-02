#' Compare Multiple CellChat Objects
#'
#' This function compares multiple CellChat objects to identify global changes in
#' communication patterns across conditions.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (NULL for all objects)
#' @param reference Integer or character, index or name of the reference object
#' @param pathways Vector of pathway names to include (NULL for all)
#' @param measure.type Character, type of measure to calculate: "sender", "receiver", "both", or "influence"
#' @param thresh P-value threshold for significant interactions
#' @param norm.method Method for normalizing values: "none", "z-score", "minmax", or "relative"
#' @param show.pathway Logical, whether to return pathway contribution plot (default: FALSE)
#' @param show.heatmap Logical, whether to return heatmap visualization (default: FALSE)
#' @param return.data Logical, whether to return data along with visualizations (default: FALSE)
#' @param label.timeline Logical, whether to add labels to the timeline plot (default: FALSE)
#' @param label.pathway Logical, whether to add labels to the pathway plot (default: FALSE)
#' @param label.top How many top values to label (default: 5, set to Inf for all)
#' @param label.size Numeric, size of the labels (default: 3)
#'
#' @return By default, returns the timeline plot. If show.pathway or show.heatmap are TRUE,
#'   returns a list of plots. If return.data is TRUE, also includes data in the returned list.
#' @export
compareCellChatsM <- function(object.list, comparison = NULL, reference = NULL,
                              pathways = NULL, measure.type = c("sender", "receiver", "both", "influence"),
                              thresh = 0.05, norm.method = c("none", "z-score", "minmax", "relative"),
                              show.pathway = FALSE, show.heatmap = FALSE, return.data = FALSE,
                              label.timeline = FALSE, label.pathway = FALSE, label.top = 5,
                              label.size = 3) {

  # Match arguments
  measure.type <- match.arg(measure.type)
  norm.method <- match.arg(norm.method)

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

  # Extract all cell types
  all_cell_types <- unique(unlist(lapply(
    object.list[comparison],
    function(obj) rownames(obj@net$weight)
  )))

  # Find common cell types
  common_cell_types <- Reduce(
    intersect,
    lapply(object.list[comparison], function(obj) rownames(obj@net$weight))
  )

  if (length(common_cell_types) == 0) {
    stop("No common cell types found across all conditions")
  }

  # Find all pathways
  if (is.null(pathways)) {
    all_pathways <- unique(unlist(lapply(
      object.list[comparison],
      function(obj) unlist(obj@netP$pathways)
    )))

    # Find common pathways
    pathways <- Reduce(
      intersect,
      lapply(object.list[comparison], function(obj) unlist(obj@netP$pathways))
    )

    if (length(pathways) == 0) {
      warning("No common pathways found. Using pathways from reference object.")
      pathways <- unlist(object.list[[reference]]@netP$pathways)
    }
  }

  # Create data structures to store results
  all_scores <- list()
  pathway_contributions <- list()

  # For each condition, calculate scores
  for (i in 1:length(comparison)) {
    idx <- comparison[i]
    obj <- object.list[[idx]]

    # Calculate selected measure
    if (measure.type %in% c("sender", "receiver", "both")) {
      scores <- calculateSignaling(obj, common_cell_types, pathways, measure.type)
    } else if (measure.type == "influence") {
      # Calculate both sender and receiver, then multiply
      sender_scores <- calculateSignaling(obj, common_cell_types, pathways, "sender")
      receiver_scores <- calculateSignaling(obj, common_cell_types, pathways, "receiver")
      scores <- sender_scores * receiver_scores
    }

    all_scores[[i]] <- scores

    # Only calculate pathway contributions if needed
    if (show.pathway || return.data) {
      pathways_obj <- intersect(pathways, unlist(obj@netP$pathways))
      contributions <- numeric(length(pathways))
      names(contributions) <- pathways

      for (pathway in pathways_obj) {
        # Search for all L-R pairs in this pathway
        lr_pairs <- CellChat::searchPair(
          signaling = pathway,
          pairLR.use = obj@LR$LRsig,
          key = "pathway_name",
          matching.exact = TRUE,
          pair.only = TRUE
        )

        if (nrow(lr_pairs) > 0) {
          # Sum probabilities for all L-R pairs in this pathway
          prob_sum <- 0
          for (lr_pair in rownames(lr_pairs)) {
            if (lr_pair %in% dimnames(obj@net$prob)[[3]]) {
              prob <- obj@net$prob[, , lr_pair]
              pval <- obj@net$pval[, , lr_pair]
              prob[pval > thresh] <- 0
              prob_sum <- prob_sum + sum(prob)
            }
          }
          contributions[pathway] <- prob_sum
        }
      }

      pathway_contributions[[i]] <- contributions
    }
  }

  # Normalize scores if requested
  if (norm.method != "none") {
    normalized_scores <- list()

    if (norm.method == "z-score") {
      # Z-score normalization
      combined_scores <- do.call(cbind, all_scores)
      row_means <- rowMeans(combined_scores)
      row_sds <- apply(combined_scores, 1, sd)

      for (i in 1:length(all_scores)) {
        normalized_scores[[i]] <- (all_scores[[i]] - row_means) / row_sds
      }
    } else if (norm.method == "minmax") {
      # Min-max normalization
      combined_scores <- do.call(cbind, all_scores)
      row_mins <- apply(combined_scores, 1, min)
      row_maxs <- apply(combined_scores, 1, max)
      row_range <- row_maxs - row_mins

      for (i in 1:length(all_scores)) {
        normalized_scores[[i]] <- (all_scores[[i]] - row_mins) / row_range
      }
    } else if (norm.method == "relative") {
      # Relative to reference
      ref_idx <- which(comparison == reference)
      ref_scores <- all_scores[[ref_idx]]

      for (i in 1:length(all_scores)) {
        normalized_scores[[i]] <- all_scores[[i]] / (ref_scores + 1e-10)
      }
    }

    all_scores <- normalized_scores
  }

  # Create data frame for timeline visualization
  score_df <- data.frame()
  for (i in 1:length(comparison)) {
    temp_df <- data.frame(
      celltype = names(all_scores[[i]]),
      score = all_scores[[i]],
      condition = condition_names[comparison[i]]
    )
    score_df <- rbind(score_df, temp_df)
  }

  # Assign global colors for cell types
  celltype_colors <- assignColors(unique(score_df$celltype))

  # Create main visualization (timeline plot)
  p_timeline <- ggplot2::ggplot(score_df, ggplot2::aes(x = condition, y = score, group = celltype, color = celltype)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = celltype_colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("Changes in", measure.type, "signaling across conditions"),
      x = "Condition",
      y = paste(measure.type, "signaling score")
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Add labels to timeline plot if requested
  if (label.timeline) {
    # Find top cell types per condition based on score
    top_labels <- score_df %>%
      dplyr::group_by(condition) %>%
      dplyr::top_n(label.top, score) %>%
      dplyr::ungroup()

    # Add labels to the plot
    p_timeline <- p_timeline +
      ggplot2::geom_text(
        data = top_labels,
        ggplot2::aes(label = celltype),
        hjust = -0.1,
        vjust = 0.5,
        size = label.size,
        check_overlap = TRUE,
        show.legend = FALSE
      )
  }

  # If only timeline plot is needed and no data
  if (!show.pathway && !show.heatmap && !return.data) {
    return(p_timeline)
  }

  # Initialize optional plots
  p_pathway <- NULL
  p_heatmap <- NULL

  # Create pathway contribution visualization if requested
  if (show.pathway || return.data) {
    # Create pathway data frame
    pathway_df <- data.frame()
    for (i in 1:length(comparison)) {
      temp_df <- data.frame(
        pathway = names(pathway_contributions[[i]]),
        contribution = pathway_contributions[[i]],
        condition = condition_names[comparison[i]]
      )
      pathway_df <- rbind(pathway_df, temp_df)
    }

    # Assign global colors for pathways
    pathway_colors <- assignColors(unique(pathway_df$pathway))

    # Create pathway plot
    p_pathway <- ggplot2::ggplot(pathway_df, ggplot2::aes(x = condition, y = contribution, group = pathway, color = pathway)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_color_manual(values = pathway_colors) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = "Pathway contribution changes across conditions",
        x = "Condition",
        y = "Pathway contribution"
      ) +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.minor = ggplot2::element_blank()
      )

    # Add labels to pathway plot if requested
    if (label.pathway) {
      # Find top pathways per condition based on contribution
      top_pathways <- pathway_df %>%
        dplyr::group_by(condition) %>%
        dplyr::top_n(label.top, contribution) %>%
        dplyr::ungroup()

      # Add labels to the plot
      p_pathway <- p_pathway +
        ggplot2::geom_text(
          data = top_pathways,
          ggplot2::aes(label = pathway),
          hjust = -0.1,
          vjust = 0.5,
          size = label.size,
          check_overlap = TRUE,
          show.legend = FALSE
        )
    }
  }

  # Create heatmap visualization if requested and pheatmap is available
  if ((show.heatmap || return.data) && requireNamespace("pheatmap", quietly = TRUE)) {
    # Create matrix for heatmap
    score_matrix <- matrix(0, nrow = length(common_cell_types), ncol = length(comparison))
    rownames(score_matrix) <- common_cell_types
    colnames(score_matrix) <- condition_names[comparison]

    for (i in 1:length(comparison)) {
      score_matrix[, i] <- all_scores[[i]]
    }

    # Use global colors for the annotation
    ann_colors <- list(
      cell_type = celltype_colors[common_cell_types]
    )

    # Create row annotation
    row_ann <- data.frame(cell_type = common_cell_types)
    rownames(row_ann) <- common_cell_types

    # Create heatmap
    p_heatmap <- pheatmap::pheatmap(
      score_matrix,
      scale = "row",
      clustering_method = "ward.D2",
      main = paste("Heatmap of", measure.type, "signaling across conditions"),
      annotation_row = row_ann,
      annotation_colors = ann_colors,
      fontsize = 10,
      silent = TRUE
    )
  }

  # Prepare return list
  result <- list(timeline = p_timeline)

  # Add optional plots
  if (show.pathway && !is.null(p_pathway)) {
    result$pathway <- p_pathway
  }

  if (show.heatmap && !is.null(p_heatmap)) {
    result$heatmap <- p_heatmap
  }

  # Add data if requested
  if (return.data) {
    result$data <- list(
      scores = all_scores,
      condition_names = condition_names[comparison],
      cell_types = common_cell_types,
      celltype_colors = celltype_colors
    )

    # Add pathway data if calculated
    if (!is.null(p_pathway)) {
      result$data$pathway_contributions <- pathway_contributions
      result$data$pathway_colors <- pathway_colors
    }
  }

  return(result)
}


#' Find Common and Unique Signaling Patterns Across Multiple Conditions
#'
#' This function identifies signaling patterns that are common across all conditions,
#' as well as patterns that are unique to specific conditions.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (NULL for all objects)
#' @param pathways Vector of pathway names to include (NULL for all)
#' @param thresh P-value threshold for significant interactions
#' @param min_overlap Minimum number of conditions a pattern must appear in to be considered common
#' @param return_networks Whether to return the full network data
#'
#' @return List of common and unique signaling patterns
#' @export
findCommonUniquePatterns <- function(object.list, comparison = NULL, pathways = NULL,
                                     thresh = 0.05, min_overlap = 2, return_networks = FALSE) {

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

  # Get condition names
  condition_names <- names(object.list)
  if (is.null(condition_names) || any(condition_names == "")) {
    condition_names <- paste0("Condition_", 1:length(object.list))
    names(object.list) <- condition_names
  }

  # Extract all cell types
  all_cell_types <- unique(unlist(lapply(
    object.list[comparison],
    function(obj) rownames(obj@net$weight)
  )))

  # Find common cell types
  common_cell_types <- Reduce(
    intersect,
    lapply(object.list[comparison], function(obj) rownames(obj@net$weight))
  )

  if (length(common_cell_types) == 0) {
    stop("No common cell types found across all conditions")
  }

  # Find all pathways
  if (is.null(pathways)) {
    all_pathways <- unique(unlist(lapply(
      object.list[comparison],
      function(obj) unlist(obj@netP$pathways)
    )))

    # Find common pathways
    pathways <- Reduce(
      intersect,
      lapply(object.list[comparison], function(obj) unlist(obj@netP$pathways))
    )

    if (length(pathways) == 0) {
      pathways <- all_pathways
    }
  }

  # Extract all L-R pairs across all conditions
  all_lr_pairs <- list()

  for (i in 1:length(comparison)) {
    idx <- comparison[i]
    obj <- object.list[[idx]]

    # Get all L-R pairs in this object
    lr_pairs_obj <- list()

    for (pathway in pathways) {
      if (pathway %in% obj@netP$pathways) {
        # Search for L-R pairs in this pathway
        lr_pairs <- CellChat::searchPair(
          signaling = pathway,
          pairLR.use = obj@LR$LRsig,
          key = "pathway_name",
          matching.exact = TRUE,
          pair.only = TRUE
        )

        if (nrow(lr_pairs) > 0) {
          for (lr_pair in rownames(lr_pairs)) {
            if (lr_pair %in% dimnames(obj@net$prob)[[3]]) {
              # Extract sender-target interactions
              prob <- obj@net$prob[, , lr_pair]
              pval <- obj@net$pval[, , lr_pair]
              prob[pval > thresh] <- 0

              # Find significant interactions
              interactions <- which(prob > 0, arr.ind = TRUE)

              if (nrow(interactions) > 0) {
                # Store each interaction
                for (j in 1:nrow(interactions)) {
                  sender <- rownames(prob)[interactions[j, 1]]
                  target <- colnames(prob)[interactions[j, 2]]

                  # Create unique interaction ID
                  interaction_id <- paste(sender, target, pathway, lr_pair, sep = "|")

                  # Store interaction details
                  lr_pairs_obj[[interaction_id]] <- list(
                    sender = sender,
                    target = target,
                    pathway = pathway,
                    lr_pair = lr_pair,
                    weight = prob[interactions[j, 1], interactions[j, 2]]
                  )
                }
              }
            }
          }
        }
      }
    }

    all_lr_pairs[[i]] <- lr_pairs_obj
  }

  # Find common and unique patterns
  interaction_counts <- table(unlist(lapply(all_lr_pairs, names)))

  # Common patterns appear in at least min_overlap conditions
  common_patterns_ids <- names(interaction_counts[interaction_counts >= min_overlap])

  # Unique patterns appear in only one condition
  unique_patterns_ids <- names(interaction_counts[interaction_counts == 1])

  # Create data frames for common and unique patterns
  common_patterns <- data.frame(
    interaction_id = character(),
    sender = character(),
    target = character(),
    pathway = character(),
    lr_pair = character(),
    conditions = character(),
    count = integer(),
    stringsAsFactors = FALSE
  )

  unique_patterns <- data.frame(
    interaction_id = character(),
    sender = character(),
    target = character(),
    pathway = character(),
    lr_pair = character(),
    condition = character(),
    weight = numeric(),
    stringsAsFactors = FALSE
  )

  # Process common patterns
  for (id in common_patterns_ids) {
    # Find which conditions have this pattern
    conditions_with_pattern <- c()

    for (i in 1:length(comparison)) {
      if (id %in% names(all_lr_pairs[[i]])) {
        conditions_with_pattern <- c(conditions_with_pattern, condition_names[comparison[i]])
      }
    }

    # Get pattern details from first occurrence
    for (i in 1:length(comparison)) {
      if (id %in% names(all_lr_pairs[[i]])) {
        pattern <- all_lr_pairs[[i]][[id]]

        common_patterns <- rbind(
          common_patterns,
          data.frame(
            interaction_id = id,
            sender = pattern$sender,
            target = pattern$target,
            pathway = pattern$pathway,
            lr_pair = pattern$lr_pair,
            conditions = paste(conditions_with_pattern, collapse = ","),
            count = length(conditions_with_pattern),
            stringsAsFactors = FALSE
          )
        )

        break  # We only need one instance
      }
    }
  }

  # Process unique patterns
  for (id in unique_patterns_ids) {
    # Find which condition has this pattern
    for (i in 1:length(comparison)) {
      if (id %in% names(all_lr_pairs[[i]])) {
        pattern <- all_lr_pairs[[i]][[id]]

        unique_patterns <- rbind(
          unique_patterns,
          data.frame(
            interaction_id = id,
            sender = pattern$sender,
            target = pattern$target,
            pathway = pattern$pathway,
            lr_pair = pattern$lr_pair,
            condition = condition_names[comparison[i]],
            weight = pattern$weight,
            stringsAsFactors = FALSE
          )
        )

        break  # This pattern should only be in one condition
      }
    }
  }

  # Return results
  if (return_networks) {
    return(list(
      common_patterns = common_patterns,
      unique_patterns = unique_patterns,
      all_interactions = all_lr_pairs
    ))
  } else {
    return(list(
      common_patterns = common_patterns,
      unique_patterns = unique_patterns
    ))
  }
}

#' Extract Network From CellChat Object for Multiple Groups Analysis
#'
#' Helper function to extract network data from a CellChat object for a specific signaling pathway.
#'
#' @param obj A CellChat object
#' @param pathway Signaling pathway name
#' @param idx Index identifier (for error messages)
#' @param all_cellTypes Vector of all cell types to include
#' @param sources.use Optional vector of source cell types to filter
#' @param targets.use Optional vector of target cell types to filter
#' @param thresh P-value threshold for significant interactions
#'
#' @return A matrix of network connections
#' @keywords internal
extractNetworkFromNet <- function(obj, pathway, idx, all_cellTypes = NULL,
                                  sources.use = NULL, targets.use = NULL, thresh = 0.05) {

  # Get network data from net slot
  net <- obj@net

  # Get L-R pairs related to the pathway
  pairLR <- CellChat::searchPair(signaling = pathway,
                                 pairLR.use = obj@LR$LRsig,
                                 key = "pathway_name",
                                 matching.exact = TRUE,
                                 pair.only = TRUE)

  if (is.null(pairLR) || nrow(pairLR) == 0) {
    message(paste("Pathway not found in condition", idx, "- using empty matrix"))
    if (is.null(all_cellTypes)) {
      all_cellTypes <- rownames(net$weight)
    }
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
    if (is.null(all_cellTypes)) {
      all_cellTypes <- rownames(net$weight)
    }
    mat <- matrix(0, length(all_cellTypes), length(all_cellTypes))
    rownames(mat) <- colnames(mat) <- all_cellTypes
    return(mat)
  }

  # Filter L-R pairs to those that are in both the pathway and the network
  pairLR <- pairLR[pairLR.name, , drop = FALSE]

  # Extract probability and p-value matrices
  prob <- net$prob
  pval <- net$pval

  # Apply p-value threshold (set probabilities to 0 where p-value > threshold)
  prob[pval > thresh] <- 0

  # Process sources.use parameter if specified
  if (!is.null(sources.use)) {
    if (is.character(sources.use)) {
      if (all(sources.use %in% rownames(prob))) {
        sources.idx <- match(sources.use, rownames(prob))
      } else {
        sources.idx <- intersect(1:nrow(prob), sources.use)
      }
    } else {
      sources.idx <- sources.use
    }
    excluded_sources <- setdiff(1:nrow(prob), sources.idx)
    if (length(excluded_sources) > 0) {
      prob[excluded_sources, , ] <- 0
    }
  }

  # Process targets.use parameter if specified
  if (!is.null(targets.use)) {
    if (is.character(targets.use)) {
      if (all(targets.use %in% colnames(prob))) {
        targets.idx <- match(targets.use, colnames(prob))
      } else {
        targets.idx <- intersect(1:ncol(prob), targets.use)
      }
    } else {
      targets.idx <- targets.use
    }
    excluded_targets <- setdiff(1:ncol(prob), targets.idx)
    if (length(excluded_targets) > 0) {
      prob[, excluded_targets, ] <- 0
    }
  }

  # Check which L-R pairs have non-zero communications after thresholding
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name, drop = FALSE], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name, drop = FALSE]) != 0]
  }

  if (length(pairLR.name.use) == 0) {
    message(paste("No significant communication for pathway in condition", idx))
    if (is.null(all_cellTypes)) {
      all_cellTypes <- rownames(net$weight)
    }
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

  return(agg_prob)
}

#' Create Standardized Matrix With Consistent Dimensions
#'
#' Helper function to create a standardized matrix with all cell types.
#'
#' @param mat Original matrix
#' @param all_cells Vector of all cell types to include
#'
#' @return A standardized matrix with all cells
#' @keywords internal
createStandardizedMatrix <- function(mat, all_cells) {
  # Create standardized matrix
  std_mat <- matrix(0, length(all_cells), length(all_cells))
  rownames(std_mat) <- colnames(std_mat) <- all_cells

  # Fill with values from original matrix
  common_rows <- intersect(rownames(mat), all_cells)
  common_cols <- intersect(colnames(mat), all_cells)

  if (length(common_rows) > 0 && length(common_cols) > 0) {
    for (r in common_rows) {
      for (c in common_cols) {
        std_mat[r, c] <- mat[r, c]
      }
    }
  }

  return(std_mat)
}

#' Normalize Matrix by Maximum or Sum
#'
#' Helper function to normalize a matrix by its maximum value or sum.
#'
#' @param mat Matrix to normalize
#' @param method Normalization method: "max", "sum", or "none"
#'
#' @return Normalized matrix
#' @keywords internal
normalizeMatrix <- function(mat, method = c("max", "sum", "none")) {
  method <- match.arg(method)

  if (method == "max") {
    max_val <- max(mat, na.rm = TRUE)
    if (max_val > 0) {
      return(mat / max_val)
    }
  } else if (method == "sum") {
    sum_val <- sum(mat, na.rm = TRUE)
    if (sum_val > 0) {
      return(mat / sum_val)
    }
  }

  # Default: return original matrix
  return(mat)
}

#' Get Common Cell Types Across CellChat Objects
#'
#' This function extracts the cell types that are common across multiple CellChat objects.
#'
#' @param object.list List of CellChat objects
#' @param indices Vector of indices to consider
#'
#' @return Vector of common cell type names
#' @export
getCommonCellTypes <- function(object.list, indices = NULL) {
  if (is.null(indices)) {
    indices <- seq_along(object.list)
  }

  # Extract cell types from each object
  cell_types <- lapply(object.list[indices], function(x) rownames(x@net$weight))

  # Find intersection
  common_cells <- Reduce(intersect, cell_types)

  return(common_cells)
}

#' Get Common Pathways Across CellChat Objects
#'
#' This function extracts the signaling pathways that are common across multiple CellChat objects.
#'
#' @param object.list List of CellChat objects
#' @param indices Vector of indices to consider
#'
#' @return Vector of common pathway names
#' @export
getCommonPathways <- function(object.list, indices = NULL) {
  if (is.null(indices)) {
    indices <- seq_along(object.list)
  }

  # Extract pathways from each object
  pathways <- lapply(object.list[indices], function(x) {
    if ("pathways" %in% names(x@netP)) {
      return(names(x@netP$pathways))
    } else {
      return(character(0))
    }
  })

  # Find intersection
  common_pathways <- Reduce(intersect, pathways)

  return(common_pathways)
}

#' Global color palette for CellDiff visualizations
#'
#' A consistent color palette used across all CellDiff visualization functions
#' for cell types and other categorical variables.
#'
#' @export
global_colors <- c(
  "#2670B8", "#B3D9CE", "#16A16E","#36BD9C", "#206A9A", "#5499B1",
  "#876AAA", "#A5CDA7", "#9E1B8E", "#B85F81",
  "#405BA0", "#C0F5D7", "#223C4D", "#4758A0", "#3F578B", "#C2DF7E", "#71B9A3",
  "#BDFBDA", "#D7E7D1", "#F8C4C9", "#E79292", "#6C9D8E", "#58B390", "#B6DDCB", "#E0F0E5",
  "#16416E", "#5499B2", "#83BBB7", "#B5D6B9", "#DAE6D3", "#6FB1B9", "#A8D7D5", "#B7CDD1",
  "#82BCD8", "#AFE8F3", "#B5F8F1", "#494574", "#427794", "#4DA8A0", "#94CFB2", "#6D8C87",
  "#87B6AA", "#C9E1BE", "#4D85A0", "#82BED8", "#C2DFE7", "#91BE93", "#BFD8D4", "#D35400",
  "#E67E22", "#F1C40F", "#FBC02D", "#EC407A", "#F06292", "#00ACC1", "#26C6DA", "#C0392B"
)

#' Assign global colors to items
#'
#' This function assigns colors from the global_colors palette to a vector of items.
#' Colors will be recycled if there are more items than colors.
#'
#' @param items Vector of item names to assign colors to
#' @param shuffle Whether to shuffle the colors, default is FALSE
#'
#' @return A named vector of colors, with names matching the input items
#' @export
assignColors <- function(items, shuffle = FALSE) {
  n_items <- length(items)
  n_colors <- length(global_colors)

  # Use all colors or recycle them if needed
  idx <- 1:n_items
  idx[idx > n_colors] <- ((idx[idx > n_colors] - 1) %% n_colors) + 1

  # Shuffle colors if requested
  if (shuffle) {
    set.seed(42)  # For reproducibility
    idx <- sample(idx)
  }

  # Assign colors to items
  assigned_colors <- global_colors[idx]
  names(assigned_colors) <- items

  return(assigned_colors)
}
