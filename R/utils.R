#' Extract Network From CellChat Object
#'
#' This function extracts network data from a CellChat object for specific signaling pathways.
#'
#' @param object A CellChat object
#' @param signaling Signaling pathway name
#' @param slot.name Slot name containing the network information, default is "netP"
#' @param thresh Significance threshold, default is 0.05
#'
#' @return A matrix of network connections
#' @importFrom CellChat searchPair
#' @export
extractNet <- function(object, signaling = NULL, slot.name = "netP", thresh = 0.05) {
  pairLR <- CellChat::searchPair(signaling = signaling,
                                 pairLR.use = object@LR$LRsig,
                                 key = "pathway_name",
                                 matching.exact = TRUE,
                                 pair.only = FALSE)

  net <- slot(object, slot.name)
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]

  prob <- net$prob
  pval <- net$pval

  # Filter by significance threshold
  prob[pval > thresh] <- 0
  pairLR.name.use <- pairLR.name

  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }

  # Extract probability matrix for selected pathways
  prob <- prob[,,pairLR.name.use]

  # Handle single pathway case
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
  }

  # Sum probabilities across pathways
  net <- apply(prob, c(1,2), sum)
  return(net)
}

#' Calculate Significance Between Groups
#'
#' This function calculates significance of differences between two groups.
#'
#' @param x Matrix of values
#' @param group Vector of group names
#'
#' @return p-value from comparison
#' @importFrom stats wilcox.test
#' @keywords internal
calculateSignificance <- function(x, group) {
  d <- x[, group]
  d <- as.data.frame(d)
  d <- d[rowSums(d, na.rm = TRUE) != 0, , drop = FALSE]

  if(nrow(d) >= 3 & sum(is.na(d)) == 0) {
    p <- wilcox.test(d[,1], d[,2], paired = TRUE)$p.value
  } else if(nrow(d) < 3) {
    p <- 1  # Not enough data for test
  } else {
    p <- NA  # Missing values
  }
  return(p)
}

#' Calculate Relative Changes Between Groups
#'
#' This function calculates relative changes between two groups.
#'
#' @param df Data frame of contributions
#' @param group Vector of two group names
#'
#' @return Data frame with relative contribution values
#' @importFrom tidyr spread
#' @importFrom dplyr filter select
#' @keywords internal
calculateRelativeChanges <- function(df, group) {
  dd <- df %>%
    dplyr::filter(Group %in% group) %>%
    dplyr::select(Group, name, contribution) %>%
    tidyr::spread(Group, contribution)

  dd <- as.data.frame(dd)
  dd$contribution.relative <- dd[, group[2]] / dd[, group[1]]
  return(dd)
}

#' Filter Groups From Data Frame
#'
#' This function filters a data frame by group and removes empty rows.
#'
#' @param x Data frame
#' @param select Vector of group names to select
#'
#' @return Filtered data frame
#' @importFrom dplyr filter
#' @keywords internal
filterGroups <- function(x, select, pair.name.all) {
  dx <- x %>% dplyr::filter(Group %in% select)

  for (i in 1:length(pair.name.all)) {
    df.t <- dx[dx$name == pair.name.all[i], "contribution"]
    if (sum(df.t) == 0) {
      dx <- dx[-which(dx$name == pair.name.all[i]), ]
    }
  }
  return(dx)
}

#' Prepare CellChat List for Comparison
#'
#' This function prepares a list of CellChat objects for comparison.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare
#' @param slot.name Slot name containing network information
#' @param measure Measure to use for comparison ("weight" or "count")
#' @param thresh Significance threshold
#'
#' @return List of prepared data for comparison
#' @importFrom methods slot
#' @keywords internal
prepareCellChatComparison <- function(object.list, comparison, slot.name = "netP",
                                      measure = "weight", thresh = 0.05) {
  prob.list <- list()
  pSum <- list()
  pSum.original <- list()
  pair.name <- list()
  idx <- list()
  pSum.original.all <- c()
  object.names.comparison <- names(object.list)[comparison]

  for (i in 1:length(comparison)) {
    object.name <- names(object.list)[comparison[i]]
    object.data <- slot(object.list[[comparison[i]]], slot.name)

    prob <- object.data$prob
    prob[object.data$pval > thresh] <- 0

    if (measure == "count") {
      prob <- 1 * (prob > 0)
    }

    prob.list[[i]] <- prob

    if (sum(prob) == 0) {
      stop("No inferred communications for the input!")
    }

    pSum.original[[i]] <- apply(prob, 3, sum)

    if (measure == "weight") {
      pSum[[i]] <- -1/log(pSum.original[[i]])
      pSum[[i]][is.na(pSum[[i]])] <- 0
      idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)
      pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
    } else if (measure == "count") {
      pSum[[i]] <- pSum.original[[i]]
    }

    pair.name[[i]] <- names(pSum.original[[i]])
  }

  # Replace values with max pSum for infinity cases
  if (measure == "weight") {
    values.assign <- seq(max(unlist(pSum))*1.1, max(unlist(pSum))*1.5, length.out = length(unlist(idx)))
    position <- sort(pSum.original.all, index.return = TRUE)$ix

    for (i in 1:length(comparison)) {
      if (i == 1) {
        pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), position)]
      } else {
        pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i-1]))+1:length(unlist(idx[1:i])), position)]
      }
    }
  }

  # Extract all pathway names
  pair.name.all <- as.character(unique(unlist(pair.name)))

  return(list(
    prob.list = prob.list,
    pSum = pSum,
    pSum.original = pSum.original,
    pair.name = pair.name,
    pair.name.all = pair.name.all,
    object.names.comparison = object.names.comparison
  ))
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

#' Create a continuous color palette
#'
#' This function creates a continuous color palette for use in heatmaps
#' and other visualizations with continuous data.
#'
#' @param low Color for low values, default is "darkblue"
#' @param mid Color for middle values, default is "white"
#' @param high Color for high values, default is "red"
#' @param n Number of colors in the palette, default is 100
#'
#' @return A vector of colors forming a continuous palette
#' @importFrom grDevices colorRampPalette
#' @export
continuousColors <- function(low = "darkblue", mid = "white", high = "red", n = 100) {
  colorRampPalette(c(low, mid, high))(n)
}

#' Check for CellChat objects and validate them
#'
#' This function checks if the provided objects are valid CellChat objects
#' with all required slots and data structures.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare
#'
#' @return TRUE if all objects are valid, throws an error otherwise
#' @keywords internal
validateCellChatObjects <- function(object.list, comparison = NULL) {
  # Check if object.list is a list
  if (!is.list(object.list)) {
    stop("object.list must be a list of CellChat objects", call. = FALSE)
  }

  # If comparison is provided, validate those specific indices
  if (!is.null(comparison)) {
    if (!all(comparison %in% seq_along(object.list))) {
      stop(sprintf("comparison indices must be between 1 and %d", length(object.list)),
           call. = FALSE)
    }
    indices <- comparison
  } else {
    indices <- seq_along(object.list)
  }

  # Check each object
  for (i in indices) {
    obj <- object.list[[i]]

    # Check if it's a CellChat object
    if (!inherits(obj, "CellChat")) {
      stop(sprintf("Object at index %d is not a CellChat object", i), call. = FALSE)
    }

    # Check for required slots
    required_slots <- c("net", "netP")
    for (slot in required_slots) {
      if (!slot %in% names(obj)) {
        stop(sprintf("Object at index %d is missing required slot '%s'", i, slot),
             call. = FALSE)
      }
    }

    # Check for weight matrix
    if (!"weight" %in% names(obj$net)) {
      stop(sprintf("Object at index %d is missing the weight matrix in the net slot", i),
           call. = FALSE)
    }
  }

  # Check that objects share common cell types if comparing
  if (length(indices) > 1) {
    cell_types <- lapply(object.list[indices], function(x) rownames(x@net$weight))
    common_cells <- Reduce(intersect, cell_types)

    if (length(common_cells) == 0) {
      warning("No common cell types found between the objects", call. = FALSE)
    }
  }

  return(TRUE)
}

#' Get common cell types across CellChat objects
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

  validateCellChatObjects(object.list, indices)

  cell_types <- lapply(object.list[indices], function(x) rownames(x@net$weight))
  common_cells <- Reduce(intersect, cell_types)

  return(common_cells)
}

#' Get common pathways across CellChat objects
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

  validateCellChatObjects(object.list, indices)

  pathways <- lapply(object.list[indices], function(x) {
    if ("pathways" %in% names(x@netP)) {
      return(names(x@netP$pathways))
    } else {
      return(character(0))
    }
  })

  common_pathways <- Reduce(intersect, pathways)

  return(common_pathways)
}

#' Simplified Aggregate Cell-Cell Communication Network
#'
#' This function creates an aggregated cell-cell communication network by
#' combining information across multiple ligand-receptor pairs or signaling pathways.
#'
#' @param object A CellChat object
#' @param signaling Vector of pathway names to include (NULL for all)
#' @param thresh P-value threshold for significant interactions (default: 0.05)
#' @param slot.name Slot name containing the network information (default: "net")
#' @param remove.isolate Whether to remove isolated cell types with no interactions
#' @param return.object Whether to return updated CellChat object or just the network
#'
#' @return Updated CellChat object or network structure, depending on return.object
#' @importFrom methods slot
#' @importFrom dplyr group_by summarize
#' @export
aggregateNetSimple <- function(object, signaling = NULL, thresh = 0.05,
                               slot.name = "net", remove.isolate = TRUE,
                               return.object = TRUE) {

  # Get the network from object
  net <- slot(object, slot.name)

  if (is.null(signaling)) {
    # Process all pathways
    prob <- net$prob
    pval <- net$pval

    # Apply threshold
    pval[prob == 0] <- 1  # Set p-value to 1 where probability is 0
    prob[pval >= thresh] <- 0  # Set probability to 0 where p-value >= threshold

    # Calculate count and weight matrices
    net$count <- apply(prob > 0, c(1, 2), sum)
    net$weight <- apply(prob, c(1, 2), sum)

  } else {
    # Filter communication based on signaling pathways
    df.net <- subsetCommunication(object, slot.name = slot.name,
                                  signaling = signaling, thresh = thresh)

    # Create source_target combined key
    df.net$source_target <- paste(df.net$source, df.net$target, sep = "_")

    # Group by source_target and calculate count and sum of probabilities
    df.counts <- df.net %>%
      dplyr::group_by(source_target) %>%
      dplyr::summarize(count = n(), prob = sum(prob), .groups = "drop")

    # Split source_target back into source and target
    source_target_split <- stringr::str_split(df.counts$source_target, "_", simplify = TRUE)
    df.counts$source <- source_target_split[, 1]
    df.counts$target <- source_target_split[, 2]

    # Get cell levels from object
    cells.level <- levels(object@idents)

    # Handle isolate cells
    if (remove.isolate) {
      # Keep only cell types involved in interactions
      cells.keep <- unique(c(df.counts$source, df.counts$target))
      cells.level <- cells.level[cells.level %in% cells.keep]
    }

    # Convert source and target to factors with proper levels
    df.counts$source <- factor(df.counts$source, levels = cells.level)
    df.counts$target <- factor(df.counts$target, levels = cells.level)

    # Create count and weight matrices
    count <- tapply(df.counts$count, list(df.counts$source, df.counts$target), sum)
    weight <- tapply(df.counts$prob, list(df.counts$source, df.counts$target), sum)

    # Update network with new matrices
    net$count <- count
    net$weight <- weight
  }

  # Replace NA values with 0
  net$weight[is.na(net$weight)] <- 0
  net$count[is.na(net$count)] <- 0

  # Return results
  if (return.object) {
    slot(object, slot.name) <- net
    return(object)
  } else {
    return(net)
  }
}

#' Get Interaction Counts for Visualization
#'
#' This function computes interaction counts (num.link) similar to those used in netVisual_scatter
#' for visualizing cell communication roles.
#'
#' @param object A CellChat object
#' @param signaling Vector of pathway names to include (NULL for all)
#' @param slot.name Slot name containing the network information, default is "netP"
#' @param thresh P-value threshold for significant interactions (default: 0.05)
#' @param remove.isolate Whether to remove isolated cell types with no interactions
#'
#' @return A named vector of interaction counts by cell type
#' @importFrom methods slot
#' @export
getInteractionCounts <- function(object, signaling = NULL, slot.name = "net",
                                 thresh = 0.05, remove.isolate = FALSE) {
  # Validate input
  if (!inherits(object, "CellChat")) {
    stop("Object must be a CellChat object")
  }

  # Extract network and validate
  if (is.null(signaling)) {
    message("Computing interaction counts from all signaling pathways")
  } else {
    message("Computing interaction counts from specified signaling pathways")
    if (!is.null(object@net$pathways)) {
      signaling <- signaling[signaling %in% object@net$pathways]
      if (length(signaling) == 0) {
        stop("No significant communication for the input signaling. All significant signaling are shown in `object@net$pathways`")
      }
    }
  }

  # Use aggregateNet to get interaction counts
  counts_matrix <- aggregateNetSimple(object, signaling = signaling,
                                      slot.name = slot.name, thresh = thresh,
                                      return.object = FALSE, remove.isolate = FALSE)$count

  # Calculate total connections for each cell type (incoming + outgoing - self)
  cell_types <- rownames(counts_matrix)
  counts <- rowSums(counts_matrix) + colSums(counts_matrix) - diag(counts_matrix)
  names(counts) <- cell_types

  # Remove isolated cell types if requested
  if (remove.isolate) {
    counts <- counts[counts > 0]
  }

  return(counts)
}


#' Calculate signaling scores for cell types using centrality metrics
#'
#' This function calculates signaling scores for specified cell types and pathways
#' based on the centrality metrics in a CellChat object.
#'
#' @param obj A CellChat object containing centrality metrics
#' @param cell_types Character vector of cell type names to analyze
#' @param pathways Character vector of pathway names to analyze
#' @param measure Character, specifies the signaling role to analyze: "sender" (outgoing signaling),
#'   "receiver" (incoming signaling), "both" (overall signaling), or "influence" (product of sender and receiver)
#' @param slot.name Character, name of the slot to extract data from (default: "netP")
#'
#' @return A named numeric vector of signaling scores for the specified cell types
#'
#' @examples
#' # Calculate sender scores for all cell types in the first pathway
#' sender_scores <- calculateSignaling(cellchat_obj,
#'                                    cell_types = levels(cellchat_obj@idents),
#'                                    pathways = cellchat_obj@netP$pathways[1],
#'                                    measure = "sender")
#'
#' @export
calculateSignaling <- function(obj, cell_types, pathways, measure = c("sender", "receiver", "both", "influence"), slot.name = "netP") {
  # Match argument
  measure <- match.arg(measure)

  # Check if centrality is computed
  if (length(slot(obj, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores!")
  }

  # Get centrality metrics
  centr <- slot(obj, slot.name)$centr

  # Initialize result vector
  result <- rep(0, length(cell_types))
  names(result) <- cell_types

  # Process each pathway
  for (pathway in pathways) {
    if (pathway %in% names(centr)) {
      # Get outdegree (sender) and indegree (receiver) for this pathway
      outdeg <- centr[[pathway]]$outdeg
      indeg <- centr[[pathway]]$indeg

      # For each cell type
      for (cell in cell_types) {
        # Get sender and receiver scores (0 if not present)
        sender_score <- if (cell %in% names(outdeg)) outdeg[cell] else 0
        receiver_score <- if (cell %in% names(indeg)) indeg[cell] else 0

        # Calculate scores based on measure type
        if (measure == "sender") {
          result[cell] <- result[cell] + sender_score
        } else if (measure == "receiver") {
          result[cell] <- result[cell] + receiver_score
        } else if (measure == "both") {
          result[cell] <- result[cell] + sender_score + receiver_score
        } else if (measure == "influence") {
          result[cell] <- result[cell] + (sender_score * receiver_score)
        }
      }
    } else {
      message(paste("Pathway", pathway, "not found in centrality metrics. Skipping."))
    }
  }

  return(result)
}

#' Calculate sender and receiver scores for cell types using centrality metrics
#'
#' This function calculates both sender and receiver signaling scores for specified
#' cell types and pathways based on the centrality metrics in a CellChat object.
#'
#' @param obj A CellChat object containing centrality metrics
#' @param cell_types Character vector of cell type names to analyze
#' @param pathways Character vector of pathway names to analyze
#' @param filter_vals Logical, whether to filter values based on significance (deprecated, kept for compatibility)
#' @param slot.name Character, name of the slot to extract data from (default: "netP")
#'
#' @return A list with two components:
#'   \itemize{
#'     \item sender: Named numeric vector of sender scores for each cell type
#'     \item receiver: Named numeric vector of receiver scores for each cell type
#'   }
#'
#' @examples
#' # Calculate sender and receiver scores for all cell types in all pathways
#' scores <- calculateScores(cellchat_obj,
#'                          cell_types = levels(cellchat_obj@idents),
#'                          pathways = cellchat_obj@netP$pathways)
#'
#' @export
calculateScores <- function(obj, cell_types, pathways, filter_vals = FALSE, slot.name = "netP") {
  # Check if centrality is computed
  if (length(slot(obj, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores!")
  }

  # Get centrality metrics
  centr <- slot(obj, slot.name)$centr

  # Initialize output vectors
  sender_scores <- rep(0, length(cell_types))
  receiver_scores <- rep(0, length(cell_types))
  names(sender_scores) <- names(receiver_scores) <- cell_types

  # Count pathways found for reporting
  pathways_found <- 0

  # Process each pathway
  for (pathway in pathways) {
    if (pathway %in% names(centr)) {
      pathways_found <- pathways_found + 1

      # Get outdegree (sender) and indegree (receiver) for this pathway
      outdeg <- centr[[pathway]]$outdeg
      indeg <- centr[[pathway]]$indeg

      # For each cell type, accumulate scores
      for (cell in cell_types) {
        # Add sender score if exists
        if (cell %in% names(outdeg)) {
          sender_scores[cell] <- sender_scores[cell] + outdeg[cell]
        }

        # Add receiver score if exists
        if (cell %in% names(indeg)) {
          receiver_scores[cell] <- receiver_scores[cell] + indeg[cell]
        }
      }
    }
  }

  # Report how many pathways were processed
  if (pathways_found < length(pathways)) {
    message(paste("Found centrality metrics for", pathways_found, "out of", length(pathways), "pathways."))
  }

  return(list(sender = sender_scores, receiver = receiver_scores))
}
