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

    # Check for required slots - use slotNames() for S4 objects
    required_slots <- c("net", "netP")
    obj_slots <- slotNames(obj)
    for (slot_name in required_slots) {
      if (!slot_name %in% obj_slots) {
        stop(sprintf("Object at index %d is missing required slot '%s'", i, slot_name),
             call. = FALSE)
      }
    }

    # Check for weight matrix in net slot
    net_slot <- methods::slot(obj, "net")
    if (is.null(net_slot) || !"weight" %in% names(net_slot)) {
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
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(cellchatlist)
#'
#' # Get common cell types across all objects
#' common_cells <- getCommonCellTypes(cellchatlist)
#'
#' # Get common cell types for specific conditions
#' common_cells <- getCommonCellTypes(cellchatlist, indices = c(1, 2))
#' }
#'
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

#' Get union of cell types across CellChat objects
#'
#' This function extracts all unique cell types present in any of the CellChat objects.
#'
#' @param object.list List of CellChat objects
#' @param indices Vector of indices to consider
#'
#' @return Vector of all unique cell type names
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(cellchatlist)
#'
#' # Get all unique cell types across all objects
#' all_cells <- getUnionCellTypes(cellchatlist)
#'
#' # Get all unique cell types for specific conditions
#' all_cells <- getUnionCellTypes(cellchatlist, indices = c(1, 2))
#'
#' # Compare with common cell types
#' common_cells <- getCommonCellTypes(cellchatlist)
#' union_cells <- getUnionCellTypes(cellchatlist)
#' lost_cells <- setdiff(union_cells, common_cells)
#' cat("Cell types lost with 'shared' strategy:", paste(lost_cells, collapse = ", "))
#' }
#'
#' @export
getUnionCellTypes <- function(object.list, indices = NULL) {
  if (is.null(indices)) {
    indices <- seq_along(object.list)
  }

  validateCellChatObjects(object.list, indices)

  cell_types <- lapply(object.list[indices], function(x) rownames(x@net$weight))
  union_cells <- Reduce(union, cell_types)

  return(union_cells)
}

#' Align cell types across CellChat objects with flexible strategy
#'
#' This function aligns cell types across multiple CellChat objects using either
#' "shared" (intersection) or "union" strategy. When using "union", missing cell
#' types in an object will be assigned zero values.
#'
#' @param object.list List of CellChat objects
#' @param indices Vector of indices to consider
#' @param strategy Character, either "shared" (intersection, default) or "union" (all cell types)
#'
#' @return List with components:
#'   \itemize{
#'     \item cell_types: Vector of aligned cell type names
#'     \item strategy: The strategy used
#'     \item missing_by_object: List showing which cell types are missing in each object
#'     \item n_missing: Count of how many objects miss each cell type
#'   }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(cellchatlist)
#'
#' # Align with shared strategy (default, intersection)
#' alignment <- alignCellTypes(cellchatlist,
#'                            indices = c(1, 2),
#'                            strategy = "shared")
#' print(alignment$cell_types)
#'
#' # Align with union strategy (all cell types)
#' alignment <- alignCellTypes(cellchatlist,
#'                            indices = c(1, 2),
#'                            strategy = "union")
#' print(alignment$cell_types)
#' print(alignment$missing_by_object)
#'
#' # Check which cell types would be lost with 'shared' strategy
#' alignment_shared <- alignCellTypes(cellchatlist, strategy = "shared")
#' alignment_union <- alignCellTypes(cellchatlist, strategy = "union")
#' lost_cells <- setdiff(alignment_union$cell_types, alignment_shared$cell_types)
#' if (length(lost_cells) > 0) {
#'   cat("Cell types lost with 'shared':", paste(lost_cells, collapse = ", "), "\n")
#' }
#'
#' # Use alignment results in analysis
#' alignment <- alignCellTypes(cellchatlist, strategy = "union")
#' result <- heatDiff(cellchatlist,
#'                   comparison = c(1, 2),
#'                   cell.type.strategy = "union")
#' }
#'
#' @export
alignCellTypes <- function(object.list, indices = NULL, strategy = c("shared", "union")) {
  strategy <- match.arg(strategy)

  if (is.null(indices)) {
    indices <- seq_along(object.list)
  }

  validateCellChatObjects(object.list, indices)

  # Get cell types from each object
  cell_types_list <- lapply(object.list[indices], function(x) rownames(x@net$weight))

  # Determine cell types based on strategy
  if (strategy == "shared") {
    aligned_cells <- Reduce(intersect, cell_types_list)

    if (length(aligned_cells) == 0) {
      stop("No common cell types found across objects. Consider using strategy='union'")
    }
  } else {
    # For union strategy, get all unique cell types preserving input order
    aligned_cells <- Reduce(union, cell_types_list)
  }

  # Track which cell types are missing in each object
  missing_by_object <- lapply(cell_types_list, function(ct) {
    setdiff(aligned_cells, ct)
  })
  names(missing_by_object) <- names(object.list)[indices]

  # Count how many objects are missing each cell type
  n_missing <- sapply(aligned_cells, function(cell) {
    sum(sapply(cell_types_list, function(ct) !cell %in% ct))
  })

  # Provide informative messages
  if (strategy == "union" && any(n_missing > 0)) {
    n_complete <- sum(n_missing == 0)
    n_partial <- sum(n_missing > 0 & n_missing < length(indices))
    n_sparse <- sum(n_missing >= length(indices) - 1)

    message(sprintf("Using union strategy: %d cell types aligned", length(aligned_cells)))
    message(sprintf("  - %d cell types present in all objects", n_complete))
    if (n_partial > 0) {
      message(sprintf("  - %d cell types present in some objects (will use 0 for missing)", n_partial))
    }
    if (n_sparse > 0) {
      message(sprintf("  - %d cell types present in only one object", n_sparse))
    }
  } else if (strategy == "shared") {
    n_excluded <- sum(sapply(cell_types_list, length)) - length(aligned_cells) * length(indices)
    message(sprintf("Using shared strategy: %d common cell types found", length(aligned_cells)))
    if (n_excluded > 0) {
      message(sprintf("  - %d cell types excluded (not present in all objects)",
                     length(unique(unlist(cell_types_list))) - length(aligned_cells)))
    }
  }

  return(list(
    cell_types = aligned_cells,
    strategy = strategy,
    missing_by_object = missing_by_object,
    n_missing = n_missing
  ))
}

#' Get union of pathways across CellChat objects
#'
#' This function extracts all unique pathways present in any of the CellChat objects.
#'
#' @param object.list List of CellChat objects
#' @param indices Vector of indices to consider
#'
#' @return Vector of all unique pathway names
#' @export
getUnionPathways <- function(object.list, indices = NULL) {
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

  union_pathways <- Reduce(union, pathways)

  return(union_pathways)
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

    # Check if any communications were found
    if (is.null(df.net) || nrow(df.net) == 0) {
      warning("No communications found for the specified pathways. Returning zero matrices.")
      # Get cell levels from object
      cells.level <- levels(object@idents)

      # Create empty matrices
      net$count <- matrix(0, nrow = length(cells.level), ncol = length(cells.level),
                         dimnames = list(cells.level, cells.level))
      net$weight <- matrix(0, nrow = length(cells.level), ncol = length(cells.level),
                          dimnames = list(cells.level, cells.level))
    } else {
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

    # Check if signaling pathways exist in netP slot (pathway-specific data)
    netP_pathways <- NULL
    if (!is.null(dimnames(object@netP$prob)[[3]])) {
      netP_pathways <- dimnames(object@netP$prob)[[3]]
    }

    # Check if pathways exist in net slot
    net_pathways <- NULL
    if (!is.null(object@net$pathways)) {
      net_pathways <- object@net$pathways
    }

    # Determine which slot to use based on pathway availability
    available_pathways <- NULL
    use_netP <- FALSE

    if (!is.null(netP_pathways)) {
      # Check if requested pathways are in netP
      signaling_in_netP <- signaling[signaling %in% netP_pathways]
      if (length(signaling_in_netP) > 0) {
        available_pathways <- netP_pathways
        use_netP <- TRUE
        slot.name <- "netP"  # Override to use netP slot
      }
    }

    if (!use_netP && !is.null(net_pathways)) {
      available_pathways <- net_pathways
    }

    # Filter to valid pathways
    if (!is.null(available_pathways)) {
      signaling_valid <- signaling[signaling %in% available_pathways]
      signaling_invalid <- signaling[!signaling %in% available_pathways]

      if (length(signaling_invalid) > 0) {
        warning(paste("Some pathways not found in object:",
                     paste(signaling_invalid, collapse = ", ")))
      }

      if (length(signaling_valid) == 0) {
        warning("No valid signaling pathways found. Using all available pathways.")
        signaling <- NULL
      } else {
        signaling <- signaling_valid
      }
    }
  }

  # Get counts matrix based on slot
  if (is.null(signaling)) {
    # Use aggregateNet for all pathways
    counts_matrix <- aggregateNetSimple(object, signaling = NULL,
                                        slot.name = slot.name, thresh = thresh,
                                        return.object = FALSE, remove.isolate = FALSE)$count
  } else if (slot.name == "netP") {
    # For netP slot, manually aggregate the specified pathways
    prob <- object@netP$prob
    pval <- object@netP$pval

    # Filter to specified pathways
    if (!is.null(signaling) && length(signaling) > 0) {
      pathway_indices <- which(dimnames(prob)[[3]] %in% signaling)
      if (length(pathway_indices) > 0) {
        prob <- prob[, , pathway_indices, drop = FALSE]
        pval <- pval[, , pathway_indices, drop = FALSE]
      } else {
        # No valid pathways found
        cell_types <- levels(object@idents)
        counts <- rep(0, length(cell_types))
        names(counts) <- cell_types
        warning("No valid pathways found in netP slot. Returning zero counts.")
        return(counts)
      }
    }

    # Apply threshold
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0

    # Aggregate across pathways
    counts_matrix <- apply(prob > 0, c(1, 2), sum)
  } else {
    # Use aggregateNet for net slot with pathway filtering
    counts_matrix <- aggregateNetSimple(object, signaling = signaling,
                                        slot.name = slot.name, thresh = thresh,
                                        return.object = FALSE, remove.isolate = FALSE)$count
  }

  # Check if counts_matrix is valid
  if (is.null(counts_matrix) || all(dim(counts_matrix) == 0)) {
    warning("No interaction counts found. Returning zero counts for all cell types.")
    # Return zero counts for all cell types in the object
    cell_types <- levels(object@idents)
    counts <- rep(0, length(cell_types))
    names(counts) <- cell_types
    return(counts)
  }

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
