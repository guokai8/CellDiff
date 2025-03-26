#' Create Integrated Network Visualization of Pathways, Cells, and Ligand-Receptor Pairs
#'
#' This function creates a network visualization that integrates pathways, cell types, and
#' ligand-receptor interactions, highlighting differences between two conditions.
#'
#' @param object.list List of CellChat objects
#' @param comparison Vector of indices to compare (length 2)
#' @param reference Index of the reference object (default is first in comparison)
#' @param pathways Vector of pathway names to include (if NULL, top differential pathways are selected)
#' @param top.LR Number of top ligand-receptor pairs to show, default is 5 per pathway
#' @param thresh Significance threshold for interactions, default is 0.05
#' @param pThresh P-value threshold for significant differences, default is 0.05
#' @param slot.name Slot name containing network information (default: "net")
#' @param color.cell Vector of colors for cell types
#' @param color.pathway Vector of colors for pathways
#' @param color.LR Vector of colors for ligand-receptor pairs
#' @param color.diff Colors for differential edges (inc/dec)
#' @param remove.isolate Whether to remove isolated nodes, default is TRUE
#' @param layout Layout algorithm to use, default is "fr" (Fruchterman-Reingold)
#' @param title Plot title
#' @param node.size.factor Scaling factor for node sizes
#' @param edge.width.factor Scaling factor for edge widths
#' @param label.size Size of labels
#' @param max.pathways Maximum number of pathways to show if pathways=NULL
#' @param sources.use Optional vector of source cell types to filter
#' @param targets.use Optional vector of target cell types to filter
#' @param node.label.repel Whether to use ggrepel for node labels to avoid overlaps
#' @param show.edge.arrows Whether to show arrows on edges
#' @param save.to.pdf File path to save PDF, or NULL to not save
#' @param pdf.width Width of PDF in inches if saving
#' @param pdf.height Height of PDF in inches if saving
#' @param use.all.pathways Whether to include pathways that only appear in one condition
#' @param show.cell.info Whether to show cell type information (like cluster, etc.)
#' @param cell.info.source Source of cell information (e.g., "idents", "meta", "labels")
#' @param cell.info.field Field name in meta.data to use for cell information
#' @param cell.color.by How to color cells ("type", "idents", "meta", "group")
#' @param debug Whether to print debug information
#'
#' @return A ggplot object showing the integrated network
#' @export
networkLRDiff <- function(object.list, comparison = c(1, 2), reference = NULL,
                          pathways = NULL, top.LR = 5, thresh = 0.05, pThresh = 0.05,
                          slot.name = "net", color.cell = NULL, color.pathway = NULL,
                          color.LR = NULL, color.diff = c("#FF4500", "#4169E1"),
                          remove.isolate = TRUE, layout = "fr",
                          title = "Differential LR Network",
                          node.size.factor = 1, edge.width.factor = 1,
                          label.size = 3, max.pathways = 5,
                          sources.use = NULL, targets.use = NULL,
                          node.label.repel = TRUE, show.edge.arrows = TRUE,
                          save.to.pdf = NULL, pdf.width = 12, pdf.height = 10,
                          use.all.pathways = FALSE, show.cell.info = FALSE,
                          cell.info.source = "idents", cell.info.field = NULL,
                          cell.color.by = "type", debug = FALSE) {

  if (debug) message("Starting networkLRDiff...")

  # Basic validation
  if (length(comparison) != 2) {
    stop("Please provide exactly 2 indices for comparison")
  }

  # Determine reference object
  if (is.null(reference)) {
    reference <- comparison[1]
    message("Using first object in comparison as reference")
  } else if (!(reference %in% comparison)) {
    warning("Reference index not in comparison. Using first object in comparison as reference.")
    reference <- comparison[1]
  }

  # Get condition names
  condition_names <- names(object.list)
  if (is.null(condition_names) || any(condition_names == "")) {
    condition_names <- paste0("Condition_", 1:length(object.list))
    names(object.list) <- condition_names
  }

  # Ensure reference is always first in comparison for direction consistency
  if (reference == comparison[2]) {
    comparison <- rev(comparison)
    message("Swapped comparison order to ensure reference is first")
  }

  if (debug) message("Finding pathways...")

  # If pathways not specified, either find differential pathways or get all pathways
  if (is.null(pathways)) {
    if (use.all.pathways) {
      # Get all pathways from both objects
      pathways1 <- unique(object.list[[comparison[1]]]@netP$pathways)
      pathways2 <- unique(object.list[[comparison[2]]]@netP$pathways)

      # Combine all unique pathways
      all_pathways <- unique(c(pathways1, pathways2))

      # Limit number of pathways if needed
      if (length(all_pathways) > max.pathways) {
        # Prioritize common pathways
        common_pathways <- intersect(pathways1, pathways2)
        unique_pathways <- setdiff(all_pathways, common_pathways)

        # Calculate weights for unique pathways
        unique_weights <- rep(0, length(unique_pathways))
        names(unique_weights) <- unique_pathways

        # Weight pathways in first object
        for (i in 1:length(unique_pathways)) {
          if (unique_pathways[i] %in% pathways1) {
            unique_weights[i] <- unique_weights[i] + 1
          }
        }

        # Weight pathways in second object
        for (i in 1:length(unique_pathways)) {
          if (unique_pathways[i] %in% pathways2) {
            unique_weights[i] <- unique_weights[i] + 1
          }
        }

        # Sort unique pathways by weight
        sorted_unique <- names(sort(unique_weights, decreasing = TRUE))

        # Combine common and top unique pathways
        pathways <- c(common_pathways, sorted_unique)
        pathways <- pathways[1:min(max.pathways, length(pathways))]
      } else {
        pathways <- all_pathways
      }

      message(paste("Using all", length(pathways), "pathways from both conditions"))
    } else {
      message("No pathways specified, finding top differential pathways...")

      # Use rankDiff to get significant pathways
      tryCatch({
        pathway_results <- rankDiff(
          object.list = object.list,
          comparison = comparison,
          pThresh = pThresh,
          slot.name = slot.name,
          show.pval = FALSE,
          show.heatmap = FALSE,
          sources.use = sources.use,
          targets.use = targets.use
        )

        if (is.null(pathway_results) || length(pathway_results$significant_paths) == 0) {
          stop("No significant pathways found. Try increasing pThresh or specify pathways manually.")
        }

        # Use top pathways
        num_pathways <- min(max.pathways, length(pathway_results$significant_paths))
        pathways <- pathway_results$significant_paths[1:num_pathways]

        message(paste("Selected top", num_pathways, "differential pathways:",
                      paste(pathways, collapse = ", ")))
      }, error = function(e) {
        stop("Error finding differential pathways: ", e$message)
      })
    }
  }

  # Extract cell types from both objects
  all_cell_types <- union(
    rownames(slot(object.list[[comparison[1]]], slot.name)$weight),
    rownames(slot(object.list[[comparison[2]]], slot.name)$weight)
  )

  if (debug) message("Found ", length(all_cell_types), " total cell types across both conditions")

  # Apply source/target filters if provided
  if (!is.null(sources.use)) {
    if (is.character(sources.use)) {
      sources.use <- intersect(sources.use, all_cell_types)
      if (length(sources.use) == 0) {
        stop("None of the specified source cell types found in the objects")
      }
    }
  } else {
    sources.use <- all_cell_types
  }

  if (!is.null(targets.use)) {
    if (is.character(targets.use)) {
      targets.use <- intersect(targets.use, all_cell_types)
      if (length(targets.use) == 0) {
        stop("None of the specified target cell types found in the objects")
      }
    }
  } else {
    targets.use <- all_cell_types
  }

  cellTypes <- union(sources.use, targets.use)

  # Get cell type information if requested
  cell_info <- data.frame(
    id = cellTypes,
    label = cellTypes,
    type = "cell",
    group = NA,
    info = NA,
    stringsAsFactors = FALSE
  )

  if (show.cell.info) {
    for (i in comparison) {
      if (cell.info.source == "idents") {
        # Get cell identities
        cells_in_obj <- intersect(cellTypes, rownames(slot(object.list[[i]], slot.name)$weight))
        for (cell in cells_in_obj) {
          idx <- which(cell_info$id == cell)
          if (length(idx) > 0) {
            if (is.na(cell_info$group[idx])) {
              cell_info$group[idx] <- condition_names[i]
            } else {
              cell_info$group[idx] <- paste(cell_info$group[idx], condition_names[i], sep = "/")
            }

            # Update label with ident information
            cell_label <- cell
            if (is.na(cell_info$info[idx])) {
              cell_info$info[idx] <- cell
            }
            cell_info$label[idx] <- paste0(cell, " (", cell_info$group[idx], ")")
          }
        }
      } else if (cell.info.source == "meta" && !is.null(cell.info.field)) {
        # Try to get information from meta.data if available
        if ("meta.data" %in% slotNames(object.list[[i]]) &&
            cell.info.field %in% colnames(object.list[[i]]@meta.data)) {

          meta_data <- object.list[[i]]@meta.data
          cells_in_obj <- intersect(cellTypes, rownames(meta_data))

          for (cell in cells_in_obj) {
            idx <- which(cell_info$id == cell)
            if (length(idx) > 0) {
              if (is.na(cell_info$group[idx])) {
                cell_info$group[idx] <- condition_names[i]
              } else {
                cell_info$group[idx] <- paste(cell_info$group[idx], condition_names[i], sep = "/")
              }

              # Update info with meta data
              cell_meta <- meta_data[cell, cell.info.field]
              if (is.na(cell_info$info[idx])) {
                cell_info$info[idx] <- cell_meta
              } else {
                cell_info$info[idx] <- paste(cell_info$info[idx], cell_meta, sep = "/")
              }

              # Update label
              cell_info$label[idx] <- paste0(cell, " (", cell_info$info[idx], ")")
            }
          }
        }
      }
    }
  }

  # Initialize node and edge data frames with predefined structure
  nodes <- data.frame(
    id = character(),
    label = character(),
    type = character(),
    pathway = character(),
    size = numeric(),
    group = character(),
    info = character(),
    stringsAsFactors = FALSE
  )

  edges <- data.frame(
    source = character(),
    target = character(),
    weight = numeric(),
    direction = character(),
    type = character(),
    significance = numeric(),
    stringsAsFactors = FALSE
  )

  # Add cell nodes
  for (i in 1:nrow(cell_info)) {
    nodes <- rbind(nodes, data.frame(
      id = cell_info$id[i],
      label = cell_info$label[i],
      type = "cell",
      pathway = NA,
      size = 0.8 * node.size.factor,
      group = cell_info$group[i],
      info = cell_info$info[i],
      stringsAsFactors = FALSE
    ))
  }

  # Gather pathway information
  pathway_info <- data.frame(
    id = pathways,
    label = pathways,
    in_obj1 = FALSE,
    in_obj2 = FALSE,
    stringsAsFactors = FALSE
  )

  # Check which pathways exist in which objects
  for (i in 1:nrow(pathway_info)) {
    if (pathway_info$id[i] %in% object.list[[comparison[1]]]@netP$pathways) {
      pathway_info$in_obj1[i] <- TRUE
    }
    if (pathway_info$id[i] %in% object.list[[comparison[2]]]@netP$pathways) {
      pathway_info$in_obj2[i] <- TRUE
    }
  }

  # Add pathway nodes with presence information
  for (i in 1:nrow(pathway_info)) {
    pathway <- pathway_info$id[i]
    if (pathway_info$in_obj1[i] && pathway_info$in_obj2[i]) {
      pathway_group <- "both"
    } else if (pathway_info$in_obj1[i]) {
      pathway_group <- condition_names[comparison[1]]
    } else {
      pathway_group <- condition_names[comparison[2]]
    }

    nodes <- rbind(nodes, data.frame(
      id = pathway,
      label = paste0(pathway, " (", pathway_group, ")"),
      type = "pathway",
      pathway = pathway,
      size = 1.5 * node.size.factor,
      group = pathway_group,
      info = pathway_group,
      stringsAsFactors = FALSE
    ))
  }

  # Extract ligand-receptor pairs for all pathways
  lr_data <- list()
  for (pathway in pathways) {
    if (debug) message("Processing pathway: ", pathway)

    pathway_lr_pairs <- list()

    # Extract LR pairs for pathway from first object if it exists there
    if (pathway %in% object.list[[comparison[1]]]@netP$pathways) {
      tryCatch({
        lr_pairs1 <- CellChat::searchPair(
          signaling = pathway,
          pairLR.use = object.list[[comparison[1]]]@LR$LRsig,
          key = "pathway_name",
          matching.exact = TRUE,
          pair.only = FALSE
        )

        if (nrow(lr_pairs1) > 0) {
          for (lr_pair in rownames(lr_pairs1)) {
            if (!(lr_pair %in% names(pathway_lr_pairs))) {
              pathway_lr_pairs[[lr_pair]] <- list(
                in_obj1 = TRUE,
                in_obj2 = FALSE,
                diff_score = 0
              )
            }
          }
        }
      }, error = function(e) {
        message(paste("Error getting LR pairs for pathway", pathway, "in object 1:", e$message))
      })
    }

    # Extract LR pairs for pathway from second object if it exists there
    if (pathway %in% object.list[[comparison[2]]]@netP$pathways) {
      tryCatch({
        lr_pairs2 <- CellChat::searchPair(
          signaling = pathway,
          pairLR.use = object.list[[comparison[2]]]@LR$LRsig,
          key = "pathway_name",
          matching.exact = TRUE,
          pair.only = FALSE
        )

        if (nrow(lr_pairs2) > 0) {
          for (lr_pair in rownames(lr_pairs2)) {
            if (!(lr_pair %in% names(pathway_lr_pairs))) {
              pathway_lr_pairs[[lr_pair]] <- list(
                in_obj1 = FALSE,
                in_obj2 = TRUE,
                diff_score = 0
              )
            } else {
              pathway_lr_pairs[[lr_pair]]$in_obj2 <- TRUE
            }
          }
        }
      }, error = function(e) {
        message(paste("Error getting LR pairs for pathway", pathway, "in object 2:", e$message))
      })
    }

    # Calculate differential scores for common LR pairs
    for (lr_pair in names(pathway_lr_pairs)) {
      if (pathway_lr_pairs[[lr_pair]]$in_obj1 && pathway_lr_pairs[[lr_pair]]$in_obj2) {
        if (lr_pair %in% dimnames(slot(object.list[[comparison[1]]], slot.name)$prob)[[3]] &&
            lr_pair %in% dimnames(slot(object.list[[comparison[2]]], slot.name)$prob)[[3]]) {

          # Extract probability matrices
          prob1 <- slot(object.list[[comparison[1]]], slot.name)$prob[,,lr_pair]
          prob2 <- slot(object.list[[comparison[2]]], slot.name)$prob[,,lr_pair]

          # Apply significance threshold
          prob1[slot(object.list[[comparison[1]]], slot.name)$pval[,,lr_pair] > thresh] <- 0
          prob2[slot(object.list[[comparison[2]]], slot.name)$pval[,,lr_pair] > thresh] <- 0

          # Find common sources and targets
          common_sources <- intersect(sources.use, union(rownames(prob1), rownames(prob2)))
          common_targets <- intersect(targets.use, union(colnames(prob1), colnames(prob2)))

          if (length(common_sources) > 0 && length(common_targets) > 0) {
            # Ensure both matrices have the same dimensions
            full_prob1 <- matrix(0, nrow = length(common_sources), ncol = length(common_targets))
            full_prob2 <- matrix(0, nrow = length(common_sources), ncol = length(common_targets))
            rownames(full_prob1) <- rownames(full_prob2) <- common_sources
            colnames(full_prob1) <- colnames(full_prob2) <- common_targets

            # Fill in values from the original matrices
            for (src in intersect(common_sources, rownames(prob1))) {
              for (tgt in intersect(common_targets, colnames(prob1))) {
                full_prob1[src, tgt] <- prob1[src, tgt]
              }
            }

            for (src in intersect(common_sources, rownames(prob2))) {
              for (tgt in intersect(common_targets, colnames(prob2))) {
                full_prob2[src, tgt] <- prob2[src, tgt]
              }
            }

            # Calculate difference
            prob_diff <- abs(full_prob2 - full_prob1)

            # Store maximum difference as the score
            pathway_lr_pairs[[lr_pair]]$diff_score <- max(prob_diff)
          }
        }
      } else {
        # For LR pairs only in one condition, assign a moderate score
        pathway_lr_pairs[[lr_pair]]$diff_score <- 0.5
      }
    }

    # Convert to a more usable format and sort by diff_score
    lr_scores <- sapply(pathway_lr_pairs, function(x) x$diff_score)
    lr_presence <- lapply(pathway_lr_pairs, function(x) c(x$in_obj1, x$in_obj2))

    # Select top LR pairs
    if (length(lr_scores) > 0) {
      top_indices <- order(lr_scores, decreasing = TRUE)[1:min(top.LR, length(lr_scores))]
      top_lr_pairs <- names(lr_scores)[top_indices]

      # Store pathway, LR pairs, and presence information
      lr_data[[pathway]] <- list(
        pairs = top_lr_pairs,
        presence = lr_presence[top_lr_pairs]
      )

      if (debug) {
        message(paste("Selected", length(top_lr_pairs), "top LR pairs for pathway:", pathway))
      }
    }
  }

  # Check if we have any data
  if (length(lr_data) == 0) {
    stop("No ligand-receptor data found for the specified pathways.")
  }

  # Add LR nodes with unique IDs and presence information
  for (pathway in names(lr_data)) {
    for (i in 1:length(lr_data[[pathway]]$pairs)) {
      lr_pair <- lr_data[[pathway]]$pairs[i]
      presence <- lr_data[[pathway]]$presence[[i]]

      # Create a unique ID for this LR pair in this pathway
      lr_id <- paste0(pathway, "_", lr_pair)

      # Split LR pair into ligand and receptor
      lr_components <- strsplit(lr_pair, "_")[[1]]
      if (length(lr_components) >= 2) {
        ligand <- lr_components[1]
        receptor <- lr_components[2]
        lr_label <- paste(ligand, receptor, sep = "-")
      } else {
        lr_label <- lr_pair
      }

      # Determine which conditions this LR pair is present in
      if (presence[1] && presence[2]) {
        lr_group <- "both"
      } else if (presence[1]) {
        lr_group <- condition_names[comparison[1]]
      } else {
        lr_group <- condition_names[comparison[2]]
      }

      # Add information to the label
      lr_label <- paste0(lr_label, " (", lr_group, ")")

      # Add LR node
      nodes <- rbind(nodes, data.frame(
        id = lr_id,
        label = lr_label,
        type = "lr",
        pathway = pathway,
        size = 1.2 * node.size.factor,
        group = lr_group,
        info = lr_group,
        stringsAsFactors = FALSE
      ))

      # Add edge from pathway to LR
      edges <- rbind(edges, data.frame(
        source = pathway,
        target = lr_id,
        weight = 1,
        direction = "none",
        type = "pathway_lr",
        significance = 1,  # Always show these edges
        stringsAsFactors = FALSE
      ))
    }
  }

  # Add cell-to-LR and LR-to-cell edges for all LR pairs
  for (pathway in names(lr_data)) {
    for (i in 1:length(lr_data[[pathway]]$pairs)) {
      lr_pair <- lr_data[[pathway]]$pairs[i]
      presence <- lr_data[[pathway]]$presence[[i]]

      # Get the unique ID for this LR pair
      lr_id <- paste0(pathway, "_", lr_pair)

      # Process interactions in first object
      if (presence[1] &&
          lr_pair %in% dimnames(slot(object.list[[comparison[1]]], slot.name)$prob)[[3]]) {

        # Extract probability matrix and p-values
        prob1 <- slot(object.list[[comparison[1]]], slot.name)$prob[,,lr_pair]
        pval1 <- slot(object.list[[comparison[1]]], slot.name)$pval[,,lr_pair]

        # Apply threshold
        prob1[pval1 > thresh] <- 0

        # Filter for source and target cells
        valid_sources <- intersect(rownames(prob1), sources.use)
        valid_targets <- intersect(colnames(prob1), targets.use)

        # Add edges for significant interactions
        for (src in valid_sources) {
          for (tgt in valid_targets) {
            if (prob1[src, tgt] > 0) {
              # Add edge from source cell to LR
              edges <- rbind(edges, data.frame(
                source = src,
                target = lr_id,
                weight = prob1[src, tgt] * edge.width.factor * 0.2,
                direction = "ref",  # Indicating interaction in reference
                type = "cell_lr",
                significance = 1,
                stringsAsFactors = FALSE
              ))

              # Add edge from LR to target cell
              edges <- rbind(edges, data.frame(
                source = lr_id,
                target = tgt,
                weight = prob1[src, tgt] * edge.width.factor *0.2,
                direction = "ref",
                type = "lr_cell",
                significance = 1,
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }

      # Process interactions in second object
      if (presence[2] &&
          lr_pair %in% dimnames(slot(object.list[[comparison[2]]], slot.name)$prob)[[3]]) {

        # Extract probability matrix and p-values
        prob2 <- slot(object.list[[comparison[2]]], slot.name)$prob[,,lr_pair]
        pval2 <- slot(object.list[[comparison[2]]], slot.name)$pval[,,lr_pair]

        # Apply threshold
        prob2[pval2 > thresh] <- 0

        # Filter for source and target cells
        valid_sources <- intersect(rownames(prob2), sources.use)
        valid_targets <- intersect(colnames(prob2), targets.use)

        # Add edges for significant interactions
        for (src in valid_sources) {
          for (tgt in valid_targets) {
            if (prob2[src, tgt] > 0) {
              # Check if there's already an edge from the reference
              if (presence[1]) {
                # Find if there's a corresponding interaction in reference
                ref_value <- 0
                if (src %in% rownames(prob1) && tgt %in% colnames(prob1)) {
                  ref_value <- prob1[src, tgt]
                }

                # Calculate difference
                diff_value <- prob2[src, tgt] - ref_value

                if (abs(diff_value) > 0.01) {  # Small threshold to avoid tiny differences
                  # Add differential edges
                  edges <- rbind(edges, data.frame(
                    source = src,
                    target = lr_id,
                    weight = abs(diff_value) * 0.4 * edge.width.factor,
                    direction = ifelse(diff_value > 0, "increase", "decrease"),
                    type = "cell_lr_diff",
                    significance = min(1, abs(diff_value) * 10),  # Scale significance by difference
                    stringsAsFactors = FALSE
                  ))

                  edges <- rbind(edges, data.frame(
                    source = lr_id,
                    target = tgt,
                    weight = abs(diff_value) * 0.4 * edge.width.factor,
                    direction = ifelse(diff_value > 0, "increase", "decrease"),
                    type = "lr_cell_diff",
                    significance = min(1, abs(diff_value) * 10),
                    stringsAsFactors = FALSE
                  ))
                }
              } else {
                # If not in reference, add as a new interaction
                edges <- rbind(edges, data.frame(
                  source = src,
                  target = lr_id,
                  weight = prob2[src, tgt] * edge.width.factor,
                  direction = "new",  # Indicating new interaction
                  type = "cell_lr",
                  significance = 1,
                  stringsAsFactors = FALSE
                ))

                edges <- rbind(edges, data.frame(
                  source = lr_id,
                  target = tgt,
                  weight = prob2[src, tgt] * edge.width.factor,
                  direction = "new",
                  type = "lr_cell",
                  significance = 1,
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
        }
      }
    }
  }

  # Remove isolated nodes if requested
  if (remove.isolate) {
    connected_nodes <- unique(c(edges$source, edges$target))
    nodes <- nodes[nodes$id %in% connected_nodes, ]
    if (debug) message("After removing isolates: ", nrow(nodes), " nodes remain")
  }

  # Generate network layout
  if (nrow(edges) == 0 || nrow(nodes) == 0) {
    stop("No significant interactions found for visualization.")
  }

  if (debug) message("Generating network layout...")

  # Create igraph network and layout
  network <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = TRUE
  )

  # Choose layout algorithm
  if (layout == "fr") {
    layout_coords <- igraph::layout_with_fr(network)
  } else if (layout == "circle") {
    layout_coords <- igraph::layout_in_circle(network)
  } else {
    layout_coords <- igraph::layout_nicely(network)
  }

  # Convert layout to data frame
  layout_df <- data.frame(
    id = igraph::V(network)$name,
    x = layout_coords[,1],
    y = layout_coords[,2],
    stringsAsFactors = FALSE
  )

  # Join layout with node attributes
  node_data <- merge(nodes, layout_df, by = "id")

  # Join edge data with node positions
  edge_data <- merge(edges, layout_df, by.x = "source", by.y = "id")
  names(edge_data)[names(edge_data) == "x"] <- "x_source"
  names(edge_data)[names(edge_data) == "y"] <- "y_source"
  edge_data <- merge(edge_data, layout_df, by.x = "target", by.y = "id")
  names(edge_data)[names(edge_data) == "x"] <- "x_target"
  names(edge_data)[names(edge_data) == "y"] <- "y_target"

  # Prepare node colors based on selected coloring scheme
  if (cell.color.by == "type") {
    # Color by node type (default)
    if (is.null(color.cell)) {
      color.cell <- "skyblue"
    }
    if (is.null(color.pathway)) {
      color.pathway <- "orange"
    }
    if (is.null(color.LR)) {
      color.LR <- "lightgreen"
    }

    node_fill_values <- c(
      "cell" = color.cell,
      "pathway" = color.pathway,
      "lr" = color.LR
    )
  } else if (cell.color.by == "group") {
    # Color by group (condition presence)
    # Create a color palette for conditions
    cond_colors <- c("both" = "purple")
    cond_colors[condition_names[comparison[1]]] <- "darkcyan"
    cond_colors[condition_names[comparison[2]]] <- "darkorange"

    # If global_colors exists, use it
    if (exists("global_colors") && length(global_colors) >= 3) {
      cond_colors <- c("both" = global_colors[1])
      cond_colors[condition_names[comparison[1]]] <- global_colors[2]
      cond_colors[condition_names[comparison[2]]] <- global_colors[3]
    }

    # Assign colors based on group
    node_fill_values <- cond_colors
  } else if (cell.color.by == "idents" || cell.color.by == "meta") {
    # Color by cell identity or metadata
    unique_info <- unique(na.omit(node_data$info))

    if (exists("global_colors")) {
      info_colors <- global_colors[1:length(unique_info)]
    } else {
      info_colors <- scales::hue_pal()(length(unique_info))
    }
    names(info_colors) <- unique_info

    # Manually create a color vector for each node
    node_colors <- rep(NA, nrow(node_data))

    for (i in 1:nrow(node_data)) {
      if (!is.na(node_data$info[i])) {
        # For cell types with info, use info-based color
        node_colors[i] <- info_colors[node_data$info[i]]
      } else {
        # For other node types, use default colors
        if (node_data$type[i] == "pathway") {
          node_colors[i] <- ifelse(is.null(color.pathway), "orange", color.pathway)
        } else if (node_data$type[i] == "lr") {
          node_colors[i] <- ifelse(is.null(color.LR), "lightgreen", color.LR)
        } else {
          node_colors[i] <- ifelse(is.null(color.cell), "skyblue", color.cell)
        }
      }
    }

    # Create dummy values to set up the fill scale
    node_fill_values <- info_colors
  } else {
    # Default coloring
    node_fill_values <- c(
      "cell" = "skyblue",
      "pathway" = "orange",
      "lr" = "lightgreen"
    )
  }

  # Update title to include comparison conditions
  if (title == "Differential LR Network") {
    title <- paste0("Differential LR Network: ",
                    condition_names[comparison[2]], " vs ",
                    condition_names[comparison[1]])
  }

  # Set up enhanced edge colors for different types of interactions
  edge_colors <- c(
    "increase" = color.diff[1],    # Increased from reference
    "decrease" = color.diff[2],    # Decreased from reference
    "none" = "grey50",             # Connection without direction
    "ref" = "grey30",              # Present in reference
    "new" = "green4"               # New in comparison condition
  )

  # Create the plot
  p <- ggplot2::ggplot() +
    # Add edges
    ggplot2::geom_segment(
      data = edge_data,
      ggplot2::aes(
        x = x_source, y = y_source,
        xend = x_target, yend = y_target,
        color = direction,
        size = weight,
        alpha = significance *0.5
      ),
      arrow = if(show.edge.arrows) ggplot2::arrow(length = ggplot2::unit(0.1, "inches"), type = "closed") else NULL
    )

  # Add nodes (different approach depending on coloring scheme)
  if (cell.color.by == "type" || cell.color.by == "group") {
    # Use fill aesthetic with categories
    p <- p + ggplot2::geom_point(
      data = node_data,
      ggplot2::aes(
        x = x, y = y,
        size = size,
        fill = if(cell.color.by == "type") type else group,
        shape = type
      ),
      color = "black",
      stroke = 0.5
    )
  } else {
    # Use individual colors directly
    p <- p + ggplot2::geom_point(
      data = node_data,
      ggplot2::aes(
        x = x, y = y,
        size = size,
        shape = type
      ),
      fill = node_colors,
      color = "black",
      stroke = 0.5
    )
  }

  # Add labels - either regular or repelled
  if (node.label.repel && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = node_data,
      ggplot2::aes(
        x = x, y = y,
        label = label
      ),
      size = label.size,
      box.padding = 0.5,
      point.padding = 0.2,
      force = 1,
      segment.color = "grey50",
      max.overlaps = 20
    )
  } else {
    p <- p + ggplot2::geom_text(
      data = node_data,
      ggplot2::aes(
        x = x, y = y,
        label = label
      ),
      size = label.size,
      hjust = 0.5, vjust = -1.2
    )
  }

  # Add scales based on coloring scheme
  p <- p + ggplot2::scale_color_manual(
    values = edge_colors,
    name = "Direction",
    labels = c("Increased", "Decreased", "Connection", "Reference", "New")
  )

  if (cell.color.by == "type") {
    p <- p + ggplot2::scale_fill_manual(
      values = node_fill_values,
      name = "Node Type",
      labels = c("Cell Type", "Pathway", "LR Pair")
    )
  } else if (cell.color.by == "group") {
    p <- p + ggplot2::scale_fill_manual(
      values = node_fill_values,
      name = "Presence",
      labels = names(node_fill_values)
    )
  } else if (cell.color.by == "idents" || cell.color.by == "meta") {
    # For idents or meta coloring, we already handled the colors directly
    p <- p + ggplot2::scale_fill_manual(
      values = node_fill_values,
      name = "Cell Info",
      guide = "legend"
    )
  }

  # Finalize the plot
  p <- p +
    ggplot2::scale_shape_manual(
      values = c("cell" = 21, "pathway" = 22, "lr" = 23),
      name = "Node Type",
      labels = c("Cell Type", "Pathway", "LR Pair")
    ) +
    ggplot2::scale_size(range = c(1, 5), guide = "none") +
    ggplot2::scale_alpha(range = c(0.1, 0.6), guide = "none") +
    # Set theme
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    ggplot2::ggtitle(title)

  # Save to PDF if requested
  if (!is.null(save.to.pdf)) {
    tryCatch({
      grDevices::pdf(save.to.pdf, width = pdf.width, height = pdf.height)
      print(p)
      grDevices::dev.off()
      message("Plot saved to ", save.to.pdf)
    }, error = function(e) {
      warning("Failed to save PDF: ", e$message)
    })
  }

  if (debug) message("Plot created successfully")

  return(p)
}
