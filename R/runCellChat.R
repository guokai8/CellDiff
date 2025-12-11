#' Create CellChat Objects from Seurat Object
#'
#' @description
#' A wrapper function that takes a Seurat object and creates CellChat objects for each
#' condition. This simplifies the CellChat workflow by automating data extraction,
#' preprocessing, and communication inference for multiple conditions.
#'
#' @param seurat_object A Seurat object containing single-cell RNA-seq data
#' @param group.by Column name in metadata to split conditions (e.g., "condition", "treatment", "genotype")
#' @param conditions Vector of condition names to include (NULL for all conditions in group.by)
#' @param species Species for CellChatDB: "human" (default) or "mouse"
#' @param min.cells Minimum number of cells required per cell type for filtering (default: 10)
#' @param cell.type.column Column name in metadata containing cell type annotations (default: "cell_type")
#'
#' @section CellChat Parameters:
#' @param population.size Logical. Whether to consider cell population size in communication strength (default: TRUE)
#' @param type Method for computing communication probability: "triMean" (default), "truncatedMean", "thresholdedMean", or "median"
#' @param trim Trim value for truncatedMean method (default: 0.1)
#' @param raw.use Logical. Whether to use raw data for communication probability (default: TRUE)
#' @param thresh Threshold for mean expression to define expressed genes in identifyOverExpressedGenes (default: 0.05)
#' @param thresh.centrality Threshold for computing centrality in netAnalysis_computeCentrality (default: 0.05)
#' @param compute.centrality Logical. Whether to compute network centrality metrics (default: TRUE)
#'
#' @section Differential Analysis Parameters:
#' @param run.analysis Logical. If TRUE, automatically runs differential analysis after creating CellChat objects (default: FALSE)
#' @param reference Reference condition for differential analysis (required if run.analysis = TRUE)
#' @param comparison_method Comparison method: "all_vs_ref" (default), "all_vs_all", or "custom_pairs"
#' @param custom_comparisons List of custom comparison pairs (for custom_pairs method)
#' @param cell.type.strategy Strategy for handling missing cell types: "shared" (default) or "union"
#' @param show.all Logical. If TRUE, shows all pathways with visual indicators (default: FALSE)
#' @param pThresh P-value threshold for significance (default: 0.05)
#' @param use_log2fc Logical. Use log2 fold change (default: TRUE)
#' @param show_plots Character vector specifying which analyses to run when run.analysis = TRUE.
#'   Options: "barplot" (rankDiffM), "heatmap" (heatDiffM), "slope" (rankDiffM slope plot),
#'   "scatter" (scatterDiff2DM), "network" (networkLRDiff).
#'   Default: c("barplot", "heatmap"). Use NULL to run only rankDiffM without plots.
#' @param top.n Number of top pathways to highlight (default: NULL for all)
#' @param verbose Logical. If TRUE, prints progress messages (default: TRUE)
#'
#' @return If run.analysis = FALSE: A named list of CellChat objects, one for each condition
#' @return If run.analysis = TRUE: A CellDiffAnalysis object with CellChat objects and analysis results
#'
#' @details
#' This wrapper function automates the CellChat workflow:
#' \itemize{
#'   \item Validates Seurat object and splits by condition
#'   \item Creates CellChat objects for each condition
#'   \item Runs CellChat preprocessing, communication inference, and aggregation
#'   \item Returns a named list of processed CellChat objects
#' }
#'
#' After creating CellChat objects, you can use CellDiff functions for differential analysis:
#' - Use \code{compareCell()} for comprehensive differential analysis
#' - Use \code{rankDiffM()} for pathway ranking
#' - Use \code{heatDiffM()} for heatmaps
#' - And other CellDiff functions as needed
#'
#' @examples
#' \dontrun{
#' # Option 1: Just create CellChat objects (run.analysis = FALSE, default)
#' cellchat_list <- runCellChat(
#'   seurat_object = pbmc,
#'   group.by = "condition",
#'   species = "human"
#' )
#'
#' # Then manually run analysis functions
#' rankDiffM(object.list = cellchat_list, reference = "WT")
#' heatDiffM(object.list = cellchat_list, reference = "WT")
#'
#' # Option 2: Run everything automatically (run.analysis = TRUE)
#' results <- runCellChat(
#'   seurat_object = pbmc,
#'   group.by = "condition",
#'   species = "human",
#'   run.analysis = TRUE,
#'   reference = "WT",
#'   show_plots = c("barplot", "heatmap", "scatter")  # Choose which analyses
#' )
#'
#' # Access CellChat objects and results
#' results$cellchat_objects[["WT"]]
#' results$pathway_analysis$comparison_barplot
#' results$heatmap  # Only if "heatmap" in show_plots
#' results$scatter  # Only if "scatter" in show_plots
#' results$summary
#'
#' # Run only specific analyses
#' results <- runCellChat(
#'   seurat_object = pbmc,
#'   group.by = "condition",
#'   species = "human",
#'   run.analysis = TRUE,
#'   reference = "WT",
#'   show_plots = c("barplot")  # Only pathway ranking
#' )
#'
#' # Run all available analyses
#' results <- runCellChat(
#'   seurat_object = pbmc,
#'   group.by = "condition",
#'   species = "human",
#'   run.analysis = TRUE,
#'   reference = "WT",
#'   show_plots = c("barplot", "heatmap", "slope", "scatter", "network")
#' )
#'
#' # Customize CellChat parameters
#' results <- runCellChat(
#'   seurat_object = pbmc,
#'   group.by = "condition",
#'   species = "human",
#'   type = "truncatedMean",        # Use truncated mean instead of triMean
#'   trim = 0.2,                    # Increase trim value
#'   thresh = 0.1,                  # Higher threshold for expressed genes
#'   thresh.centrality = 0.1,       # Centrality threshold
#'   population.size = FALSE,       # Don't consider population size
#'   compute.centrality = TRUE,     # Compute centrality metrics
#'   min.cells = 5                  # Lower minimum cell requirement
#' )
#' }
#'
#' @export
runCellChat <- function(seurat_object,
                        group.by,
                        conditions = NULL,
                        species = c("human", "mouse"),
                        min.cells = 10,
                        cell.type.column = "cell_type",
                        # CellChat parameters
                        population.size = TRUE,
                        type = c("triMean", "truncatedMean", "thresholdedMean", "median"),
                        trim = 0.1,
                        raw.use = TRUE,
                        thresh = 0.05,
                        thresh.centrality = 0.05,
                        compute.centrality = TRUE,
                        # Differential analysis parameters
                        run.analysis = FALSE,
                        reference = NULL,
                        comparison_method = "all_vs_ref",
                        custom_comparisons = NULL,
                        cell.type.strategy = c("shared", "union"),
                        show.all = FALSE,
                        pThresh = 0.05,
                        use_log2fc = TRUE,
                        show_plots = c("barplot", "heatmap"),
                        top.n = NULL,
                        verbose = TRUE) {

  # Match arguments
  species <- match.arg(species)
  type <- match.arg(type)
  cell.type.strategy <- match.arg(cell.type.strategy)

  # Validate run.analysis parameter
  if (run.analysis && is.null(reference) && comparison_method == "all_vs_ref") {
    stop("reference must be specified when run.analysis = TRUE and comparison_method = 'all_vs_ref'")
  }

  # Load CellChat
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("CellChat package is required. Install it with: devtools::install_github('sqjin/CellChat')")
  }

  # Validate Seurat object
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }

  # Check Seurat version compatibility
  version_info <- checkSeuratVersion(seurat_object, verbose = verbose)

  if (verbose) {
    cat("\n=======================================================\n")
    cat("Creating CellChat Objects from Seurat\n")
    cat("=======================================================\n\n")
  }

  # Check group.by column exists
  if (!group.by %in% colnames(seurat_object@meta.data)) {
    stop("group.by column '", group.by, "' not found in Seurat object metadata")
  }

  # Check cell type column exists
  if (!cell.type.column %in% colnames(seurat_object@meta.data)) {
    stop("cell.type.column '", cell.type.column, "' not found in Seurat object metadata")
  }

  # Get available conditions
  available_conditions <- unique(seurat_object@meta.data[[group.by]])
  if (is.null(conditions)) {
    conditions <- available_conditions
  } else {
    # Validate specified conditions exist
    missing <- setdiff(conditions, available_conditions)
    if (length(missing) > 0) {
      stop("Conditions not found in data: ", paste(missing, collapse = ", "))
    }
  }

  if (verbose) {
    cat("Conditions:", paste(conditions, collapse = ", "), "\n")
    cat("Number of conditions:", length(conditions), "\n\n")
  }

  # Create CellChat objects for each condition
  if (verbose) {
    cat("Creating CellChat objects...\n")
  }

  cellchat_list <- list()

  for (cond in conditions) {
    if (verbose) cat("  Processing condition:", cond, "\n")

    # Subset Seurat object for this condition
    cells_idx <- seurat_object@meta.data[[group.by]] == cond
    seurat_subset <- seurat_object[, cells_idx]

    # Extract data using compatibility function
    data_input <- extractSeuratData(seurat_subset,
                                     assay = "RNA",
                                     layer_or_slot = "data")
    meta_data <- seurat_subset@meta.data

    # Ensure cell type column is a factor
    if (!is.factor(meta_data[[cell.type.column]])) {
      meta_data[[cell.type.column]] <- factor(meta_data[[cell.type.column]])
    }

    # Create CellChat object
    cellchat <- CellChat::createCellChat(
      object = data_input,
      meta = meta_data,
      group.by = cell.type.column
    )

    # Set ligand-receptor database
    if (species == "human") {
      CellChatDB <- CellChat::CellChatDB.human
    } else {
      CellChatDB <- CellChat::CellChatDB.mouse
    }

    cellchat@DB <- CellChatDB

    # Preprocessing
    cellchat <- CellChat::subsetData(cellchat)
    cellchat <- CellChat::identifyOverExpressedGenes(cellchat, thresh.p = thresh)
    cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

    # Inference
    cellchat <- CellChat::computeCommunProb(
      cellchat,
      type = type,
      trim = trim,
      raw.use = raw.use,
      population.size = population.size
    )
    cellchat <- CellChat::filterCommunication(cellchat, min.cells = min.cells)

    # Compute pathway level
    cellchat <- CellChat::computeCommunProbPathway(cellchat)

    # Aggregate network
    cellchat <- CellChat::aggregateNet(cellchat)

    # Compute network centrality (needed for some differential analyses)
    if (compute.centrality) {
      cellchat <- CellChat::netAnalysis_computeCentrality(
        cellchat,
        slot.name = "netP",
        thresh = thresh.centrality
      )
    }

    # Store
    cellchat_list[[cond]] <- cellchat

    if (verbose) {
      cat("    Cells:", ncol(seurat_subset), "\n")
      cat("    Cell types:", length(unique(meta_data[[cell.type.column]])), "\n")
    }
  }

  # Print summary
  if (verbose && !run.analysis) {
    cat("\n=======================================================\n")
    cat("CellChat Objects Created\n")
    cat("=======================================================\n\n")
    cat("Created", length(cellchat_list), "CellChat objects\n")
    cat("Conditions:", paste(names(cellchat_list), collapse = ", "), "\n\n")

    cat("Next steps:\n")
    cat("  Use CellDiff functions for differential analysis:\n")
    cat("  - compareCell() for comprehensive analysis\n")
    cat("  - rankDiffM() for pathway ranking\n")
    cat("  - heatDiffM() for heatmaps\n")
    cat("  - scatterDiff2DM() for scatter plots\n\n")
  }

  # Optionally run differential analysis
  if (run.analysis) {
    if (verbose) {
      cat("\n=======================================================\n")
      cat("Running Differential Analysis\n")
      cat("=======================================================\n\n")
    }

    # Run compareCell
    analysis_results <- compareCell(
      object.list = cellchat_list,
      comparison = NULL,
      reference = reference,
      comparison_method = comparison_method,
      custom_comparisons = custom_comparisons,
      cell.type.strategy = cell.type.strategy,
      show.all = show.all,
      pThresh = pThresh,
      use_log2fc = use_log2fc,
      show_plots = show_plots,
      top.n = top.n,
      verbose = verbose
    )

    # Add CellChat objects to results
    analysis_results$cellchat_objects <- cellchat_list

    if (verbose) {
      cat("\n=======================================================\n")
      cat("Analysis Complete\n")
      cat("=======================================================\n\n")
      cat("Results include:\n")
      cat("  $cellchat_objects - CellChat objects\n")
      cat("  $pathway_analysis - Pathway ranking results\n")
      if (!is.null(analysis_results$heatmap)) cat("  $heatmap - Heatmap\n")
      if (!is.null(analysis_results$scatter)) cat("  $scatter - Scatter plot\n")
      cat("  $summary - Summary statistics\n\n")
    }

    return(analysis_results)
  } else {
    return(cellchat_list)
  }
}
