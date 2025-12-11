#' Comprehensive Differential Analysis of CellChat Objects
#'
#' A wrapper function that performs a complete differential cell-cell communication
#' analysis workflow, including pathway ranking, heatmaps, and visualizations.
#'
#' @param object.list A list of CellChat objects to compare
#' @param comparison Vector of condition indices or names to compare
#' @param reference Reference condition (index or name). Required for "all_vs_ref" method
#' @param comparison_method Comparison method: "all_vs_ref" (default), "all_vs_all", or "custom_pairs"
#' @param custom_comparisons List of custom comparison pairs (for custom_pairs method)
#' @param cell.type.strategy Strategy for handling missing cell types: "shared" (default) or "union"
#' @param show.all Logical. If TRUE, shows all pathways with visual indicators in comparison plots
#' @param pThresh P-value threshold for significance (default: 0.05)
#' @param use_log2fc Logical. Use log2 fold change instead of regular fold change
#' @param show_plots Character vector of plots to generate. Options: "barplot", "heatmap", "slope", "scatter", "network"
#' @param top.n Number of top pathways to highlight (default: NULL for all)
#' @param verbose Logical. If TRUE, prints progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item{pathway_analysis}{Results from rankDiffM including plots and significant pathways}
#'   \item{heatmap}{Results from heatDiffM if requested}
#'   \item{scatter}{Results from scatterDiff2DM if requested}
#'   \item{summary}{Summary statistics of the analysis}
#' }
#'
#' @details
#' This wrapper function streamlines the differential analysis workflow by:
#' \itemize{
#'   \item Automatically handling cell type alignment
#'   \item Running multiple analysis functions in sequence
#'   \item Organizing results in a structured format
#'   \item Providing informative progress messages
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage - compare all conditions to reference
#' results <- compareCell(
#'   object.list = cellchatlist,
#'   reference = "WT"
#' )
#'
#' # Comprehensive analysis with all plots
#' results <- compareCell(
#'   object.list = cellchatlist,
#'   reference = "WT",
#'   cell.type.strategy = "union",
#'   show.all = TRUE,
#'   show_plots = c("barplot", "heatmap", "slope", "scatter")
#' )
#'
#' # Custom comparisons
#' results <- compareCell(
#'   object.list = cellchatlist,
#'   comparison_method = "custom_pairs",
#'   custom_comparisons = list(
#'     c("WT", "KO"),
#'     c("KO", "DKO")
#'   ),
#'   show_plots = c("barplot", "heatmap")
#' )
#'
#' # Access results
#' results$pathway_analysis$comparison_barplot
#' results$pathway_analysis$all_significant_paths_full
#' results$summary
#' }
#'
#' @export
compareCell <- function(object.list,
                        comparison = NULL,
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

  cell.type.strategy <- match.arg(cell.type.strategy)

  # Initialize result list
  results <- list()

  # Validate inputs
  if (!is.list(object.list)) {
    stop("object.list must be a list of CellChat objects")
  }

  # Set up comparison indices
  if (is.null(comparison)) {
    comparison <- seq_along(object.list)
    if (verbose) {
      cat("Using all", length(object.list), "conditions in the analysis\n")
    }
  }

  # Validate reference for all_vs_ref
  if (comparison_method == "all_vs_ref" && is.null(reference)) {
    stop("Reference condition must be specified for 'all_vs_ref' comparison method")
  }

  # Print analysis summary
  if (verbose) {
    cat("\n=======================================================\n")
    cat("CellDiff Comprehensive Analysis\n")
    cat("=======================================================\n\n")
    cat("Settings:\n")
    cat("  Comparison method:", comparison_method, "\n")
    if (comparison_method == "all_vs_ref") {
      ref_name <- if (is.character(reference)) reference else names(object.list)[reference]
      cat("  Reference condition:", ref_name, "\n")
    }
    cat("  Cell type strategy:", cell.type.strategy, "\n")
    cat("  Show all pathways:", show.all, "\n")
    cat("  P-value threshold:", pThresh, "\n")
    cat("  Use log2 fold change:", use_log2fc, "\n")
    cat("  Plots to generate:", paste(show_plots, collapse = ", "), "\n\n")
  }

  # 1. Pathway Ranking Analysis
  if (verbose) cat("Step 1/", length(show_plots) + 1, ": Running pathway ranking analysis...\n")

  pathway_results <- rankDiffM(
    object.list = object.list,
    comparison = comparison,
    reference = reference,
    comparison_method = comparison_method,
    custom_comparisons = custom_comparisons,
    cell.type.strategy = cell.type.strategy,
    use_log2fc = use_log2fc,
    show_comparison_barplot = "barplot" %in% show_plots,
    show_comparison_heatmap = "heatmap" %in% show_plots,
    show_comparison_slope = "slope" %in% show_plots,
    show.all = show.all,
    pThresh = pThresh,
    top.n = top.n,
    return_top_paths = TRUE
  )

  results$pathway_analysis <- pathway_results

  step_counter <- 2

  # 2. Detailed Heatmap (if requested and not already in pathway analysis)
  if ("heatmap" %in% show_plots && comparison_method != "custom_pairs") {
    if (verbose) cat("Step", step_counter, "/", length(show_plots) + 1, ": Generating detailed heatmap...\n")

    heatmap_result <- tryCatch({
      heatDiffM(
        object.list = object.list,
        comparison = comparison,
        reference = reference,
        comparison_method = comparison_method,
        cell.type.strategy = cell.type.strategy,
        measure = "both",
        use_log2fc = use_log2fc,
        big_heatmap = TRUE,
        show_values = FALSE
      )
    }, error = function(e) {
      if (verbose) cat("  Warning: Heatmap generation failed:", e$message, "\n")
      NULL
    })

    results$heatmap <- heatmap_result
    step_counter <- step_counter + 1
  }

  # 3. Scatter Plot Analysis (if requested)
  if ("scatter" %in% show_plots && comparison_method == "all_vs_ref") {
    if (verbose) cat("Step", step_counter, "/", length(show_plots) + 1, ": Generating scatter plot...\n")

    scatter_result <- tryCatch({
      scatterDiff2DM(
        object.list = object.list,
        comparison = comparison,
        reference = reference,
        comparison_method = comparison_method,
        cell.type.strategy = cell.type.strategy,
        show.group.legend = TRUE,
        convex.hull = TRUE
      )
    }, error = function(e) {
      if (verbose) cat("  Warning: Scatter plot generation failed:", e$message, "\n")
      NULL
    })

    results$scatter <- scatter_result
    step_counter <- step_counter + 1
  }

  # 4. Network visualization (if requested)
  if ("network" %in% show_plots && !is.null(pathway_results$top_paths)) {
    if (verbose) cat("Step", step_counter, "/", length(show_plots) + 1, ": Generating network visualization...\n")

    # Use top pathways for network
    top_pathways <- head(pathway_results$top_paths, min(3, length(pathway_results$top_paths)))

    network_result <- tryCatch({
      networkLRDiff(
        object.list = object.list,
        comparison = comparison[1:min(2, length(comparison))],  # Network works best with 2 conditions
        pathways = top_pathways,
        node.size.factor = 1.2,
        edge.width.factor = 1.5,
        node.label.repel = TRUE
      )
    }, error = function(e) {
      if (verbose) cat("  Warning: Network visualization failed:", e$message, "\n")
      NULL
    })

    results$network <- network_result
  }

  # Generate summary statistics
  if (verbose) cat("\nGenerating summary statistics...\n")

  summary_stats <- list()

  # Count significant pathways per comparison
  if (!is.null(pathway_results$all_significant_paths_full)) {
    sig_counts <- sapply(pathway_results$all_significant_paths_full, length)
    summary_stats$significant_pathways_per_comparison <- sig_counts

    # Total unique significant pathways
    all_sig <- unique(unlist(pathway_results$all_significant_paths_full))
    summary_stats$total_unique_significant_pathways <- length(all_sig)

    # Common significant pathways (in all comparisons)
    if (length(pathway_results$all_significant_paths_full) > 1) {
      common_paths <- Reduce(intersect, pathway_results$all_significant_paths_full)
      summary_stats$common_significant_pathways <- common_paths
      summary_stats$n_common_pathways <- length(common_paths)
    }
  }

  # Cell type alignment info
  summary_stats$cell_type_strategy <- cell.type.strategy

  # Analysis parameters
  summary_stats$parameters <- list(
    comparison_method = comparison_method,
    pThresh = pThresh,
    use_log2fc = use_log2fc,
    show.all = show.all,
    top.n = top.n
  )

  results$summary <- summary_stats

  # Print summary
  if (verbose) {
    cat("\n=======================================================\n")
    cat("Analysis Summary\n")
    cat("=======================================================\n\n")

    if (!is.null(summary_stats$significant_pathways_per_comparison)) {
      cat("Significant pathways per comparison:\n")
      for (i in seq_along(summary_stats$significant_pathways_per_comparison)) {
        cat(sprintf("  %s: %d pathways\n",
                   names(summary_stats$significant_pathways_per_comparison)[i],
                   summary_stats$significant_pathways_per_comparison[i]))
      }
      cat("\nTotal unique significant pathways:", summary_stats$total_unique_significant_pathways, "\n")

      if (!is.null(summary_stats$n_common_pathways)) {
        cat("Common pathways (in all comparisons):", summary_stats$n_common_pathways, "\n")
      }
    }

    cat("\nResults stored in:\n")
    cat("  $pathway_analysis - Pathway ranking results\n")
    if (!is.null(results$heatmap)) cat("  $heatmap - Detailed heatmap\n")
    if (!is.null(results$scatter)) cat("  $scatter - Scatter plot analysis\n")
    if (!is.null(results$network)) cat("  $network - Network visualization\n")
    cat("  $summary - Summary statistics\n\n")
  }

  class(results) <- c("CellDiffAnalysis", "list")
  return(results)
}

#' Print method for CellDiffAnalysis
#'
#' @param x A CellDiffAnalysis object
#' @param ... Additional arguments (not used)
#' @export
print.CellDiffAnalysis <- function(x, ...) {
  cat("CellDiff Comprehensive Analysis Results\n")
  cat("========================================\n\n")

  if (!is.null(x$summary$significant_pathways_per_comparison)) {
    cat("Significant pathways per comparison:\n")
    print(x$summary$significant_pathways_per_comparison)
    cat("\nTotal unique:", x$summary$total_unique_significant_pathways, "\n")
  }

  cat("\nAvailable components:\n")
  cat("  $pathway_analysis - Pathway ranking results\n")
  if (!is.null(x$heatmap)) cat("  $heatmap - Detailed heatmap\n")
  if (!is.null(x$scatter)) cat("  $scatter - Scatter plot analysis\n")
  if (!is.null(x$network)) cat("  $network - Network visualization\n")
  cat("  $summary - Summary statistics\n\n")

  cat("Access plots:\n")
  cat("  results$pathway_analysis$comparison_barplot\n")
  cat("  results$pathway_analysis$comparison_heatmap\n")
  if (!is.null(x$heatmap)) cat("  results$heatmap\n")
  if (!is.null(x$scatter)) cat("  results$scatter\n\n")

  invisible(x)
}
