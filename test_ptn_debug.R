#!/usr/bin/env Rscript
# Debug script to trace PTN significance issue
# Run this with your actual data

# Enable debug mode
options(rankDiffM.debug = TRUE)

# Load the package
devtools::load_all(".")

# Replace 'your_cellchat_list' with your actual data
# cellchatlist <- your_cellchat_list

cat("=================================================================\n")
cat("DEBUG SCRIPT FOR PTN SIGNIFICANCE TRACKING\n")
cat("=================================================================\n\n")

cat("INSTRUCTIONS:\n")
cat("1. Load your CellChat data into 'cellchatlist'\n")
cat("2. Update the comparison and reference below to match your analysis\n")
cat("3. Run this script\n\n")

cat("Expected output from your description:\n")
cat("  - PTN should be in all_significant_paths_full[[1]] (5AC1-/- vs C1+/+)\n")
cat("  - PTN should NOT be in all_significant_paths_full[[2]] (5AC1+/+ vs C1+/+)\n")
cat("  - PTN should NOT be in all_significant_paths_full[[3]] (C-/- vs C1+/+)\n\n")

# Run the analysis
results <- rankDiffM(
  object.list = cellchatlist,  # YOUR DATA HERE
  comparison = c(1, 2, 3, 4),  # UPDATE THIS
  reference = 1,               # UPDATE THIS (position of C1+/+)
  comparison_method = "all_vs_ref",
  show_comparison_barplot = TRUE,
  show_comparison_heatmap = TRUE,
  show.all = TRUE,
  pThresh = 0.05
)

cat("\n\n=================================================================\n")
cat("CHECKING RESULTS\n")
cat("=================================================================\n\n")

cat("1. Checking all_significant_paths_full for PTN:\n")
for (i in 1:length(results$all_significant_paths_full)) {
  has_ptn <- "PTN" %in% results$all_significant_paths_full[[i]]
  cat(sprintf("   Comparison %d: PTN present = %s\n", i, has_ptn))
}

cat("\n2. Checking comparison_barplot data for PTN:\n")
if (!is.null(results$comparison_barplot)) {
  ptn_data <- results$comparison_barplot$data[
    results$comparison_barplot$data$pathway == "PTN",
  ]
  if (nrow(ptn_data) > 0) {
    cat("   PTN found in barplot:\n")
    print(ptn_data[, c("pathway", "comparison", "significant", "bar_alpha")])
  } else {
    cat("   PTN NOT found in barplot\n")
  }
} else {
  cat("   comparison_barplot is NULL\n")
}

cat("\n3. Checking star labels for PTN:\n")
if (!is.null(results$comparison_barplot)) {
  for (layer in results$comparison_barplot$layers) {
    if (inherits(layer$geom, "GeomText")) {
      ptn_stars <- layer$data[layer$data$pathway == "PTN", ]
      if (nrow(ptn_stars) > 0) {
        cat("   PTN stars found:\n")
        print(ptn_stars)
      } else {
        cat("   No PTN stars (correct if PTN not significant in any shown comparison)\n")
      }
    }
  }
}

cat("\n\nIf you see debug output above starting with 'STAR DEBUG [PTN in ...]',\n")
cat("that shows exactly which comparisons PTN was found in and whether stars were added.\n")
