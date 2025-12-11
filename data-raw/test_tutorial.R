# Test Script for CellDiff Tutorial Data
# Demonstrates the complete workflow with runCellChat and compareCell

library(CellDiff)
library(Seurat)

# Load tutorial data
load("data/pbmc_tutorial.rda")

cat("\n=======================================================\n")
cat("Testing CellDiff Workflow with Tutorial Data\n")
cat("=======================================================\n\n")

# Examine the data
cat("Tutorial Data Summary:\n")
cat("  Conditions:", paste(unique(pbmc_tutorial$condition), collapse = ", "), "\n")
cat("  Cell types:", paste(unique(pbmc_tutorial$cell_type), collapse = ", "), "\n")
cat("  Total cells:", ncol(pbmc_tutorial), "\n\n")

# Step 1: Create CellChat objects from Seurat
cat("Step 1: Creating CellChat objects...\n")
cat("---------------------------------------\n")

cellchat_list <- runCellChat(
  seurat_object = pbmc_tutorial,
  group.by = "condition",
  species = "human",
  cell.type.column = "cell_type",
  min.cells = 10,
  verbose = TRUE
)

cat("\nCellChat objects created:\n")
print(names(cellchat_list))

# Step 2: Run comprehensive differential analysis
cat("\n\nStep 2: Running differential analysis...\n")
cat("---------------------------------------\n")

results <- compareCell(
  object.list = cellchat_list,
  reference = "WT",
  cell.type.strategy = "union",
  show.all = TRUE,
  pThresh = 0.05,
  use_log2fc = TRUE,
  show_plots = c("barplot", "heatmap"),
  verbose = TRUE
)

# Examine results
cat("\n\nResults Summary:\n")
cat("---------------------------------------\n")
print(results)

cat("\n\nSignificant pathways per comparison:\n")
if (!is.null(results$summary$significant_pathways_per_comparison)) {
  print(results$summary$significant_pathways_per_comparison)
} else {
  cat("No significant pathways detected (this is expected with synthetic data)\n")
}

# Access individual components
cat("\n\nAvailable result components:\n")
cat("  - pathway_analysis: ", !is.null(results$pathway_analysis), "\n")
cat("  - heatmap: ", !is.null(results$heatmap), "\n")
cat("  - summary: ", !is.null(results$summary), "\n")

cat("\n\nTest completed successfully!\n")
cat("=======================================================\n\n")

cat("Next steps:\n")
cat("  # View comparison barplot\n")
cat("  results$pathway_analysis$comparison_barplot\n\n")
cat("  # View significant pathways\n")
cat("  results$pathway_analysis$all_significant_paths_full\n\n")
cat("  # Or use individual CellDiff functions:\n")
cat("  rankDiffM(cellchat_list, reference = 'WT')\n")
cat("  heatDiffM(cellchat_list, reference = 'WT')\n\n")
