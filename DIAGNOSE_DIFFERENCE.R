#!/usr/bin/env Rscript

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Diagnosing Differences Between CellChat Objects\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

library(CellDiff)
library(CellChat)

# Load the datasets
data(celllist)

cat("Please load your 'cells' object (manual pipeline) into the environment.\n")
cat("Then compare the following:\n\n")

cat("1. Check pathway communication strengths:\n")
cat("───────────────────────────────────────────────────────────────\n")
cat("For celllist (runCellChat):\n")
cat("  pathways_celllist_wt <- celllist$WT@netP$prob\n")
cat("  # Check specific pathway\n")
cat("  sum(pathways_celllist_wt[,,'PTN'])  # Sum of PTN pathway\n\n")

cat("For cells (manual):\n")
cat("  pathways_cells_wt <- cells$WT@netP$prob\n")
cat("  sum(pathways_cells_wt[,,'PTN'])  # Sum of PTN pathway\n\n")

cat("2. Check number of identified interactions:\n")
cat("───────────────────────────────────────────────────────────────\n")
cat("For celllist:\n")
cat("  nrow(celllist$WT@LR$LRsig)  # Significant L-R pairs\n\n")

cat("For cells:\n")
cat("  nrow(cells$WT@LR$LRsig)  # Significant L-R pairs\n\n")

cat("3. Check data.project slot:\n")
cat("───────────────────────────────────────────────────────────────\n")
cat("For celllist:\n")
cat("  dim(celllist$WT@data.project)  # Should be NULL if not projected\n\n")

cat("For cells:\n")
cat("  dim(cells$WT@data.project)  # Check if PPI projection was done\n\n")

cat("4. Check specific parameters used:\n")
cat("───────────────────────────────────────────────────────────────\n")
cat("Parameters are stored in @options slot:\n")
cat("  celllist$WT@options\n")
cat("  cells$WT@options\n\n")

# If objects exist, do comparison
if (exists("celllist") && exists("cells")) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  Automatic Comparison (both objects loaded)\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  cat("Comparing WT condition:\n\n")

  # Compare pathways
  cat("Number of pathways:\n")
  cat(sprintf("  celllist: %d\n", length(celllist$WT@netP$pathways)))
  cat(sprintf("  cells: %d\n", length(cells$WT@netP$pathways)))

  # Compare specific pathway values
  if ("PTN" %in% celllist$WT@netP$pathways && "PTN" %in% cells$WT@netP$pathways) {
    ptn_celllist <- sum(celllist$WT@netP$prob[,,'PTN'])
    ptn_cells <- sum(cells$WT@netP$prob[,,'PTN'])
    cat("\nPTN pathway total strength:\n")
    cat(sprintf("  celllist: %.6e\n", ptn_celllist))
    cat(sprintf("  cells: %.6e\n", ptn_cells))
    cat(sprintf("  Ratio (cells/celllist): %.2f\n", ptn_cells/ptn_celllist))
  }

  # Compare L-R pairs
  cat("\nSignificant L-R pairs:\n")
  cat(sprintf("  celllist: %d\n", nrow(celllist$WT@LR$LRsig)))
  cat(sprintf("  cells: %d\n", nrow(cells$WT@LR$LRsig)))

  # Check options
  cat("\nParameters used (celllist):\n")
  if (!is.null(celllist$WT@options)) {
    print(celllist$WT@options)
  } else {
    cat("  No @options slot found\n")
  }

  cat("\nParameters used (cells):\n")
  if (!is.null(cells$WT@options)) {
    print(cells$WT@options)
  } else {
    cat("  No @options slot found\n")
  }

} else {
  cat("\nNote: Load both 'celllist' and 'cells' to enable automatic comparison\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  To Fix the Issue\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Once you identify the difference, update runCellChat to match your\n")
cat("manual pipeline exactly. Key parameters to check:\n\n")

cat("1. identifyOverExpressedGenes:\n")
cat("   - thresh.p (default: 0.05)\n")
cat("   - thresh.pc (default: 0.1)\n")
cat("   - thresh.fc (default: 0.1)\n\n")

cat("2. computeCommunProb:\n")
cat("   - type (default: 'triMean')\n")
cat("   - trim (default: 0.1)\n")
cat("   - population.size (default: TRUE)\n")
cat("   - raw.use (default: FALSE or TRUE - CHECK THIS!)\n\n")

cat("3. netAnalysis_computeCentrality:\n")
cat("   - thresh (default: 0.05)\n\n")

cat("The most likely culprit is 'raw.use' or 'type' in computeCommunProb!\n\n")
