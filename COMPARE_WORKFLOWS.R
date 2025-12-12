#!/usr/bin/env Rscript

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Comparing: runCellChat vs Manual CellChat Pipeline\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

library(CellDiff)
library(CellChat)

# Load both datasets
data(celllist)  # Created with runCellChat
# Assuming you have 'cells' in your environment (manual pipeline)

cat("Dataset 1: celllist (runCellChat)\n")
cat("───────────────────────────────────────────────────────────────\n")

if (exists("celllist")) {
  cat("Number of conditions:", length(celllist), "\n")
  cat("Condition names:", names(celllist), "\n\n")

  for (cond in names(celllist)) {
    obj <- celllist[[cond]]
    cat(sprintf("Condition: %s\n", cond))
    cat(sprintf("  Class: %s\n", class(obj)[1]))
    cat(sprintf("  Cells: %d\n", ncol(obj@data.signaling)))
    cat(sprintf("  Cell groups: %d\n", length(levels(obj@idents))))
    cat(sprintf("  Cell group names: %s\n", paste(levels(obj@idents), collapse=", ")))

    # Check database
    if (!is.null(obj@DB)) {
      cat(sprintf("  Database interactions: %d\n", nrow(obj@DB$interaction)))
    }

    # Check network
    if (!is.null(obj@net$count)) {
      cat(sprintf("  Network interactions (count): %d non-zero\n", sum(obj@net$count > 0)))
    }
    if (!is.null(obj@netP$pathways)) {
      cat(sprintf("  Pathways identified: %d\n", length(obj@netP$pathways)))
      cat(sprintf("  Top 5 pathways: %s\n", paste(head(obj@netP$pathways, 5), collapse=", ")))
    }
    cat("\n")
  }
}

cat("\n")
cat("Dataset 2: cells (Manual pipeline)\n")
cat("───────────────────────────────────────────────────────────────\n")

# If you have the manual pipeline results, examine them
if (exists("cells")) {
  cat("Number of conditions:", length(cells), "\n")
  cat("Condition names:", names(cells), "\n\n")

  for (cond in names(cells)) {
    obj <- cells[[cond]]
    cat(sprintf("Condition: %s\n", cond))
    cat(sprintf("  Class: %s\n", class(obj)[1]))
    cat(sprintf("  Cells: %d\n", ncol(obj@data.signaling)))
    cat(sprintf("  Cell groups: %d\n", length(levels(obj@idents))))
    cat(sprintf("  Cell group names: %s\n", paste(levels(obj@idents), collapse=", ")))

    # Check database
    if (!is.null(obj@DB)) {
      cat(sprintf("  Database interactions: %d\n", nrow(obj@DB$interaction)))
    }

    # Check network
    if (!is.null(obj@net$count)) {
      cat(sprintf("  Network interactions (count): %d non-zero\n", sum(obj@net$count > 0)))
    }
    if (!is.null(obj@netP$pathways)) {
      cat(sprintf("  Pathways identified: %d\n", length(obj@netP$pathways)))
      cat(sprintf("  Top 5 pathways: %s\n", paste(head(obj@netP$pathways, 5), collapse=", ")))
    }
    cat("\n")
  }
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  Key Differences to Check\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Please verify the following:\n\n")

cat("1. Data Input:\n")
cat("   - runCellChat uses: GetAssayData(seurat, assay='RNA', slot='data')\n")
cat("   - Your pipeline uses: GetAssayData(wt, assay='RNA', slot='data')\n")
cat("   ✓ Should be the same\n\n")

cat("2. Meta Data:\n")
cat("   - runCellChat: Creates meta from metadata column\n")
cat("   - Your pipeline: Uses Idents(wt) for labels\n")
cat("   ? Could be different if Idents != metadata column\n\n")

cat("3. Database:\n")
cat("   - runCellChat: CellChatDB.mouse or CellChatDB.human\n")
cat("   - Your pipeline: sddB.use <- CellChatDB.mouse\n")
cat("   ? Check if database was modified or subset\n\n")

cat("4. Subsetting:\n")
cat("   - Both use: subsetData()\n")
cat("   ✓ Should be the same\n\n")

cat("5. OverExpressed Genes:\n")
cat("   - runCellChat: identifyOverExpressedGenes(thresh.p=0.05, thresh.pc=0.1, thresh.fc=0.1)\n")
cat("   - Your pipeline: identifyOverExpressedGenes() with defaults\n")
cat("   ? Check default parameter values\n\n")

cat("6. Communication Probability:\n")
cat("   - runCellChat: computeCommunProb(type='triMean', trim=0.1, population.size=TRUE)\n")
cat("   - Your pipeline: computeCommunProb() with defaults\n")
cat("   ? Check parameter differences\n\n")

cat("7. Filtering:\n")
cat("   - runCellChat: filterCommunication(min.cells=10)\n")
cat("   - Your pipeline: filterCommunication(min.cells=10)\n")
cat("   ✓ Should be the same\n\n")

cat("8. Missing Step:\n")
cat("   - Your pipeline: NO projectData(wt, PPI.mouse) - commented out\n")
cat("   - runCellChat: Does NOT include projectData\n")
cat("   ✓ Both skip this step\n\n")

cat("═══════════════════════════════════════════════════════════════\n\n")

cat("To debug further, please run:\n")
cat("  1. Compare number of pathways identified in each\n")
cat("  2. Compare specific pathway values (e.g., PTN, GAS, etc.)\n")
cat("  3. Check if parameters differ in key steps\n\n")

cat("Run this with both 'celllist' and 'cells' loaded in your environment\n")
cat("to see the actual differences.\n\n")
