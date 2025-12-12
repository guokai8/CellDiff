#!/usr/bin/env Rscript

library(CellChat)

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Checking CellChat Function Default Parameters\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("1. identifyOverExpressedGenes() defaults:\n")
cat("───────────────────────────────────────────────────────────────\n")
# Check function signature
if (exists("identifyOverExpressedGenes", where = "package:CellChat")) {
  args_overexp <- formals(CellChat::identifyOverExpressedGenes)
  cat("  thresh.p:", args_overexp$thresh.p, "\n")
  cat("  thresh.pc:", args_overexp$thresh.pc, "\n")
  cat("  thresh.fc:", args_overexp$thresh.fc, "\n")
}

cat("\nrunCellChat uses:\n")
cat("  thresh.p: 0.05 (parameter: thresh)\n")
cat("  thresh.pc: 0.1\n")
cat("  thresh.fc: 0.1\n")

cat("\nYour manual pipeline uses:\n")
cat("  identifyOverExpressedGenes(wt)  # defaults\n")

cat("\n2. computeCommunProb() defaults:\n")
cat("───────────────────────────────────────────────────────────────\n")
if (exists("computeCommunProb", where = "package:CellChat")) {
  args_prob <- formals(CellChat::computeCommunProb)
  cat("  type:", ifelse(is.null(args_prob$type), "NULL", args_prob$type), "\n")
  cat("  trim:", args_prob$trim, "\n")
  cat("  population.size:", args_prob$population.size, "\n")
  cat("  raw.use:", args_prob$raw.use, "\n")
}

cat("\nrunCellChat uses:\n")
cat("  type: 'triMean' (default in function signature)\n")
cat("  trim: 0.1\n")
cat("  population.size: TRUE\n")
cat("  raw.use: TRUE\n")

cat("\nYour manual pipeline uses:\n")
cat("  computeCommunProb(wt)  # defaults\n")

cat("\n3. filterCommunication() defaults:\n")
cat("───────────────────────────────────────────────────────────────\n")
if (exists("filterCommunication", where = "package:CellChat")) {
  args_filter <- formals(CellChat::filterCommunication)
  cat("  min.cells:", args_filter$min.cells, "\n")
}

cat("\nBoth use:\n")
cat("  min.cells: 10\n")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  KEY FINDING\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("The difference is likely in:\n")
cat("  1. identifyOverExpressedGenes parameters\n")
cat("  2. computeCommunProb parameters\n\n")

cat("RECOMMENDATION:\n")
cat("  Check if CellChat's default values differ from what runCellChat uses.\n")
cat("  If defaults are different, that explains the different results!\n\n")

cat("To fix: Modify runCellChat to match CellChat defaults exactly,\n")
cat("OR explicitly set parameters in your manual pipeline.\n\n")
