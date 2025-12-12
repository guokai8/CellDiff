#!/usr/bin/env Rscript

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Identifying Parameter Differences\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

library(CellChat)

cat("Checking CellChat function defaults...\n\n")

cat("1. computeCommunProb defaults:\n")
cat("───────────────────────────────────────────────────────────────\n")

# Get the function
func <- getFromNamespace("computeCommunProb", "CellChat")
defaults <- formals(func)

cat("CellChat defaults:\n")
for (param in c("type", "trim", "LR.use", "raw.use", "population.size",
                "distance.use", "interaction.range", "scale.distance", "contact.dependent")) {
  if (param %in% names(defaults)) {
    val <- defaults[[param]]
    if (is.null(val)) {
      cat(sprintf("  %s: NULL\n", param))
    } else if (is.logical(val)) {
      cat(sprintf("  %s: %s\n", param, val))
    } else if (is.numeric(val)) {
      cat(sprintf("  %s: %g\n", param, val))
    } else if (is.character(val)) {
      cat(sprintf("  %s: '%s'\n", param, val))
    } else {
      cat(sprintf("  %s: %s\n", param, deparse(val)))
    }
  }
}

cat("\nrunCellChat uses (from R/runCellChat.R line 169-176):\n")
cat("  type: 'triMean' (matched from defaults)\n")
cat("  trim: 0.1\n")
cat("  raw.use: TRUE  <-- CHECK IF THIS DIFFERS!\n")
cat("  population.size: TRUE\n")
cat("  distance.use: TRUE\n")
cat("  interaction.range: 250\n")
cat("  scale.distance: 0.01\n")
cat("  contact.dependent: TRUE\n")

cat("\n2. identifyOverExpressedGenes defaults:\n")
cat("───────────────────────────────────────────────────────────────\n")

func2 <- getFromNamespace("identifyOverExpressedGenes", "CellChat")
defaults2 <- formals(func2)

cat("CellChat defaults:\n")
for (param in c("thresh.p", "thresh.pc", "thresh.fc")) {
  if (param %in% names(defaults2)) {
    val <- defaults2[[param]]
    if (is.null(val)) {
      cat(sprintf("  %s: NULL\n", param))
    } else {
      cat(sprintf("  %s: %g\n", param, val))
    }
  }
}

cat("\nrunCellChat uses (from R/runCellChat.R line 178-180):\n")
cat("  thresh: 0.05 (maps to thresh.p)\n")
cat("  thresh.pc: 0.1\n")
cat("  thresh.fc: 0.1\n")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  KEY FINDINGS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Compare the values above!\n\n")

cat("If any default differs from runCellChat's values, that's the cause.\n")
cat("Most likely culprit: raw.use parameter\n\n")

cat("To fix: Edit R/runCellChat.R and change the parameter to match\n")
cat("CellChat's default shown above.\n\n")

cat("═══════════════════════════════════════════════════════════════\n")
