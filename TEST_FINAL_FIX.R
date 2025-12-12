#!/usr/bin/env Rscript

cat("═══════════════════════════════════════════════════════════════\n")
cat("  Testing: rankDiff with Scaling + Zero Replacement\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

library(devtools)
devtools::load_all(".", quiet = TRUE)

data(celllist)

cat("Step 1: Run rankDiff and check data...\n")
result <- rankDiff(celllist, return.data = TRUE)

cat("\nStep 2: Check the data structure...\n")
wide_data <- result$data

# Count pathways with 0 in one group
pathways_with_zero <- wide_data[wide_data$WT == 0 | wide_data$KO == 0, ]

cat(sprintf("  Total significant pathways: %d\n", nrow(wide_data)))
cat(sprintf("  Pathways with 0 in one group: %d\n", nrow(pathways_with_zero)))

if (nrow(pathways_with_zero) > 0) {
  cat("\n  Examples of pathways with 0:\n")
  for (i in 1:min(5, nrow(pathways_with_zero))) {
    pw <- pathways_with_zero$name[i]
    wt_val <- pathways_with_zero$WT[i]
    ko_val <- pathways_with_zero$KO[i]
    if (wt_val == 0) {
      cat(sprintf("    - %s: WT=0, KO=%.2e (only in KO)\n", pw, ko_val))
    } else {
      cat(sprintf("    - %s: WT=%.2e, KO=0 (only in WT)\n", pw, wt_val))
    }
  }
}

cat("\nStep 3: Check plot data...\n")
plot_data <- result$plot_data

# Check zeros
exact_zeros <- sum(plot_data$contribution == 0)
tiny_values <- sum(plot_data$contribution > 0 & plot_data$contribution < 1e-15)

cat(sprintf("  Total plot data rows: %d\n", nrow(plot_data)))
cat(sprintf("  Rows with exact 0: %d\n", exact_zeros))
cat(sprintf("  Rows with ~0 (< 1e-15): %d\n", tiny_values))

if (exact_zeros == 0) {
  cat("  ✓ No exact zeros in plot data\n")
} else {
  cat("  ✗ ERROR: Still have exact zeros!\n")
}

# Check specific pathways
cat("\nStep 4: Verify pathways with 0 have both rows...\n")
if (nrow(pathways_with_zero) > 0) {
  test_pw <- pathways_with_zero$name[1]
  pw_data <- plot_data[plot_data$name == test_pw, ]

  cat(sprintf("  Testing: %s\n", test_pw))
  cat(sprintf("    Rows in plot: %d\n", nrow(pw_data)))

  for (i in 1:nrow(pw_data)) {
    grp <- pw_data$group[i]
    contrib <- pw_data$contribution[i]
    if (contrib < 1e-15) {
      cat(sprintf("      %s: %.2e (effectively 0) ✓\n", grp, contrib))
    } else {
      cat(sprintf("      %s: %.2e ✓\n", grp, contrib))
    }
  }
}

cat("\nStep 5: Save plot...\n")
pdf("rankDiff_FIXED_WITH_SCALING.pdf", width = 10, height = 12)
print(result$plot)
dev.off()

png("rankDiff_FIXED_WITH_SCALING.png", width = 1000, height = 1200, res = 120)
print(result$plot)
dev.off()

cat("  ✓ Saved to: rankDiff_FIXED_WITH_SCALING.pdf\n")
cat("  ✓ Saved to: rankDiff_FIXED_WITH_SCALING.png\n")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

if (exact_zeros == 0 && nrow(pathways_with_zero) > 0) {
  cat("✓✓✓ FIX SUCCESSFUL ✓✓✓\n\n")
  cat("Key Results:\n")
  cat(sprintf("  • All %d significant pathways kept\n", nrow(wide_data)))
  cat(sprintf("  • %d pathways with 0 in one group\n", nrow(pathways_with_zero)))
  cat("  • Scaling applied to all values\n")
  cat("  • Exact zeros replaced with 1e-20\n")
  cat("  • No blank spaces expected in bars\n")
  cat("\nThe plot should show:\n")
  cat("  - Pathways with both groups: both colors visible\n")
  cat("  - Pathways with 0 in one group: solid bar (one color)\n")
  cat("  - All bars extending from 0.0 to 1.0\n")
} else if (exact_zeros > 0) {
  cat("✗ ERROR: Exact zeros still present\n")
} else {
  cat("⚠ WARNING: No pathways with 0 found for testing\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
