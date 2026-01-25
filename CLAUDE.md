# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CellDiff is an R package for differential analysis of cell-cell communication between experimental groups. It extends CellChat functionality to compare communication patterns across multiple biological conditions (e.g., normal vs. disease, treated vs. untreated).

## Build and Development Commands

```r
# Generate documentation from roxygen2 comments
devtools::document()

# Run package checks
devtools::check()

# Run tests
devtools::test()

# Install locally
devtools::install()

# Build vignettes
devtools::build_vignettes()
```

## Architecture

### Layer Structure

1. **Wrapper Functions (Entry Points)**
   - `runCellChat()` - Creates CellChat objects from Seurat objects (39 customizable parameters)
   - `compareCell()` - One-step comprehensive differential analysis

2. **Multi-Condition Functions (3+ groups)** - Files ending in `M.R`
   - `rankDiffM()` - Pathway ranking with comparison methods (all_vs_ref, all_vs_all, custom_pairs)
   - `heatDiffM()` - Differential heatmaps with ComplexHeatmap
   - `scatterDiff2DM()` - Sender-receiver scatter plots
   - `ContriDiffM()` - L-R pair contribution analysis

3. **Two-Condition Functions (Pairwise)**
   - `rankDiff()`, `heatDiff()`, `scatterDiff()`, `scatterDiff2D()`, `ContriDiff()`, `netDiff()`, `networkLRDiff()`

4. **Utilities**
   - `R/utils.R` - Core utility functions (alignCellTypes, assignColors, calculateScores)
   - `R/misc.R` - Miscellaneous helpers
   - `R/seurat_compat.R` - Seurat v4/v5 compatibility

### Key Design Patterns

- Multi-condition functions (suffix `M`) handle 3+ groups; base functions handle 2 conditions
- `cell.type.strategy` parameter: "shared" (conservative) vs "union" (inclusive) for missing cell types
- `comparison_method`: "all_vs_ref", "all_vs_all", or "custom_pairs"
- Named list outputs for easier result access (e.g., `results$all_significant_paths_full[["KO vs WT"]]`)

### Dependencies

Core: CellChat, Seurat, ggplot2, ComplexHeatmap, circlize, igraph, reshape2, dplyr, tidyr

## Data Files

- `data/celllist.rda` - 2-condition example (NL, LS)
- `data/cellchatlist.rda` - 3-condition example (WT, KO, DKO)
- `data/pbmc_tutorial.rda` - Tutorial Seurat object (1,500 cells, 3 conditions, 5 cell types)

## Testing Workflows

```r
# Load tutorial data
data(pbmc_tutorial)

# Quick test with automatic analysis
results <- runCellChat(pbmc_tutorial, group.by = "condition",
                       species = "human", run.analysis = TRUE, reference = "WT")

# Test with existing CellChat objects
data(cellchatlist)
rankDiffM(cellchatlist, reference = "WT")
```

## Code Conventions

- roxygen2 for documentation (RoxygenNote: 7.3.3)
- Functions return plots directly by default; use `return.data = TRUE` for underlying data
- ComplexHeatmap for all heatmaps (pheatmap deprecated)
- Base R preferred over tidyverse for core logic
