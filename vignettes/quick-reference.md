# CellDiff Quick Reference Guide

## Three Ways to Use CellDiff

### 1. Automatic (Fastest)
```r
results <- runCellChat(seurat_obj, group.by = "condition",
                       run.analysis = TRUE, reference = "WT",
                       show_plots = c("barplot", "heatmap"))
```

### 2. Manual Control
```r
cellchat_list <- runCellChat(seurat_obj, group.by = "condition")
rankDiffM(cellchat_list, reference = "WT")
heatDiffM(cellchat_list, reference = "WT")
```

### 3. From Existing CellChat Objects
```r
data(cellchatlist)
results <- compareCell(cellchatlist, reference = "WT")
```

## Essential Parameters

| Parameter | Values | Purpose |
|-----------|--------|---------|
| `run.analysis` | TRUE/FALSE | Automatic analysis? |
| `reference` | "WT" or 1 | Reference condition |
| `show_plots` | c("barplot", "heatmap", ...) | Which analyses |
| `comparison_method` | "all_vs_ref", "all_vs_all", "custom_pairs" | How to compare |
| `cell.type.strategy` | "shared", "union" | Missing cell types |
| `show.all` | TRUE/FALSE | Show all pathways |
| `species` | "human", "mouse" | CellChatDB |

## Available Analyses in show_plots

| Value | Function | What It Shows |
|-------|----------|---------------|
| `"barplot"` | rankDiffM | Pathway ranking |
| `"heatmap"` | heatDiffM | Differential heatmap |
| `"slope"` | rankDiffM | Pathway trends |
| `"scatter"` | scatterDiff2DM | Sender-receiver changes |
| `"network"` | networkLRDiff | Network visualization |

## Function Categories

### Wrapper Functions (Beginners Start Here)
- `runCellChat()` - Seurat → CellChat objects (+ optional analysis)
- `compareCell()` - Comprehensive differential analysis

### Multi-Condition Functions (3+ Groups)
- `rankDiffM()` - Pathway ranking
- `heatDiffM()` - Differential heatmaps
- `scatterDiff2DM()` - Sender-receiver scatter
- `ContriDiffM()` - L-R contributions

### Two-Condition Functions (Pairwise)
- `rankDiff()` - Pathway ranking (2 conditions)
- `heatDiff()` - Heatmap (2 conditions)
- `scatterDiff()` - Scatter (2 conditions)

### Utility Functions
- `findCommonUniquePatterns()` - Pattern discovery
- `compareCellChatsM()` - Global comparison
- `alignCellTypes()` - Cell type alignment

## Common Workflows

### Quick Analysis
```r
library(CellDiff)
data(pbmc_tutorial)

results <- runCellChat(pbmc_tutorial, group.by = "condition",
                       species = "human", run.analysis = TRUE,
                       reference = "WT")
results$pathway_analysis$comparison_barplot
```

### Custom Analysis
```r
# Step 1: Create CellChat objects
cc_list <- runCellChat(seurat_obj, group.by = "condition")

# Step 2: Explore pathways
rankDiffM(cc_list, reference = "Control", show.all = TRUE)

# Step 3: Deep dive specific pathways
heatDiffM(cc_list, reference = "Control",
          signaling = c("WNT", "TGFb"))

# Step 4: Analyze L-R pairs
ContriDiffM(cc_list, signaling = "WNT", reference = "Control")
```

### Pairwise Comparisons
```r
# All vs reference
results <- runCellChat(seurat_obj, group.by = "condition",
                       run.analysis = TRUE, reference = "WT",
                       comparison_method = "all_vs_ref")

# All vs all
results <- runCellChat(seurat_obj, group.by = "condition",
                       run.analysis = TRUE,
                       comparison_method = "all_vs_all")

# Custom pairs
custom <- list(c("WT", "KO"), c("WT", "DKO"))
results <- runCellChat(seurat_obj, group.by = "condition",
                       run.analysis = TRUE,
                       comparison_method = "custom_pairs",
                       custom_comparisons = custom)
```

## Accessing Results

```r
# From runCellChat with run.analysis = TRUE
results$cellchat_objects[["WT"]]        # CellChat objects
results$pathway_analysis$comparison_barplot  # Barplot
results$pathway_analysis$top_paths      # Top pathways
results$heatmap                         # Heatmap (if requested)
results$scatter                         # Scatter (if requested)
results$summary                         # Summary stats

# From individual functions
pathway_results <- rankDiffM(cc_list, reference = "WT")
pathway_results$comparison_barplot
pathway_results$top_paths
pathway_results$data
```

## Troubleshooting

### Error: "reference must be specified"
```r
# Add reference parameter
runCellChat(..., run.analysis = TRUE, reference = "WT")
```

### Error: "Column not found"
```r
# Check column names
colnames(seurat_obj@meta.data)

# Use exact names
runCellChat(seurat_obj, group.by = "condition",
            cell.type.column = "cell_type")
```

### No significant pathways
```r
# Show all pathways with indicators
rankDiffM(cc_list, reference = "WT", show.all = TRUE)
```

### Missing cell types across conditions
```r
# Use union strategy (includes all cell types)
rankDiffM(cc_list, reference = "WT", cell.type.strategy = "union")

# Or use shared strategy (only common types)
rankDiffM(cc_list, reference = "WT", cell.type.strategy = "shared")
```

## Tips

✅ **Start simple**: Use `run.analysis = TRUE` for quick results
✅ **Use tutorial data**: `data(pbmc_tutorial)` for learning
✅ **Check vignettes**: `vignette("getting-started", package = "CellDiff")`
✅ **Named access**: Use comparison names instead of indices
✅ **Debug mode**: `options(CellDiff.debug = TRUE)` for troubleshooting
✅ **Show all**: Use `show.all = TRUE` to see all pathways
✅ **Save results**: `saveRDS(results, "analysis.rds")`

## Getting Help

```r
# Function help
?runCellChat
?compareCell
?rankDiffM

# Vignettes
vignette("getting-started", package = "CellDiff")
vignette("CellDiff", package = "CellDiff")
browseVignettes("CellDiff")

# Tutorial data
?pbmc_tutorial
data(pbmc_tutorial)
```

## Citation

If you use CellDiff:
- CellDiff: Guo, K. (2024)
- CellChat: Jin, S., et al. (2021). Nature Communications, 12, 1088.
