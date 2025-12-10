# CellDiff

## Differential Analysis of Cell-Cell Communication

[![GitHub release](https://img.shields.io/github/release/guokai8/CellDiff.svg)](https://github.com/guokai8/CellDiff/releases)
[![R-CMD-check](https://github.com/guokai8/CellDiff/workflows/R-CMD-check/badge.svg)](https://github.com/guokai8/CellDiff/actions)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

CellDiff enhances visualization and statistical analysis of cell-cell communication differences between experimental groups.

## Overview

CellDiff is an R package that extends the functionality of [CellChat](https://github.com/sqjin/CellChat) by providing specialized tools for differential cell-cell communication analysis between multiple conditions. While CellChat enables the inference and analysis of cell-cell communication networks from single-cell RNA sequencing data, CellDiff focuses on identifying and visualizing the differences in these communication patterns across biological conditions (e.g., normal vs. disease, treated vs. untreated, multiple genotypes).

## Features

* **Multi-group comparison**: Compare cell-cell communication across multiple conditions simultaneously
* **Flexible comparison methods**: Support for all-vs-all, all-vs-reference, and custom pairwise comparisons
* **Cell type alignment**: Handle missing cell types with "shared" or "union" strategies
* **Differential pathway analysis**: Identify signaling pathways that are significantly altered between conditions
* **Advanced visualization options**:
  - Show all pathways with significance indicators (stars, fading)
  - Comparison-specific pathway display
  - Named list outputs for easier data access
* **Rich visualizations**: Venn diagrams, UpSet plots, heatmaps, network graphs, and timeline plots
* **Ligand-receptor analysis**: Compare L-R pair contributions with customizable plots and heatmaps
* **Pattern discovery**: Identify common and unique signaling patterns across conditions
* **Sender-receiver analysis**: Analyze how cells change their roles as signaling senders or receivers
* **Integrated visualizations**: ComplexHeatmap integration with consistent color schemes
* **No external dependencies on dplyr/tidyr**: Pure base R + specialized visualization packages
## Installation

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install CellDiff from GitHub
devtools::install_github("guokai8/CellDiff")
```

### Dependencies

CellDiff requires the following packages:
- CellChat (>= 1.1.0)
- ggplot2
- reshape2
- circlize
- ComplexHeatmap
- cowplot
- igraph
- ggrepel
- scales
- VennDetail (optional, for Venn diagrams)
- grid

## What's New in Version 0.1.4

### Major Improvements

1. **Enhanced Cell Type Alignment**
   - **NEW**: `cell.type.strategy` parameter for handling missing cell types across conditions
   - **"shared"** strategy: Uses only cell types present in all conditions (safer, more conservative)
   - **"union"** strategy: Uses all cell types from any condition, filling missing data with zeros
   - Automatically reports alignment statistics and warns about missing cell types
   - Available in: `rankDiffM`, `heatDiffM`, `scatterDiff2DM`, and `ContriDiffM`

2. **Advanced Comparison Visualization (rankDiffM)**
   - **NEW**: `show.all` parameter for comprehensive pathway visualization
   - **show.all = FALSE** (default): Shows only pathways significant in each specific comparison
   - **show.all = TRUE**: Shows all pathways with visual indicators:
     - Faded bars (alpha=0.3) for non-significant pathway-comparison combinations
     - Significance stars (*p<0.05, **p<0.01, ***p<0.001)
     - Informative subtitle explaining the visualization
   - **FIXED**: Comparison-specific pathway display - pathways now only appear in panels where they're actually significant
   - **FIXED**: Star assignment bug - stars now only appear for pathways passing both p-value AND fold-change thresholds

3. **Named List Output**
   - **NEW**: All output lists now have meaningful comparison names
   - Access results by comparison name instead of numeric index
   - Example: `results$all_significant_paths_full[["KO vs WT"]]` instead of `[[1]]`
   - Applies to: `all_significant_paths_full`, `all_significant_paths`, `data`, `plots`, `heatmaps`

4. **Custom Pairwise Comparisons**
   - `rankDiffM`, `heatDiffM`, `scatterDiff2DM`, and `ContriDiffM` now support `comparison_method = "custom_pairs"`
   - Specify exactly which condition pairs to compare: `custom_comparisons = list(c("WT", "KO"), c("KO", "DKO"))`
   - Generates focused visualizations for your specific comparisons of interest

5. **Enhanced findCommonUniquePatterns**
   - **NEW**: Venn diagram visualization using VennDetail (up to 7 conditions)
   - **NEW**: UpSet plot for complex multi-way intersections
   - **NEW**: Barplot showing total, unique, and common pattern counts
   - **NEW**: Network visualization of top signaling interactions
   - Identifies patterns that are shared or unique across conditions

6. **Improved compareCellChatsM**
   - **FIXED**: Replaced pheatmap with ComplexHeatmap for consistency
   - **FIXED**: Removed dplyr dependency - uses pure base R
   - **NEW**: Enhanced timeline plots with optional labeling
   - **NEW**: Pathway contribution timeline with labeling support
   - Better color consistency across all visualizations

7. **Bug Fixes**
   - Fixed RStudio display issues in `heatDiffM` (changed `invisible()` to `return()`)
   - Fixed custom_pairs index conversion in all multi-condition functions
   - Fixed missing `scPalette()` reference (replaced with `assignColors()`)
   - Fixed width/height swap in big_heatmap mode
   - Fixed pathway display in comparison plots to be comparison-specific

### Migration Notes

- **No breaking changes** - all existing code will continue to work
- `pheatmap` dependency removed in favor of `ComplexHeatmap`
- `dplyr` and `tidyr` are no longer required
- VennDetail is optional (only needed for `findCommonUniquePatterns` Venn diagrams)

## Quick Start

```r
library(CellDiff)

# Load example data
data(celllist)  # Contains two CellChat objects: NL (Normal Lung) and LS (Lung Scleroderma)

# Compare ligand-receptor pair contributions between conditions
ContriDiff(object.list = celllist, signaling = "TNF")

# Create a heatmap of differential signaling patterns
heatDiff(object.list = celllist, comparison = c(1, 2), measure = "sender")

# Visualize network differences with a chord diagram
netDiff(object.list = celllist, comparison = c(1, 2), signaling = "MIF")

# Compare signaling influence between conditions
scatterDiff(object.list = celllist, comparison = c(1, 2))

# Compare sender and receiver roles between conditions
scatterDiff2D(object.list = celllist, comparison = c(1, 2))

# Rank differential pathways
rankDiff(object.list = celllist)$plot

# Create an integrated network visualization
networkLRDiff(object.list = celllist, comparison = c(1, 2), pathways = c("MIF"))
```

## Comprehensive Tutorials

### 1. Working with Multiple Conditions

CellDiff provides specialized functions for comparing cell-cell communication across more than two conditions. Let's assume we have a list of CellChat objects for three conditions: WT, KO, and DKO.

```r
# First, create a named list of CellChat objects
cellchatlist <- data("cellchatlist")

# Set WT as the reference
reference <- 1  # Or use the name: reference = "WT"
```

#### Multi-condition Pathway Comparison

```r
# Compare signaling pathway differences across all conditions
pathway_results <- rankDiffM(
  object.list = cellchatlist,
  comparison_method = "all_vs_ref",  # Compare all conditions to the reference
  reference = "WT",                  # Set WT as reference
  cell.type.strategy = "union",      # Use all cell types (fills missing with 0)
  use_log2fc = TRUE,                 # Use log2 fold change for better visualization
  show_comparison_barplot = TRUE,    # Show faceted barplot
  show_comparison_heatmap = TRUE,    # Show heatmap
  show_comparison_slope = TRUE,      # Show slope plot
  show.all = FALSE                   # Show only significant pathways per comparison (default)
)

# View different visualization types
pathway_results$comparison_barplot
pathway_results$comparison_heatmap
pathway_results$comparison_slope_plot

# Get the top differential pathways
top_pathways <- pathway_results$top_paths

# Access results by comparison name (NEW!)
pathway_results$all_significant_paths_full[["KO vs WT"]]
pathway_results$all_significant_paths_full[["DKO vs WT"]]

# Show all pathways with visual indicators
pathway_results_all <- rankDiffM(
  object.list = cellchatlist,
  comparison_method = "all_vs_ref",
  reference = "WT",
  cell.type.strategy = "shared",     # Use only common cell types (safer)
  show_comparison_barplot = TRUE,
  show.all = TRUE                    # Show all pathways with fading and stars
)

# With show.all = TRUE:
# - All pathways appear in all comparison panels
# - Non-significant pathways are faded (alpha = 0.3)
# - Significant pathways have stars (*p<0.05, **p<0.01, ***p<0.001)
pathway_results_all$comparison_barplot

# Use custom pairs for specific comparisons
pathway_custom <- rankDiffM(
  object.list = cellchatlist,
  comparison_method = "custom_pairs",
  custom_comparisons = list(
    c("WT", "KO"),
    c("KO", "DKO"),
    c("WT", "DKO")
  ),
  show_comparison_barplot = TRUE,
  show_comparison_heatmap = TRUE
)

# View custom comparison results
pathway_custom$comparison_barplot  # Shows all custom pairs
pathway_custom$comparison_heatmap  # Heatmap of custom comparisons
```

#### Multi-condition Heatmap Visualization

```r
# Create a multi-condition heatmap comparing sender signaling
heatmap_result <- heatDiffM(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),       # Compare all three conditions
  reference = "WT",              # Set WT as reference
  measure = "sender",            # Focus on sender signaling
  use_log2fc = TRUE,             # Use log2 fold change
  show_values = TRUE,
  big_heatmap = TRUE,
  color.heatmap = c("blue", "white", "red")  # Custom colors
)

# The heatmap will display in RStudio automatically
# To save to PDF:
pdf("heatmap_output.pdf", width = 15, height = 8)
ComplexHeatmap::draw(heatmap_result)
dev.off()
```

#### Custom Comparison Pairs

```r
# Compare specific condition pairs (e.g., WT vs KO, KO vs DKO, WT vs DKO)
heatmap_custom <- heatDiffM(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),
  comparison_method = "custom_pairs",
  custom_comparisons = list(
    c("WT", "KO"),    # First comparison
    c("KO", "DKO"),   # Second comparison
    c("WT", "DKO")    # Third comparison
  ),
  measure = "sender",
  use_log2fc = TRUE,
  show_values = TRUE,
  big_heatmap = TRUE,
  color.heatmap = c("blue", "white", "red")
)
```

#### Multi-condition Scatter Plot Analysis

```r
# Compare sender-receiver roles across conditions
#scatter_result <- scatterDiffM(
 # object.list = cellchatlist,
#  comparison = c(1, 2, 3),       # Use all three conditions
#  reference = "WT",              # Set WT as reference
#  comparison_method = "all_vs_ref",  # Compare all to reference
#  show_group_labels = TRUE,      # Show group labels
#  convex_hull = FALSE             # Draw convex hulls around groups
#)

# For more detailed 2D visualization
scatter_result_2d <- scatterDiff2DM(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),       # Use all three conditions
  reference = "WT",              # Set WT as reference
  comparison_method = "all_vs_ref",  # Compare all to reference
  show.group.legend = TRUE,      # Show legend
  convex.hull = TRUE            # Draw convex hull
)
scatter_result_2d <- scatterDiff2DM(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),       # Use all three conditions
  reference = "WT",              # Set WT as reference
  comparison_method = "all_vs_ref",  # Compare all to reference
  show.group.legend = TRUE,      # Show legend
  convex.hull = TRUE,            # Draw convex hulls
  highlight.cells = c("Endo", "Mac", "ImmSC")  # Highlight specific cell types
)
```

#### L-R Pair Contribution Analysis

```r
# Compare L-R pair contributions for the WNT pathway across conditions
lr_result <- ContriDiffM(
  object.list = cellchatlist,
  signaling = "WNT",
  reference = "WT",
  stack.method = "side-by-side",  # Display conditions side by side
  show.heatmap = TRUE,            # Show heatmap visualization
  top.n = 10                      # Limit to top 10 L-R pairs
)

# The function returns a list of plots
lr_result$barplot
lr_result$heatmap

# Use custom pairs for L-R contribution analysis
lr_custom <- ContriDiffM(
  object.list = cellchatlist,
  signaling = "WNT",
  comparison_method = "custom_pairs",
  custom_comparisons = list(
    c("WT", "KO"),
    c("KO", "DKO")
  ),
  show.heatmap = TRUE,
  top.n = 10
)
```

#### Finding Common and Unique Signaling Patterns

```r
# Identify signaling patterns that are common or unique across conditions
# With comprehensive visualizations
patterns <- findCommonUniquePatterns(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),       # Compare all three conditions
  min_overlap = 2,               # Minimum number of conditions to consider a pattern common
  show.venn = TRUE,              # Show Venn diagram (default)
  show.upset = TRUE,             # Show UpSet plot
  show.barplot = TRUE,           # Show pattern count barplot
  show.network = TRUE,           # Show network visualization
  return_networks = TRUE,        # Return full network data
  return.data = TRUE,            # Return underlying data
  top.n = 20                     # Top 20 patterns for network
)

# Examine common patterns across conditions
head(patterns$common_patterns)

# Examine unique patterns for each condition
head(patterns$unique_patterns)

# View visualizations
patterns$plots$venn         # Venn diagram showing overlaps
patterns$plots$upset        # UpSet plot for complex intersections
patterns$plots$barplot      # Pattern counts per condition
patterns$plots$network_graph # Network of top interactions

# Access the VennDetail object for further analysis
venn_obj <- patterns$plots$venn_object
```

#### Comparing Multiple Cell-Cell Communications

```r
# Compare multiple CellChat objects to identify global changes
global_comparison <- compareCellChatsM(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),      # Compare all three conditions
  reference = "WT",             # Set WT as reference
  measure.type = "sender",      # Focus on sender signaling
  norm.method = "z-score",      # Normalize with z-score
  show.pathway = TRUE,          # Show pathway contributions
  show.heatmap = TRUE,          # Generate heatmap (ComplexHeatmap)
  label.timeline = TRUE,        # Add labels to timeline
  label.pathway = TRUE,         # Add labels to pathway plot
  label.top = 5,                # Label top 5 per condition
  return.data = TRUE            # Return underlying data
)

# View the results
global_comparison$timeline  # Timeline plot of signaling changes
global_comparison$pathway   # Pathway contribution changes
global_comparison$heatmap   # ComplexHeatmap with cell type annotations

# Access the underlying data
global_comparison$data$scores              # Signaling scores per condition
global_comparison$data$celltype_colors     # Consistent color scheme
```

### 2. Advanced Visualization Techniques

#### Network Difference Visualization

```r
# Visualize network differences for a specific pathway
net_diff <- netDiff(
  object.list = cellchatlist,
  comparison = c(1, 2),         # Compare first two conditions
  signaling = "WNT",            # Focus on WNT signaling
  use_normalized = TRUE,        # Normalize networks before comparison
  color.use = c("darkorange", "deepskyblue"),  # Custom colors
  highlight_sources = c("Tcell", "Bcell")  # Highlight specific sources
)

# Integrated network visualization showing pathways, cells, and LR pairs
network_viz <- networkLRDiff(
  object.list = cellchatlist,
  comparison = c(1, 2),         # Compare first two conditions
  pathways = top_pathways[1:3], # Use top differential pathways
  node.size.factor = 1.2,       # Larger nodes
  edge.width.factor = 1.5,      # Thicker edges
  label.cell = TRUE,            # Label cell types
  node.label.repel = TRUE       # Avoid label overlaps
)
```

#### Customizing Heatmaps

```r
# Create a custom heatmap with specific pathways
custom_heatmap <- heatDiffM(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),
  reference = "WT",
  measure = "both",             # Both sender and receiver
  signaling = top_pathways,     # Use top differential pathways
  transpose = TRUE,             # Transpose the heatmap
  use_log2fc = TRUE,            # Use log2 fold change
  color.heatmap = c("purple", "white", "green"),  # Custom colors
  cluster_rows = TRUE,          # Cluster rows
  cluster_cols = TRUE,          # Cluster columns
  show_values = TRUE,           # Show values in cells
  value_digits = 1,             # 1 decimal place
  show.heatmap_border = TRUE,   # Show borders
  save_plot = TRUE,             # Save the plot
  save_name = "custom_heatmap.pdf"  # Filename
)
```

### 3. Working with Large Datasets

For large datasets with many cell types or pathways, CellDiff provides options to filter and focus on the most relevant information:

```r
# Filter to include only significant pathways
filtered_results <- rankDiffM(
  object.list = cellchatlist,
  pThresh = 0.01,               # Stricter p-value threshold
  filter_min_change = 0.5,      # Higher fold change threshold
  top.n = 15                    # Limit to top 15 pathways
)

# Filter specific cell types for visualization
specific_heatmap <- heatDiffM(
  object.list = cellchatlist,
  filter_zeros = TRUE,          # Filter out pathways with no signal
  filter_thresh = 0.1           # Minimum signal threshold
)
```

### 4. Interpreting Results

#### Understanding Differential Pathways

```r
# Get detailed information about differential pathways
pathway_data <- rankDiffM(
  object.list = cellchatlist,
  comparison_method = "all_vs_ref",
  reference = "WT",
  return_top_paths = TRUE,
  return_data = TRUE
)

# Extract detailed information
pathway_stats <- pathway_data$data
significant_paths <- pathway_data$all_significant_paths

# Print statistics for top pathways
print(pathway_stats[pathway_stats$name %in% pathway_data$top_paths, ])
```

#### Cell Type Role Changes

```r
# Analyze how cell types change their roles
role_changes <- scatterDiff2DM(
  object.list = cellchatlist,
  comparison_method = "all_vs_ref",
  reference = "WT",
  return_data = TRUE
)

# Examine cells with the largest changes
arrows_data <- role_changes$arrows
top_changes <- arrows_data[order(arrows_data$total_change, decreasing = TRUE), ]
head(top_changes)

# Identify cells that change quadrants (e.g., from sender to receiver)
quadrant_changes <- arrows_data[arrows_data$quadrant %in% c("II", "IV"), ]
print(quadrant_changes)
```

## Frequently Asked Questions

### How do I handle missing cell types across conditions?
Use the `cell.type.strategy` parameter:
```r
# Use only cell types present in ALL conditions (safer, recommended)
rankDiffM(object.list = cellchatlist, cell.type.strategy = "shared")

# Use ALL cell types from any condition (fills missing with 0)
rankDiffM(object.list = cellchatlist, cell.type.strategy = "union")
```
The function will report which cell types are missing and how many:
```
Using union strategy: 10 cell types aligned
  - 8 cell types present in all objects
  - 2 cell types present in some objects (will use 0 for missing)
```

### How do I access results by comparison name instead of index?
All output lists now have meaningful names:
```r
results <- rankDiffM(object.list = cellchatlist, reference = "WT")

# Old way (still works)
results$all_significant_paths_full[[1]]

# New way (clearer!)
results$all_significant_paths_full[["KO vs WT"]]
results$all_significant_paths_full[["DKO vs WT"]]

# See all comparison names
names(results$all_significant_paths_full)
```

### How do I show all pathways with significance indicators?
Use the `show.all` parameter in `rankDiffM`:
```r
# Default: Only show pathways significant in each specific comparison
results <- rankDiffM(object.list = cellchatlist, show.all = FALSE)

# Show all pathways with visual indicators
results <- rankDiffM(object.list = cellchatlist, show.all = TRUE)
```
With `show.all = TRUE`:
- All pathways appear in all comparison panels
- Non-significant pathways are faded (alpha = 0.3)
- Significant pathways are opaque with stars (*p<0.05, **p<0.01, ***p<0.001)
- A subtitle explains the visualization

### How do I select a reference condition?
For multiple group comparisons, you can specify a reference using either the index or name:
```r
# By index (position in the list)
rankDiffM(object.list = cellchatlist, reference = 1)

# By name (if your list has named elements)
rankDiffM(object.list = cellchatlist, reference = "WT")
```

### How can I compare specific conditions within a larger dataset?
Use the `comparison` parameter to select specific conditions:
```r
# Compare only WT and DKO (assuming positions 1 and 3)
heatDiffM(object.list = cellchatlist, comparison = c(1, 3))

# Or by name
heatDiffM(object.list = cellchatlist, comparison = c("WT", "DKO"))
```

### How do I find the most significant pathway differences?
```r
# Use rankDiffM and extract top pathways
results <- rankDiffM(object.list = cellchatlist)
top_paths <- results$top_paths
```

### Can I export results to other formats?
Many functions have a `return_data = TRUE` option that provides the underlying data:
```r
# Get data for further processing
heatmap_data <- heatDiffM(object.list = cellchatlist, return_data = TRUE)
write.csv(heatmap_data$diff_matrix, "differential_heatmap_data.csv")
```

## Contact

For questions or issues, please contact: guokai8@gmail.com

Create an issue on the GitHub repository: [https://github.com/guokai8/CellDiff/issues](https://github.com/guokai8/CellDiff/issues)
