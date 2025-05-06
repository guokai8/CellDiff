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
* **Differential pathway analysis**: Identify signaling pathways that are significantly altered between conditions
* **Visualization of changes**: Compare ligand-receptor pair contributions with customizable plots
* **Heatmap visualization**: Visualize differential signaling patterns across cell types and pathways
* **Network visualization**: Generate chord diagrams showing changes in cell-cell communication
* **Sender-receiver analysis**: Analyze how cells change their roles as signaling senders or receivers
* **Integrated network visualization**: Create comprehensive networks showing pathways, cells, and ligand-receptor pairs

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
- dplyr
- tidyr
- reshape2
- circlize
- ComplexHeatmap
- cowplot
- igraph
- ggrepel

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
  use_log2fc = TRUE,                 # Use log2 fold change for better visualization
  show_comparison_barplot = TRUE,    # Show faceted barplot
  show_comparison_heatmap = TRUE,    # Show heatmap
  show_comparison_slope = TRUE       # Show slope plot
)

# View different visualization types
pathway_results$comparison_barplot
pathway_results$comparison_heatmap
pathway_results$comparison_slope_plot

# Get the top differential pathways
top_pathways <- pathway_results$top_paths
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
  big_heatmap=TRUE, 
  color.heatmap = c("blue", "white", "red")  # Custom colors
)

# The function returns a list with visualization and data
print(heatmap_result$heatmap)
```

#### Multi-condition Scatter Plot Analysis

```r
# Compare sender-receiver roles across conditions
scatter_result <- scatterDiffM(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),       # Use all three conditions
  reference = "WT",              # Set WT as reference
  comparison_method = "all_vs_ref",  # Compare all to reference
  show_group_labels = TRUE,      # Show group labels
  convex_hull = FALSE             # Draw convex hulls around groups
)

# For more detailed 2D visualization
scatter_result_2d <- scatterDiff2DM(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),       # Use all three conditions
  reference = "WT",              # Set WT as reference
  comparison_method = "all_vs_ref",  # Compare all to reference
  show.group.legend = TRUE,      # Show legend
  convex.hull = TRUE,            # Draw convex hulls
  highlight.cells = c("Tcell", "Bcell")  # Highlight specific cell types
)
```

#### L-R Pair Contribution Analysis

```r
# Compare L-R pair contributions for the TNF pathway across conditions
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
```

#### Finding Common and Unique Signaling Patterns

```r
# Identify signaling patterns that are common or unique across conditions
patterns <- findCommonUniquePatterns(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),       # Compare all three conditions
  min_overlap = 2,               # Minimum number of conditions to consider a pattern common
  return_networks = TRUE         # Return full network data
)

# Examine common patterns across conditions
head(patterns$common_patterns)

# Examine unique patterns for each condition
head(patterns$unique_patterns)
```

#### Comparing Multiple Cell-Cell Communications

```r
# Compare multiple CellChat objects to identify global changes
global_comparison <- compareCellChatsM(
  object.list = cellchatlist,
  comparison = c(1, 2, 3),      # Compare all three conditions
  reference = "WT",             # Set WT as reference
  measure.type = "sender",      # Focus on sender signaling
  show.heatmap = TRUE,          # Generate heatmap
  show.barplot = TRUE           # Generate barplot
)

# View the results
global_comparison$timeline  # Changes in signaling across conditions
global_comparison$heatmap   # Heatmap of differences
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
  link_style = "bezier",        # Use bezier curves
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
  sources.use = c("Tcell", "Bcell", "Mono"),  # Filter sender cells
  targets.use = c("Fibroblast", "Endothelial"),  # Filter receiver cells
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
