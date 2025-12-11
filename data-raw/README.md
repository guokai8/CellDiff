# Tutorial Data for CellDiff

This directory contains scripts for generating tutorial data for testing CellDiff functions.

## Files

- **create_tutorial_data.R** - Generates synthetic Seurat object with multiple conditions
- **test_tutorial.R** - Test script demonstrating the complete workflow
- **README.md** - This file

## Tutorial Dataset: pbmc_tutorial

### Quick Start

```r
# Load CellDiff and tutorial data
library(CellDiff)
data(pbmc_tutorial)

# Step 1: Create CellChat objects from Seurat
cellchat_list <- runCellChat(
  seurat_object = pbmc_tutorial,
  group.by = "condition",
  species = "human",
  cell.type.column = "cell_type"
)

# Step 2: Run comprehensive differential analysis
results <- compareCell(
  object.list = cellchat_list,
  reference = "WT",
  show_plots = c("barplot", "heatmap", "scatter")
)

# View results
results$pathway_analysis$comparison_barplot
results$summary
```

### Dataset Details

**Dimensions:**
- 1,500 cells (500 per condition)
- 2,000 genes
- 3 conditions: WT, KO, DKO
- 5 cell types: T_cells, B_cells, Monocytes, NK_cells, Dendritic

**Design:**
- 100 cells per cell type per condition
- Cell type-specific marker genes
- Synthetic ligand and receptor genes
- Simulated differential expression:
  - KO: Reduced ligand expression (30% of WT)
  - DKO: Further reduced ligands (10% of WT) + reduced receptors (40% of WT)

**Metadata columns:**
- `condition`: WT, KO, or DKO
- `cell_type`: T_cells, B_cells, Monocytes, NK_cells, or Dendritic
- `nCount_RNA`: Total UMI counts
- `nFeature_RNA`: Number of detected genes

## Regenerating the Data

To regenerate the tutorial data:

```r
# Run the data generation script
source("data-raw/create_tutorial_data.R")

# This will create data/pbmc_tutorial.rda
```

## Testing the Workflow

To test the complete CellDiff workflow:

```r
# Run the test script
source("data-raw/test_tutorial.R")

# This will:
# 1. Load the tutorial data
# 2. Create CellChat objects
# 3. Run differential analysis
# 4. Display results summary
```

## Use Cases

This tutorial data is perfect for:

1. **Learning CellDiff**: Understand the workflow without real data
2. **Testing new features**: Validate code changes
3. **Troubleshooting**: Debug issues with a known dataset
4. **Workshops**: Teach cell-cell communication analysis
5. **Documentation**: Create examples and vignettes

## Notes

- This is **synthetic data** for demonstration purposes
- Gene names are simplified (LIGAND_1, RECEPTOR_1, etc.)
- CellChat analysis may not detect many real pathways (expected behavior)
- Use real CellChatDB for actual biological analysis
- The data is designed to show the **workflow**, not biological insights

## Real Data Analysis

For real biological analysis:

1. Start with your own Seurat object with real gene names
2. Use the same workflow:
   ```r
   cellchat_list <- runCellChat(your_seurat, group.by = "condition", species = "human")
   results <- compareCell(cellchat_list, reference = "control")
   ```
3. CellChat will use the real ligand-receptor database
4. Results will reflect actual biological communication patterns

## Citation

If you use CellDiff in your research, please cite both:
- CellDiff package (this package)
- CellChat: Jin et al., Nature Communications, 2021
