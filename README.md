# CellDiff
Enhances visualization and statistical analysis of cell-cell communication differences between experimental groups.

CellDiff is an R package that extends the functionality of CellChat by providing specialized tools for differential cell-cell communication analysis between multiple conditions. While CellChat enables the inference and analysis of cell-cell communication networks from single-cell RNA sequencing data, CellDiff focuses on identifying and visualizing the differences in these communication patterns across biological conditions (e.g., normal vs. disease, treated vs. untreated).
Features

* Differential pathway analysis: Identify signaling pathways that are significantly altered between conditions
* Visualization of changes: Compare ligand-receptor pair contributions with customizable plots
* Heatmap visualization: Visualize differential signaling patterns across cell types and pathways
* Network visualization: Generate chord diagrams showing changes in cell-cell communication
* Sender-receiver analysis: Analyze how cells change their roles as signaling senders or receivers
* Integrated network visualization: Create comprehensive networks showing pathways, cells, and ligand-receptor pairs

## Installation
### Install devtools if not already installed
```{r}
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install CellDiff from GitHub
devtools::install_github("guokai8/CellDiff")
```

Dependencies
CellDiff requires the following packages: CellChat (>= 1.1.0), ggplot2, dplyr,tidyr,reshape2,circlize,ComplexHeatmap,cowplot,igraph,ggrepel

# Quick Start
```{r}
library(CellDiff)

# Load example data
data(celllist)  # Contains two CellChat objects: NL (Normal Lung) and LS (Lung Scleroderma)

# Compare ligand-receptor pair contributions between conditions
ContriDiff(object.list = celllist, signaling = "TNF")

# Create a heatmap of differential signaling patterns
heatDiff(celllist, comparison = c(1, 2), measure = "sender")

# Visualize network differences with a chord diagram
netDiff(celllist, comparison = c(1, 2), signaling = "TNF")

# Compare sender and receiver roles between conditions
scatterDiff2D(celllist, comparison = c(1, 2))

# Create an integrated network visualization
networkLRDiff(celllist, comparison = c(1, 2), pathways = c("TNF", "TGFb"))
```

### Contact
For questions or issues, please contact:guokai8@gmail.com
### Create an issue on the GitHub repository: 
https://github.com/guokai8/CellDiff/issues
