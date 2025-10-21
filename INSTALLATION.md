# CellDiff Installation Guide

## Package Ready for Installation

The CellDiff package has been cleaned and prepared for installation. All temporary files, test scripts, and debug code have been removed.

## Quick Installation

### From GitHub
```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install CellDiff from GitHub
devtools::install_github("guokai8/CellDiff")
```

### From Local Source
```r
# If you have the source code locally
devtools::install("path/to/CellDiff")

# Or from the built tar.gz file
install.packages("path/to/CellDiff_0.1.3.tar.gz", repos = NULL, type = "source")
```

## What Was Fixed

### 1. Core Function Improvements
- **heatDiffM**: Fixed return statement to display correctly in RStudio
  - Changed from `invisible(final_heatmap)` to `return(final_heatmap)`
  - Fixed custom_pairs index conversion bug
  - Fixed width/height swap in big_heatmap mode
  - Fixed legend breaks for edge cases

- **heatDiff**: Fixed missing scPalette function
  - Replaced with `assignColors()` from utils.R

- **rankDiffM**: Fixed custom_pairs index conversion

- **ContriDiffM**: Added custom_pairs comparison method support

- **scatterDiff2DM**: Fixed custom_pairs index conversion and quadrant label positioning

### 2. Files Removed
All temporary test, debug, and documentation files have been removed:
- Test scripts (test_*.R, debug_*.R, etc.)
- Temporary documentation (CHEAT_SHEET.txt, COMPLETE_GUIDE.md, etc.)
- Test PDF outputs
- Test data files (cellchat_objects_list.rds)

### 3. Files Kept
Essential package files:
- **R/** - All core function files (fixed and working)
- **man/** - Documentation files
- **data/** - Package data (cellchatlist.rda)
- **vignettes/** - Package vignettes
- **DESCRIPTION** - Package metadata
- **NAMESPACE** - Exported functions
- **README.md** - Updated with custom_pairs examples
- **LICENSE** - GPL-3 license

## Using the Fixed Functions

### Basic Usage
```r
library(CellDiff)

# Load your CellChat objects
cellchat_list <- list(WT = cellchat_wt, KO = cellchat_ko, DKO = cellchat_dko)

# The functions now work correctly in RStudio
heatDiffM(
  object.list = cellchat_list,
  comparison_method = "custom_pairs",
  custom_comparisons = list(
    c("WT", "KO"),
    c("KO", "DKO"),
    c("WT", "DKO")
  ),
  measure = "sender",
  use_log2fc = TRUE,
  show_values = TRUE,
  big_heatmap = TRUE,
  color.heatmap = c("blue", "white", "red")
)
```

### Comparison Methods
The multi-condition functions now support three comparison methods:

1. **all_vs_ref**: Compare all conditions to a reference
   ```r
   heatDiffM(object.list, comparison_method = "all_vs_ref", reference = "WT")
   ```

2. **all_vs_all**: All pairwise comparisons
   ```r
   heatDiffM(object.list, comparison_method = "all_vs_all")
   ```

3. **custom_pairs**: Specific comparison pairs (NEW - FIXED!)
   ```r
   heatDiffM(
     object.list,
     comparison_method = "custom_pairs",
     custom_comparisons = list(c("WT", "KO"), c("KO", "DKO"))
   )
   ```

## Optional Helper Functions

The package includes optional helper functions for advanced use cases:
- `view_heatmap()` - View heatmaps via PNG rendering (useful for very complex plots)
- `save_heatmap()` - Save heatmaps to PDF safely

These are not exported by default but are available in R/view_heatmap.R and R/save_heatmap.R if needed.

## Verification

To verify the installation:
```r
library(CellDiff)

# Check that functions are available
ls("package:CellDiff")

# Test with example data
data(celllist)
heatDiff(object.list = celllist, comparison = c(1, 2), measure = "sender")
```

## Support

For issues or questions:
- Email: guokai8@gmail.com
- GitHub Issues: https://github.com/guokai8/CellDiff/issues

## Version Information

**Current Version**: 0.1.3
**Date**: 2025-05-02
**Built Package**: CellDiff_0.1.3.tar.gz (23MB)
