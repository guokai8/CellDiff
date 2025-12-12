# runCellChat Parameters Fixed

## Problem Identified

`runCellChat` was using different default parameters than CellChat's own defaults, causing different results compared to manual CellChat pipelines.

## Parameters Fixed

### 1. population.size
- **CellChat default**: `FALSE`
- **Old runCellChat default**: `TRUE` ❌
- **New runCellChat default**: `FALSE` ✅

**Impact**: This affects how communication strength is calculated. When `FALSE`, communication probability is independent of cell population size.

### 2. thresh.pc
- **CellChat default**: `0`
- **Old runCellChat default**: `0.1` ❌
- **New runCellChat default**: `0` ✅

**Impact**: This is the threshold for percentage of cells expressing a gene in `identifyOverExpressedGenes()`. Setting to 0 means no filtering by percent cells.

### 3. thresh.fc
- **CellChat default**: `0`
- **Old runCellChat default**: `0.1` ❌
- **New runCellChat default**: `0` ✅

**Impact**: This is the fold change threshold for overexpressed genes. Setting to 0 means no fold change filtering.

## Changes Made

**File**: `R/runCellChat.R`

**Line 169**: Changed `population.size = TRUE` → `population.size = FALSE`
**Line 179**: Changed `thresh.pc = 0.1` → `thresh.pc = 0`
**Line 180**: Changed `thresh.fc = 0.1` → `thresh.fc = 0`

## Testing

After these changes, `runCellChat` should produce identical results to your manual CellChat pipeline:

```r
library(devtools)
load_all(".")

# Create CellChat objects
cellchat_list <- runCellChat(
  seurat_object = your_seurat,
  group.by = "condition",
  species = "mouse"
)

# Should now match your manual pipeline
# Compare pathways
cellchat_list$WT@netP$pathways
cells$WT@netP$pathways  # Your manual version

# Compare specific pathway strengths
sum(cellchat_list$WT@netP$prob[,,'PTN'])
sum(cells$WT@netP$prob[,,'PTN'])  # Should be identical now!
```

## Backwards Compatibility

**BREAKING CHANGE**: If you were using `runCellChat` before, the results will now be different (but more correct, matching CellChat defaults).

If you want the old behavior back:
```r
cellchat_list <- runCellChat(
  seurat_object = your_seurat,
  group.by = "condition",
  species = "mouse",
  population.size = TRUE,  # Explicitly set old value
  thresh.pc = 0.1,         # Explicitly set old value
  thresh.fc = 0.1          # Explicitly set old value
)
```

## Summary

✅ **Fixed**: `runCellChat` now uses CellChat's actual defaults
✅ **Result**: Should produce identical results to manual CellChat pipeline
✅ **Impact**: More accurate and consistent with CellChat documentation

---

**Date**: 2025-12-12
**Status**: Fixed and documented
