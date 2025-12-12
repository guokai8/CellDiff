# How to Fix runCellChat to Match Your Manual Pipeline

## The Problem

`runCellChat` produces different results than your manual CellChat pipeline, even though the steps look the same.

## Most Likely Cause: Parameter Defaults

The issue is probably in **`computeCommunProb`** parameters. CellChat's defaults may differ from what `runCellChat` explicitly sets.

### Current runCellChat Parameters (line 172 in runCellChat.R):

```r
raw.use = TRUE,        # <-- This is the likely culprit!
```

### Your Manual Pipeline:

```r
wt <- computeCommunProb(wt)  # Uses CellChat defaults
```

## To Diagnose the Exact Difference:

Run this R code with both objects loaded:

```r
library(CellChat)
library(CellDiff)

# Load both
data(celllist)  # From runCellChat
# load your 'cells' object (manual pipeline)

# Compare specific pathway
ptn_celllist <- sum(celllist$WT@netP$prob[,,'PTN'], na.rm=TRUE)
ptn_cells <- sum(cells$WT@netP$prob[,,'PTN'], na.rm=TRUE)

cat("PTN pathway strength:\n")
cat(sprintf("  celllist (runCellChat): %.6e\n", ptn_celllist))
cat(sprintf("  cells (manual): %.6e\n", ptn_cells))
cat(sprintf("  Ratio: %.2f\n", ptn_cells/ptn_celllist))

# If ratio is close to 1.0, parameters are the same
# If ratio is different, parameters differ!
```

## To Fix runCellChat:

### Option 1: Change raw.use default

If CellChat's default for `raw.use` is `FALSE`:

**Edit line 172 in R/runCellChat.R:**
```r
# Change FROM:
raw.use = TRUE,

# Change TO:
raw.use = FALSE,
```

### Option 2: Check other parameters

Other parameters to verify:

1. **type** (line 170): Default is `"triMean"` - check if CellChat uses this
2. **trim** (line 171): Default is `0.1` - verify this matches
3. **population.size** (line 169): Default is `TRUE` - verify this matches

## Quick Fix Steps:

1. **Find the difference:**
   ```r
   # In R, check CellChat's actual defaults
   library(CellChat)
   args(computeCommunProb)
   ```

2. **Update runCellChat:**
   ```r
   # Edit R/runCellChat.R
   # Change the default values to match CellChat's defaults
   ```

3. **Rebuild and test:**
   ```r
   library(devtools)
   document()
   load_all(".")

   # Test with your data
   cellchat_new <- runCellChat(your_seurat, group.by="condition", species="mouse")
   ```

4. **Compare results:**
   ```r
   # Should now match your manual pipeline
   sum(cellchat_new$WT@netP$prob[,,'PTN'])
   ```

## Most Likely Solution:

Based on typical CellChat usage, the fix is probably:

**Change line 172 in R/runCellChat.R:**
```r
raw.use = FALSE,  # Changed from TRUE
```

This is the most common difference between different CellChat workflows.

## Alternative: Make runCellChat Match Your Exact Workflow

If you want runCellChat to exactly match your pipeline, explicitly set all parameters:

```r
cellchat_list <- runCellChat(
  seurat_object = your_seurat,
  group.by = "condition",
  species = "mouse",
  # Explicitly set to match your pipeline
  raw.use = FALSE,  # or whatever you use
  type = "triMean",
  trim = 0.1,
  population.size = TRUE,
  min.cells = 10
)
```

## To Verify the Fix Works:

```r
# After fixing, both should give same results
result1 <- rankDiff(celllist, return.data=TRUE)  # Old
result2 <- rankDiff(cells, return.data=TRUE)     # Your manual

# Compare
all.equal(result1$data, result2$data)
```

---

**TL;DR**: The most likely fix is changing `raw.use = TRUE` to `raw.use = FALSE` in line 172 of `R/runCellChat.R`.
