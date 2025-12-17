# Changelog

All notable changes to the CellDiff package will be documented in this file.

## Version 0.1.9 (2025-12-17)

### Added

#### P-value Adjustment Control
- **NEW `adjust.pval` parameter** for optional multiple testing correction in `rankDiff` and `rankDiffM`
  - `adjust.pval = FALSE` (default): No p-value adjustment, maintains backward compatibility
  - `adjust.pval = TRUE`: Applies multiple testing correction
- **NEW `adjust.method` parameter** to specify correction method
  - Default: `"BH"` (Benjamini-Hochberg/FDR)
  - Supports all `p.adjust()` methods: "bonferroni", "holm", "hochberg", "hommel", "BY", "fdr", "none"
- Users now have full control over when and how to apply p-value adjustment
- Both `rankDiff` and `rankDiffM` have consistent API for p-value adjustment

### Changed

- **`rankDiff`**: Previously did not adjust p-values; now offers optional adjustment via `adjust.pval` parameter
- **`rankDiffM`**: Previously always adjusted p-values with BH method; now adjustment is optional via `adjust.pval` parameter
- Default behavior (`adjust.pval = FALSE`) maintains backward compatibility for both functions

### Technical Details

**Files Modified:**
- `R/rankDiff.R`:
  - Lines 20-21: Added parameter documentation for `adjust.pval` and `adjust.method`
  - Lines 31-32: Added parameters to function signature
  - Lines 255-260: Conditional p-value adjustment logic
- `R/rankDiffM.R`:
  - Lines 53-54: Added parameter documentation for `adjust.pval` and `adjust.method`
  - Lines 107-115: Added example of using p-value adjustment
  - Lines 133-134: Added parameters to function signature
  - Lines 497-502: Conditional p-value adjustment logic
- `DESCRIPTION`: Updated version to 0.1.9 and date to 2025-12-17

### Usage Examples

```r
# Without p-value adjustment (default, backward compatible)
rankDiff(cellchat.list, comparison = c(1, 2))
rankDiffM(cellchat.list, comparison_method = "all_vs_ref")

# With Benjamini-Hochberg adjustment
rankDiff(cellchat.list, comparison = c(1, 2), adjust.pval = TRUE)
rankDiffM(cellchat.list, adjust.pval = TRUE, adjust.method = "BH")

# With Bonferroni correction (more conservative)
rankDiff(cellchat.list, comparison = c(1, 2),
         adjust.pval = TRUE, adjust.method = "bonferroni")
```

## Version 0.1.4 (2025-12-10)

### Added

#### Cell Type Alignment
- **NEW `cell.type.strategy` parameter** for handling missing cell types across conditions
  - `"shared"`: Uses only cell types present in all conditions (safer, more conservative)
  - `"union"`: Uses all cell types from any condition, filling missing data with zeros
  - Automatically reports alignment statistics
  - Warns about missing cell types and their distribution
  - Available in: `rankDiffM`, `heatDiffM`, `scatterDiff2DM`, `ContriDiffM`

#### Advanced Comparison Visualization (rankDiffM)
- **NEW `show.all` parameter** for comprehensive pathway visualization
  - `show.all = FALSE` (default): Shows only pathways significant in each specific comparison
  - `show.all = TRUE`: Shows all pathways with visual indicators:
    - Faded bars (alpha=0.3) for non-significant pathway-comparison combinations
    - Opaque bars (alpha=1.0) for significant combinations
    - Significance stars (*p<0.05, **p<0.01, ***p<0.001)
    - Informative subtitle explaining the visualization

#### Named List Output
- **All output lists now have meaningful comparison names**
  - Access results by comparison name: `results$all_significant_paths_full[["KO vs WT"]]`
  - Instead of numeric index: `results$all_significant_paths_full[[1]]`
  - Applies to: `all_significant_paths_full`, `all_significant_paths`, `data`, `plots`, `heatmaps`
  - Improves code readability and reduces errors

#### Debug Mode
- Added comprehensive debug output for troubleshooting
  - Enable with: `options(CellDiff.debug = TRUE)` (package-wide option)
  - Shows pathway selection logic
  - Shows star assignment decisions
  - Shows comparison-pair matching
  - Does not affect other functions (isolated to CellDiff package)

#### Wrapper Functions
- **NEW** `runCellChat`: Create CellChat objects from Seurat
  - Takes Seurat object directly as input
  - Automatically splits by condition and creates CellChat objects
  - Runs complete CellChat workflow (preprocessing, inference, aggregation)
  - Returns named list of ready-to-use CellChat objects
  - Simplifies CellChat setup for multi-condition experiments
  - Users can then apply any CellDiff functions for differential analysis

- **NEW** `compareCell`: Comprehensive differential analysis wrapper
  - One-step workflow running `rankDiffM`, `heatDiffM`, `scatterDiff2DM`, and `networkLRDiff`
  - Automatic cell type alignment
  - Organized results with meaningful comparison names
  - Integrated summary statistics
  - Informative progress messages
  - Simplified API for common analysis workflows

### Fixed

#### Comparison-Specific Pathway Display
- **FIXED**: Pathways now only appear in comparison panels where they're actually significant when `show.all = FALSE`
- **Previous behavior**: All significant pathways from any comparison appeared in all panels
- **New behavior**: Each comparison panel shows only its own significant pathways
- Example: If PTN is only significant in "KO vs WT", it now only appears in that panel

#### Star Assignment Bug
- **FIXED**: Stars now only appear for pathways passing BOTH p-value AND fold-change thresholds
- **Previous behavior**: Stars appeared for any pathway with p-value < threshold, ignoring fold-change
- **New behavior**: Stars only appear if pathway is in `all_significant_paths_full` (passed complete criteria)
- Affects both barplot and heatmap visualizations

### Changed

- Improved code documentation with more detailed parameter descriptions
- Enhanced error messages for better troubleshooting
- Updated README.md with comprehensive examples of new features
- Added FAQ section for common questions about new features

### Technical Details

**Files Modified:**
- `R/rankDiffM.R`:
  - Lines 784-787: Added debug output
  - Lines 834-844: Fixed show.all = FALSE to use comparison-specific pathways
  - Line 993: Fixed star assignment for barplots (check `all_significant_paths_full`)
  - Line 1190: Fixed star assignment for heatmaps (check `all_significant_paths_full`)
  - Lines 758-770: Added named list creation
- `R/runCellChat.R`: **NEW** wrapper function to create CellChat objects from Seurat
  - Lines 59-197: Main wrapper function implementation
  - Returns named list of CellChat objects
- `R/compareCell.R`: **NEW** wrapper function for comprehensive differential analysis
  - Lines 71-287: Main wrapper function implementation
  - Lines 294-321: Print method for CellDiffAnalysis objects
- `README.md`: Updated with new features, wrapper function examples
- `man/rankDiffM.Rd`: Auto-generated documentation updates
- `man/runCellChat.Rd`: **NEW** auto-generated documentation for Seurat wrapper
- `man/compareCell.Rd`: **NEW** auto-generated documentation for differential analysis wrapper

### Migration Guide

**No breaking changes** - all existing code will continue to work.

#### New Features to Adopt

1. **Use cell type alignment strategies**
   ```r
   # Old (still works, uses default "shared")
   rankDiffM(object.list = cellchatlist)

   # New (explicit strategy)
   rankDiffM(object.list = cellchatlist, cell.type.strategy = "union")
   ```

2. **Access results by name**
   ```r
   # Old (still works)
   results$all_significant_paths_full[[1]]

   # New (clearer)
   results$all_significant_paths_full[["KO vs WT"]]
   ```

3. **Use show.all for comprehensive visualization**
   ```r
   # Show all pathways with significance indicators
   rankDiffM(object.list = cellchatlist, show.all = TRUE)
   ```

4. **Use wrapper functions for streamlined workflow**
   ```r
   # Create CellChat objects from Seurat
   cellchat_list <- runCellChat(
     seurat_object = pbmc,
     group.by = "condition",
     species = "human"
   )

   # Then run differential analysis
   results <- compareCell(
     object.list = cellchat_list,
     reference = "WT",
     cell.type.strategy = "union",
     show.all = TRUE,
     show_plots = c("barplot", "heatmap", "scatter")
   )

   # Or use individual CellDiff functions
   rankDiffM(object.list = cellchat_list, reference = "WT")
   heatDiffM(object.list = cellchat_list, reference = "WT")
   ```

#### Behavior Changes to Be Aware Of

1. **show.all = FALSE now shows comparison-specific pathways**
   - Each comparison panel shows different pathways
   - Pathways only appear where they're significant
   - This is more accurate but may look different from previous versions

2. **Stars are now more strict**
   - Only appear for pathways passing both p-value AND fold-change thresholds
   - You may see fewer stars, but they're more meaningful

---

## Version 0.1.3 (2025-05-02)

### Bug Fixes

#### Critical Fixes
- **heatDiffM**: Fixed RStudio display issue
  - Changed `invisible(final_heatmap)` to `return(final_heatmap)` to match working heatDiff pattern
  - Heatmaps now display correctly in RStudio without crashes or workarounds

- **heatDiff**: Fixed missing function reference
  - Replaced non-existent `scPalette()` with `assignColors()` from utils.R

#### custom_pairs Support
Fixed index conversion bugs in all multi-condition comparison functions:

- **heatDiffM**:
  - Fixed custom_pairs index conversion (lines 384-387)
  - Fixed width/height swap in big_heatmap mode (lines 617-618)
  - Fixed legend breaks for edge cases (lines 502-534)
  - Fixed duplicate heatmap names (line 605)

- **rankDiffM**:
  - Fixed custom_pairs index conversion (lines 219-222)

- **ContriDiffM**:
  - Added custom_pairs comparison method support
  - Implemented proper index conversion (lines 142-179)

- **scatterDiff2DM**:
  - Fixed custom_pairs index conversion (lines 251-254)
  - Moved quadrant labels to corners for better visibility (lines 810-850)

### Improvements

- **Color Management**:
  - Consolidated color palette functions
  - Proper use of `assignColors()` throughout package

- **Documentation**:
  - Updated README.md with custom_pairs examples
  - Added INSTALLATION.md with setup instructions
  - Cleaned up all temporary test and debug files

### Package Cleanup

Removed all temporary files:
- Test scripts (test_*.R, debug_*.R)
- Temporary documentation (CHEAT_SHEET.txt, COMPLETE_GUIDE.md, etc.)
- Test PDF outputs
- Test data files

Package is now clean and ready for installation.

### Technical Details

**Index Conversion Pattern** (Applied to all custom_pairs implementations):
```r
# BEFORE (broken):
comparison_pairs[[i]] <- c(pair[1], pair[2])  # Stored object indices

# AFTER (fixed):
pos1 <- which(comparison == pair[1])  # Convert to positions in comparison vector
pos2 <- which(comparison == pair[2])
comparison_pairs[[i]] <- c(pos1, pos2)
```

**Return Statement Fix** (heatDiffM.R line 1057):
```r
# BEFORE:
invisible(final_heatmap)  # Prevented display in RStudio

# AFTER:
return(final_heatmap)  # Displays correctly in RStudio
```

### Known Issues

None at this time.

### Upgrade Instructions

If upgrading from a previous version:
```r
# Remove old version
remove.packages("CellDiff")

# Install new version
devtools::install_github("guokai8/CellDiff")
```

---

## Previous Versions

### Version 0.1.2
- Initial release with multi-condition comparison support
- Known issues with custom_pairs comparison method
- RStudio display issues

### Version 0.1.1
- Beta release

### Version 0.1.0
- Initial development version
