# Changelog

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
