#' Save heatmap to PDF safely
#'
#' This is a helper function to safely save ComplexHeatmap objects to PDF
#' without risking RStudio crashes. Use this after creating a heatmap with
#' any of the heatDiffM, rankDiffM, etc. functions.
#'
#' @param heatmap_obj A heatmap object (HeatmapList or result from heatDiffM)
#' @param output_file Path where PDF should be saved
#' @param width Width of PDF in inches (default: 15)
#' @param height Height of PDF in inches (default: 8)
#' @param open_after Logical, whether to open PDF after saving (default: TRUE)
#'
#' @return Path to saved file (invisibly)
#'
#' @examples
#' # After creating a heatmap:
#' hm <- heatDiffM(cellchat, comparison_method = "custom_pairs", ...)
#' 
#' # Save it safely:
#' save_heatmap(hm, "my_heatmap.pdf")
#'
#' @export
save_heatmap <- function(heatmap_obj, 
                         output_file, 
                         width = 15, 
                         height = 8,
                         open_after = TRUE) {
  
  # Validate inputs
  if (missing(output_file)) {
    stop("output_file is required")
  }
  
  if (!is.character(output_file)) {
    stop("output_file must be a character string")
  }
  
  # Handle both raw heatmap objects and result lists
  if (is.list(heatmap_obj) && !is.null(heatmap_obj$heatmap)) {
    # If it's a result list from return_data = TRUE
    hm_to_draw <- heatmap_obj$heatmap
    
    # Use draw_params if available
    if (!is.null(heatmap_obj$draw_params) && length(heatmap_obj$draw_params) > 0) {
      message("Saving heatmap with custom draw parameters...")
      pdf(output_file, width = width, height = height)
      ComplexHeatmap::draw(
        hm_to_draw,
        column_title = heatmap_obj$draw_params$column_title,
        column_title_gp = heatmap_obj$draw_params$column_title_gp
      )
      dev.off()
    } else {
      message("Saving heatmap to: ", output_file)
      pdf(output_file, width = width, height = height)
      ComplexHeatmap::draw(hm_to_draw)
      dev.off()
    }
  } else {
    # Direct heatmap object
    message("Saving heatmap to: ", output_file)
    pdf(output_file, width = width, height = height)
    ComplexHeatmap::draw(heatmap_obj)
    dev.off()
  }
  
  # Get absolute path
  abs_path <- normalizePath(output_file, mustWork = FALSE)
  
  # Check if file was created successfully
  if (file.exists(abs_path)) {
    file_size <- file.size(abs_path)
    message(sprintf("✓ Success! Saved %.1f KB to: %s", file_size/1024, abs_path))
    
    # Open file if requested (Mac only for now)
    if (open_after && .Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
      system2("open", abs_path)
      message("✓ Opening PDF...")
    }
  } else {
    warning("File may not have been created successfully")
  }
  
  invisible(abs_path)
}


#' Safe wrapper for heatDiffM that automatically saves to PDF
#'
#' This function wraps heatDiffM and automatically saves the result to PDF,
#' preventing RStudio crashes from trying to display large heatmaps.
#'
#' @param output_file Path where PDF should be saved (required)
#' @param width Width of PDF in inches (default: 15)
#' @param height Height of PDF in inches (default: 8)
#' @param open_after Logical, whether to open PDF after saving (default: TRUE)
#' @param ... All other parameters passed to heatDiffM
#'
#' @return Path to saved file (invisibly)
#'
#' @examples
#' heatDiffM_safe(
#'   output_file = "my_heatmap.pdf",
#'   object.list = cellchat,
#'   comparison = names(cellchat),
#'   comparison_method = "custom_pairs",
#'   custom_comparisons = list(c("CON_NS", "CUMS_NS")),
#'   measure = "sender"
#' )
#'
#' @export
heatDiffM_safe <- function(output_file, 
                           width = 15, 
                           height = 8,
                           open_after = TRUE,
                           ...) {
  
  if (missing(output_file)) {
    stop("output_file is required for heatDiffM_safe")
  }
  
  # Create heatmap object
  message("Creating heatmap...")
  hm <- heatDiffM(...)
  
  # Save it
  save_heatmap(hm, output_file, width, height, open_after)
}
