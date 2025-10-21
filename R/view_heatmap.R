#' View ComplexHeatmap safely in RStudio Plots pane
#'
#' This function displays a ComplexHeatmap object in RStudio's Plots pane
#' without crashing. It works by rendering to PNG first, then displaying
#' the image.
#'
#' @param heatmap_obj A ComplexHeatmap object or result from heatDiffM
#' @param width Width of the plot in pixels (default: 2400)
#' @param height Height of the plot in pixels (default: 1200)
#' @param res Resolution in DPI (default: 150)
#'
#' @return NULL (displays plot as side effect)
#'
#' @examples
#' # After creating a heatmap:
#' hm <- heatDiffM(cellchat, comparison_method = "custom_pairs", ...)
#' 
#' # View it safely in RStudio:
#' view_heatmap(hm)
#' 
#' # Or with custom size:
#' view_heatmap(hm, width = 3000, height = 1500)
#'
#' @export
view_heatmap <- function(heatmap_obj, 
                         width = 2400, 
                         height = 1200,
                         res = 150) {
  
  # Check if png package is available
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package 'png' is required. Install it with: install.packages('png')")
  }
  
  # Handle both raw heatmap objects and result lists
  if (is.list(heatmap_obj) && !is.null(heatmap_obj$heatmap)) {
    hm_to_draw <- heatmap_obj$heatmap
    has_draw_params <- !is.null(heatmap_obj$draw_params) && 
                       length(heatmap_obj$draw_params) > 0
  } else {
    hm_to_draw <- heatmap_obj
    has_draw_params <- FALSE
  }
  
  # Create temporary PNG file
  temp_png <- tempfile(fileext = ".png")
  
  # Render to PNG
  png(temp_png, width = width, height = height, res = res)
  
  if (has_draw_params) {
    ComplexHeatmap::draw(
      hm_to_draw,
      column_title = heatmap_obj$draw_params$column_title,
      column_title_gp = heatmap_obj$draw_params$column_title_gp
    )
  } else {
    ComplexHeatmap::draw(hm_to_draw)
  }
  
  dev.off()
  
  # Read and display in RStudio Plots pane
  img <- png::readPNG(temp_png)
  grid::grid.newpage()
  grid::grid.raster(img)
  
  # Cleanup
  unlink(temp_png)
  
  message("✓ Heatmap displayed in Plots pane")
  message("Tip: Click 'Zoom' for full-size view")
  
  invisible(NULL)
}


#' Create and view heatmap in one step
#'
#' Wrapper that creates a heatmap with heatDiffM and immediately displays it
#' in RStudio's Plots pane.
#'
#' @param ... All parameters passed to heatDiffM
#' @param view Logical, whether to display in Plots pane (default: TRUE)
#' @param width Width for viewing in pixels (default: 2400)
#' @param height Height for viewing in pixels (default: 1200)
#'
#' @return The heatmap object (invisibly)
#'
#' @examples
#' heatDiffM_view(
#'   object.list = cellchat,
#'   comparison = names(cellchat),
#'   comparison_method = "custom_pairs",
#'   custom_comparisons = list(c("CON_NS", "CUMS_NS")),
#'   measure = "sender"
#' )
#'
#' @export
heatDiffM_view <- function(..., 
                           view = TRUE,
                           width = 2400,
                           height = 1200) {
  
  # Create heatmap with safe device
  options(device = function() { pdf(NULL) })
  
  hm <- heatDiffM(...)
  
  options(device = NULL)
  
  # View if requested
  if (view) {
    view_heatmap(hm, width = width, height = height)
  }
  
  invisible(hm)
}
