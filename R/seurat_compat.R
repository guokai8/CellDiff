#' Extract Data from Seurat Object (v4/v5 Compatible)
#'
#' @description
#' Internal helper function to extract data from Seurat objects that works with
#' both Seurat v4 and v5. Handles the differences in data access methods and
#' assay structures between versions.
#'
#' @param seurat_object A Seurat object (v4 or v5)
#' @param assay Name of assay to use (default: "RNA")
#' @param layer_or_slot For v5: layer name ("data", "counts", "scale.data").
#'   For v4: slot name ("data", "counts", "scale.data"). Default: "data"
#'
#' @return A sparse matrix of expression data
#'
#' @keywords internal
extractSeuratData <- function(seurat_object,
                               assay = "RNA",
                               layer_or_slot = "data") {

  # Check Seurat version
  seurat_version <- packageVersion("Seurat")
  is_v5 <- seurat_version >= "5.0.0"

  # Check if assay exists
  if (!assay %in% names(seurat_object@assays)) {
    stop("Assay '", assay, "' not found in Seurat object. Available assays: ",
         paste(names(seurat_object@assays), collapse = ", "))
  }

  # Get assay
  assay_obj <- seurat_object[[assay]]

  # Check assay class (Assay5 vs Assay)
  is_assay5 <- inherits(assay_obj, "Assay5")

  if (is_v5 && is_assay5) {
    # Seurat v5 with Assay5 - use layer parameter
    tryCatch({
      # Suppress deprecation warning
      suppressWarnings({
        data <- Seurat::GetAssayData(seurat_object,
                                     layer = layer_or_slot,
                                     assay = assay)
      })
      return(data)
    }, error = function(e) {
      # Fallback: try accessing layer directly
      tryCatch({
        if (layer_or_slot %in% names(assay_obj@layers)) {
          return(assay_obj@layers[[layer_or_slot]])
        } else {
          stop("Layer '", layer_or_slot, "' not found. Available layers: ",
               paste(names(assay_obj@layers), collapse = ", "))
        }
      }, error = function(e2) {
        stop("Failed to extract data from Seurat v5 object: ", e2$message)
      })
    })
  } else {
    # Seurat v4 or v5 with old Assay - use slot parameter
    tryCatch({
      # Suppress deprecation warning for v5
      suppressWarnings({
        data <- Seurat::GetAssayData(seurat_object,
                                     slot = layer_or_slot,
                                     assay = assay)
      })
      return(data)
    }, error = function(e) {
      # Fallback: try direct slot access
      tryCatch({
        slot_name <- layer_or_slot
        if (slot_name == "scale.data") slot_name <- "scale.data"
        if (slot_name %in% slotNames(assay_obj)) {
          return(slot(assay_obj, slot_name))
        } else {
          stop("Slot '", layer_or_slot, "' not found in assay")
        }
      }, error = function(e2) {
        stop("Failed to extract data from Seurat v4 object: ", e2$message)
      })
    })
  }
}


#' Check and Report Seurat Version Compatibility
#'
#' @description
#' Internal helper function that checks Seurat version and reports compatibility
#' information. Used for debugging and user information.
#'
#' @param seurat_object A Seurat object
#' @param verbose Logical. Print version information?
#'
#' @return List with version information
#'
#' @keywords internal
checkSeuratVersion <- function(seurat_object, verbose = TRUE) {

  # Get version info
  seurat_version <- packageVersion("Seurat")
  is_v5 <- seurat_version >= "5.0.0"
  is_v4 <- seurat_version >= "4.0.0" && seurat_version < "5.0.0"

  # Check assay class
  default_assay <- Seurat::DefaultAssay(seurat_object)
  assay_obj <- seurat_object[[default_assay]]
  is_assay5 <- inherits(assay_obj, "Assay5")

  # Determine compatibility mode
  if (is_v5 && is_assay5) {
    mode <- "Seurat v5 (Assay5)"
    data_access <- "layer"
  } else if (is_v5 && !is_assay5) {
    mode <- "Seurat v5 (legacy Assay)"
    data_access <- "slot"
  } else if (is_v4) {
    mode <- "Seurat v4"
    data_access <- "slot"
  } else {
    mode <- "Seurat v3 or earlier"
    data_access <- "slot"
  }

  if (verbose) {
    message("Seurat compatibility mode: ", mode)
    message("Data access method: ", data_access)
  }

  return(list(
    version = as.character(seurat_version),
    is_v5 = is_v5,
    is_v4 = is_v4,
    is_assay5 = is_assay5,
    mode = mode,
    data_access = data_access,
    default_assay = default_assay
  ))
}


#' Convert Seurat v5 to v4 Compatible Format (If Needed)
#'
#' @description
#' Internal helper to ensure Seurat object is in a format compatible with
#' downstream tools. For Seurat v5 objects with Assay5, this can optionally
#' join layers if needed.
#'
#' @param seurat_object A Seurat object
#' @param assay Assay name (default: "RNA")
#' @param join_layers Logical. For v5, join layers into single matrix? (default: FALSE)
#'
#' @return Seurat object (potentially modified)
#'
#' @keywords internal
prepareSeuratForCellChat <- function(seurat_object,
                                      assay = "RNA",
                                      join_layers = FALSE) {

  version_info <- checkSeuratVersion(seurat_object, verbose = FALSE)

  # If v5 with Assay5 and join_layers requested
  if (version_info$is_v5 && version_info$is_assay5 && join_layers) {

    # Check if JoinLayers function exists
    if (exists("JoinLayers", where = asNamespace("Seurat"))) {
      message("Joining Seurat v5 layers for compatibility...")
      tryCatch({
        seurat_object <- Seurat::JoinLayers(seurat_object, assay = assay)
      }, error = function(e) {
        warning("Failed to join layers: ", e$message,
                "\nProceeding with layered data...")
      })
    }
  }

  return(seurat_object)
}
