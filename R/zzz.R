# Package startup messages and hooks

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "CellDiff v", utils::packageVersion("CellDiff"), " loaded."
  )
}
