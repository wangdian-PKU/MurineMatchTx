#' MurineMatchTx package
#'
#' Internal package level configuration for shared imports and check notes.
#'
#' @keywords internal
#' @importFrom utils read.csv read.table
"_PACKAGE"

utils::globalVariables(
  c(
    ".",
    "FDR",
    "ID",
    "PATIENT_ID",
    "change",
    "condition",
    "contrast",
    "foldchange",
    "group",
    "logFC",
    "sig",
    "tissue",
    "x",
    "x.y",
    "y"
  )
)
