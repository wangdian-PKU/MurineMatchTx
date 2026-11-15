# Load the complete metabolism signature collection distributed with IOBR.
.get_metabolism_signatures <- function() {
  if ("load_data" %in% getNamespaceExports("IOBR")) {
    signatures <- suppressWarnings(
      getExportedValue("IOBR", "load_data")("signature_metabolism")
    )
  } else {
    data_env <- new.env(parent = emptyenv())
    utils::data("signature_metabolism", package = "IOBR", envir = data_env)
    signatures <- data_env$signature_metabolism
  }

  if (!is.list(signatures) || length(signatures) == 0L) {
    stop("Error: the IOBR metabolism signature collection could not be loaded.")
  }

  lapply(signatures, function(genes) unique(stats::na.omit(as.character(genes))))
}
