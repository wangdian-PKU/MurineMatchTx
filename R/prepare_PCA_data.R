#' Prepare PCA Data for TCGA and Mouse Models
#'
#' Merge expression files by gene and join sample annotations by `sample_id`.
#' Sample labels are never inferred from column positions when
#' `sample_metadata` is supplied.
#'
#' @param tcga_file Path to the human expression file.
#' @param mouse_files Named list of mouse expression files.
#' @param sample_counts Deprecated-compatible named count list. It is used only
#'   when `sample_metadata` is `NULL`.
#' @param batch_correction Whether to apply `sva::ComBat`. Defaults to `FALSE`
#'   so PCA reflects the uncorrected expression structure.
#' @param sample_metadata Data frame containing at least `sample_id`, `species`,
#'   `dataset`, `model`, and `tissue`. Optional columns such as `batch`,
#'   `pair_id`, and `subtype` are retained.
#'
#' @return A list containing the uncorrected merged matrix, sample labels, batch
#'   labels, an optional corrected matrix, parameter settings, and metadata in
#'   matrix-column order.
#'
#' @importFrom sva ComBat
#' @export
prepare_PCA_data <- function(tcga_file,
                             mouse_files,
                             sample_counts = NULL,
                             batch_correction = FALSE,
                             sample_metadata = NULL) {
  read_expression_file <- function(file) {
    extension <- tolower(tools::file_ext(file))
    if (extension == "csv") {
      data <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (extension == "tsv") {
      data <- utils::read.table(
        file,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      stop("Unsupported expression file format: use .csv or .tsv.", call. = FALSE)
    }
    if (!"X" %in% colnames(data)) {
      stop("Expression files must contain an 'X' gene identifier column.", call. = FALSE)
    }
    if (anyNA(data$X) || any(!nzchar(as.character(data$X))) || anyDuplicated(data$X)) {
      stop("Expression gene identifiers must be non-missing and unique within each file.", call. = FALSE)
    }
    data
  }

  if (!is.list(mouse_files) || is.null(names(mouse_files)) || any(!nzchar(names(mouse_files)))) {
    stop("'mouse_files' must be a named list.", call. = FALSE)
  }
  tcga_data <- read_expression_file(tcga_file)
  mouse_list <- lapply(mouse_files, read_expression_file)

  all_sample_ids <- c(
    setdiff(colnames(tcga_data), "X"),
    unlist(lapply(mouse_list, function(data) setdiff(colnames(data), "X")), use.names = FALSE)
  )
  duplicate_samples <- unique(all_sample_ids[duplicated(all_sample_ids)])
  if (length(duplicate_samples) > 0) {
    stop(
      "Expression files have duplicate sample names: ",
      .format_sample_ids(duplicate_samples),
      call. = FALSE
    )
  }

  merged_data <- tcga_data
  for (mouse_data in mouse_list) {
    merged_data <- merge(merged_data, mouse_data, by = "X", sort = FALSE)
  }
  rownames(merged_data) <- merged_data$X
  merged_data <- merged_data[, setdiff(colnames(merged_data), "X"), drop = FALSE]

  if (is.null(sample_metadata)) {
    if (is.null(sample_counts)) {
      stop("Supply 'sample_metadata'; 'sample_counts' is only a deprecated compatibility input.", call. = FALSE)
    }
    warning(
      "'sample_counts' is deprecated; supply explicit 'sample_metadata' keyed by sample_id.",
      call. = FALSE
    )
    datasets <- c("TCGA", names(mouse_files))
    dataset_samples <- c(
      list(TCGA = setdiff(colnames(tcga_data), "X")),
      lapply(mouse_list, function(data) setdiff(colnames(data), "X"))
    )
    metadata_parts <- lapply(datasets, function(dataset) {
      counts <- sample_counts[[dataset]]
      if (is.null(counts) || !all(c("normal", "tumor") %in% names(counts))) {
        stop("Each legacy sample count must contain named 'normal' and 'tumor' values.", call. = FALSE)
      }
      if (sum(counts[c("normal", "tumor")]) != length(dataset_samples[[dataset]])) {
        stop("Legacy sample counts do not match expression samples for dataset '", dataset, "'.", call. = FALSE)
      }
      data.frame(
        sample_id = dataset_samples[[dataset]],
        species = if (dataset == "TCGA") "human" else "mouse",
        dataset = dataset,
        model = dataset,
        tissue = c(rep("normal", counts[["normal"]]), rep("tumor", counts[["tumor"]])),
        stringsAsFactors = FALSE
      )
    })
    sample_metadata <- do.call(rbind, metadata_parts)
  }

  sample_metadata <- .validate_sample_metadata(sample_metadata, colnames(merged_data))
  group_all <- stats::setNames(
    paste(sample_metadata$dataset, sample_metadata$tissue, sep = "_"),
    sample_metadata$sample_id
  )
  batch <- stats::setNames(sample_metadata$dataset, sample_metadata$sample_id)

  correction_requested <- isTRUE(batch_correction)
  if (correction_requested &&
      length(unique(sample_metadata$species)) > 1 &&
      all(vapply(split(sample_metadata$species, sample_metadata$dataset),
        function(species) length(unique(species)) == 1L,
        logical(1)
      ))) {
    warning(
      "Batch correction requested although species and dataset are completely confounded; ComBat may remove biological species differences.",
      call. = FALSE
    )
  }

  corrected_data <- if (correction_requested) {
    corrected_data <- sva::ComBat(dat = as.matrix(merged_data), batch = unname(batch))
    attr(corrected_data, "MurineMatchTx_tpm_source") <- "combat_corrected"
    corrected_data
  } else {
    NULL
  }

  list(
    merged_data = merged_data,
    group_all = group_all,
    batch = batch,
    uncorrected_data = merged_data,
    corrected_data = corrected_data,
    parameters = list(batch_correction = correction_requested),
    sample_metadata = sample_metadata
  )
}
