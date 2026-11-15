.format_sample_ids <- function(x) {
  paste(utils::head(x, 5), collapse = ", ")
}

.validate_sample_metadata <- function(sample_metadata, sample_ids = NULL) {
  required <- c("sample_id", "species", "dataset", "model", "tissue")
  if (!is.data.frame(sample_metadata)) {
    stop("'sample_metadata' must be a data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(required, colnames(sample_metadata))
  if (length(missing_columns) > 0) {
    stop(
      "'sample_metadata' is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  sample_metadata <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)
  for (column in required) {
    sample_metadata[[column]] <- as.character(sample_metadata[[column]])
    if (anyNA(sample_metadata[[column]]) || any(!nzchar(sample_metadata[[column]]))) {
      stop("'sample_metadata$", column, "' contains missing or empty values.", call. = FALSE)
    }
  }

  duplicate_ids <- unique(sample_metadata$sample_id[duplicated(sample_metadata$sample_id)])
  if (length(duplicate_ids) > 0) {
    stop(
      "'sample_metadata' has duplicate sample_id values: ",
      .format_sample_ids(duplicate_ids),
      call. = FALSE
    )
  }

  sample_metadata$species <- tolower(sample_metadata$species)
  invalid_species <- setdiff(unique(sample_metadata$species), c("human", "mouse"))
  if (length(invalid_species) > 0) {
    stop("'species' must contain only 'human' or 'mouse'.", call. = FALSE)
  }

  sample_metadata$tissue <- tolower(sample_metadata$tissue)
  invalid_tissue <- setdiff(unique(sample_metadata$tissue), c("normal", "tumor"))
  if (length(invalid_tissue) > 0) {
    stop("'tissue' must contain only 'normal' or 'tumor'.", call. = FALSE)
  }

  if (!is.null(sample_ids)) {
    sample_ids <- as.character(sample_ids)
    if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
      stop("Expression data contain missing or empty sample names.", call. = FALSE)
    }
    duplicate_expression_ids <- unique(sample_ids[duplicated(sample_ids)])
    if (length(duplicate_expression_ids) > 0) {
      stop(
        "Expression data have duplicate sample names: ",
        .format_sample_ids(duplicate_expression_ids),
        call. = FALSE
      )
    }

    missing_ids <- setdiff(sample_ids, sample_metadata$sample_id)
    extra_ids <- setdiff(sample_metadata$sample_id, sample_ids)
    if (length(missing_ids) > 0) {
      stop(
        "'sample_metadata' is missing expression samples: ",
        .format_sample_ids(missing_ids),
        call. = FALSE
      )
    }
    if (length(extra_ids) > 0) {
      stop(
        "'sample_metadata' has extra samples not present in expression data: ",
        .format_sample_ids(extra_ids),
        call. = FALSE
      )
    }
    sample_metadata <- sample_metadata[match(sample_ids, sample_metadata$sample_id), , drop = FALSE]
  }

  rownames(sample_metadata) <- NULL
  sample_metadata
}
