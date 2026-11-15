.legacy_heatmap_metadata <- function(trom_matrix,
                                     group_all,
                                     tcga_sample_size,
                                     tcga_normal_count,
                                     condition) {
  sample_ids <- c(rownames(trom_matrix), colnames(trom_matrix))
  if (is.null(group_all) || length(group_all) != length(sample_ids)) {
    stop("Legacy heatmap inputs do not match TROM matrix sample names.", call. = FALSE)
  }
  metadata <- .legacy_trom_metadata(
    matrix(0, nrow = 1, ncol = length(sample_ids), dimnames = list("gene", sample_ids)),
    group_all,
    tcga_sample_size,
    tcga_normal_count
  )
  if (!is.null(condition)) {
    condition_models <- rep(names(condition), unlist(condition, use.names = FALSE))
    mouse_rows <- metadata$species == "mouse"
    if (length(condition_models) != sum(mouse_rows)) {
      stop("Legacy 'condition' counts do not match mouse samples.", call. = FALSE)
    }
    metadata$model[mouse_rows] <- condition_models
  }
  metadata
}

#' Plot TROM Score Heatmap for Human and Mouse Samples
#'
#' Align scores and annotations by sample ID. If the human reference contains no
#' normal samples, the empty normal-normal branch and normal mouse columns are
#' omitted while tumor-tumor scores remain available.
#'
#' @param dm_trom_del_backgroud_clinic Annotated table returned by
#'   [compute_similarity_matrix()].
#' @param trom_matrix Raw TROM matrix returned by [compute_similarity_matrix()].
#' @param tcga_sample_size,tcga_normal_count,model_group_count,condition,group_all
#'   Deprecated-compatible positional inputs used only when `sample_metadata`
#'   is `NULL`.
#' @param clinical_vars,cnv_genes Annotation columns to display.
#' @param width,height PDF dimensions.
#' @param output_path Output PDF path.
#' @param cluster_rows,cluster_cols Whether to cluster heatmap rows or columns.
#' @param sample_metadata Explicit sample metadata shared with PCA and TROM.
#'
#' @return A list containing the heatmap, model-level long data, and the exact
#'   matrix displayed.
#'
#' @importFrom ComplexHeatmap Heatmap draw HeatmapAnnotation rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit
#' @export
plot_TROMscore_heatmap <- function(dm_trom_del_backgroud_clinic,
                                   trom_matrix,
                                   tcga_sample_size = NULL,
                                   tcga_normal_count = NULL,
                                   model_group_count = NULL,
                                   condition = NULL,
                                   clinical_vars = character(),
                                   cnv_genes = character(),
                                   group_all = NULL,
                                   width = 11.84,
                                   height = 7.11,
                                   output_path = "./similarity/TROM_heatmap.pdf",
                                   cluster_rows = FALSE,
                                   cluster_cols = FALSE,
                                   sample_metadata = NULL) {
  if (!is.data.frame(dm_trom_del_backgroud_clinic)) {
    stop("'dm_trom_del_backgroud_clinic' must be a data frame.", call. = FALSE)
  }
  trom_matrix <- as.matrix(trom_matrix)
  if (!is.numeric(trom_matrix) || is.null(rownames(trom_matrix)) || is.null(colnames(trom_matrix))) {
    stop("'trom_matrix' must be a numeric matrix with human row and mouse column names.", call. = FALSE)
  }
  if (nrow(trom_matrix) == 0 || ncol(trom_matrix) == 0) {
    stop("'trom_matrix' must contain human and mouse samples.", call. = FALSE)
  }

  if (is.null(sample_metadata)) {
    warning(
      "Count- and position-based heatmap inputs are deprecated; supply explicit 'sample_metadata'.",
      call. = FALSE
    )
    sample_metadata <- .legacy_heatmap_metadata(
      trom_matrix,
      group_all,
      tcga_sample_size,
      tcga_normal_count,
      condition
    )
  }
  sample_metadata <- .validate_sample_metadata(
    sample_metadata,
    c(rownames(trom_matrix), colnames(trom_matrix))
  )
  human_metadata <- sample_metadata[match(rownames(trom_matrix), sample_metadata$sample_id), , drop = FALSE]
  mouse_metadata <- sample_metadata[match(colnames(trom_matrix), sample_metadata$sample_id), , drop = FALSE]
  if (any(human_metadata$species != "human") || any(mouse_metadata$species != "mouse")) {
    stop("TROM matrix rows must be human samples and columns must be mouse samples.", call. = FALSE)
  }

  if ("sample_id" %in% colnames(dm_trom_del_backgroud_clinic)) {
    displayed_human_ids <- as.character(dm_trom_del_backgroud_clinic$sample_id)
  } else {
    displayed_human_ids <- rownames(dm_trom_del_backgroud_clinic)
  }
  if (is.null(displayed_human_ids) || anyDuplicated(displayed_human_ids)) {
    stop("Annotated TROM data must contain unique human sample IDs.", call. = FALSE)
  }
  missing_human <- setdiff(rownames(trom_matrix), displayed_human_ids)
  extra_human <- setdiff(displayed_human_ids, rownames(trom_matrix))
  if (length(missing_human) > 0 || length(extra_human) > 0) {
    stop("Annotated TROM sample IDs do not match TROM matrix rows.", call. = FALSE)
  }
  human_metadata <- human_metadata[match(displayed_human_ids, human_metadata$sample_id), , drop = FALSE]

  available_tissues <- unique(human_metadata$tissue)
  compatible_mouse <- mouse_metadata$tissue %in% available_tissues
  mouse_metadata <- mouse_metadata[compatible_mouse, , drop = FALSE]
  if (nrow(mouse_metadata) == 0) {
    stop("No mouse samples have a tissue represented in the human reference.", call. = FALSE)
  }

  plot_matrix <- trom_matrix[
    human_metadata$sample_id,
    mouse_metadata$sample_id,
    drop = FALSE
  ]
  same_tissue <- outer(human_metadata$tissue, mouse_metadata$tissue, FUN = "==")
  plot_matrix[!same_tissue] <- NA_real_

  long_parts <- lapply(seq_len(nrow(mouse_metadata)), function(index) {
    mouse_id <- mouse_metadata$sample_id[index]
    tissue <- mouse_metadata$tissue[index]
    matching_human <- human_metadata$tissue == tissue
    data.frame(
      sample_id = mouse_id,
      x = paste(mouse_metadata$model[index], tissue, sep = "_"),
      y = mean(trom_matrix[human_metadata$sample_id[matching_human], mouse_id, drop = TRUE]),
      group = tissue,
      condition = mouse_metadata$model[index],
      stringsAsFactors = FALSE
    )
  })
  merge_long <- do.call(rbind, long_parts)
  group_means <- stats::aggregate(merge_long$y, by = list(merge_long$x), FUN = mean)
  merge_long$score <- round(group_means$x[match(merge_long$x, group_means$Group.1)], 1)
  merge_long$x <- factor(merge_long$x, levels = unique(merge_long$x))
  merge_long$condition <- factor(merge_long$condition, levels = unique(merge_long$condition))
  merge_long$group <- factor(merge_long$group, levels = intersect(c("normal", "tumor"), unique(merge_long$group)))

  color_function <- circlize::colorRamp2(
    c(0, 1, 5, 10, 20, 30),
    c("white", "#3f84a2", "#cae2a7", "#f8e4a1", "#d85b39", "#cb3025")
  )
  column_annotation <- ComplexHeatmap::HeatmapAnnotation(
    model_tissue = mouse_metadata$tissue,
    col = list(model_tissue = c(normal = "#79c2f2", tumor = "#e5616e"))
  )

  annotation_data <- data.frame(
    TCGA_tissue = human_metadata$tissue,
    stringsAsFactors = FALSE
  )
  color_mapping <- list(TCGA_tissue = c(normal = "#79c2f2", tumor = "#e5616e"))
  annotation_index <- match(human_metadata$sample_id, displayed_human_ids)
  for (variable in c(clinical_vars, cnv_genes)) {
    if (!variable %in% colnames(dm_trom_del_backgroud_clinic)) {
      stop("Annotated TROM data are missing column '", variable, "'.", call. = FALSE)
    }
    annotation_data[[variable]] <- dm_trom_del_backgroud_clinic[[variable]][annotation_index]
    if (variable %in% clinical_vars) {
      color_mapping[[variable]] <- c(Positive = "#f6d593", Negative = "#ddeaac")
    } else {
      color_mapping[[variable]] <- c("-2" = "#1010ff", "-1" = "#109eff", "0" = "#bababa", "1" = "#e274dd")
    }
  }
  row_annotation <- ComplexHeatmap::rowAnnotation(df = annotation_data, col = color_mapping)

  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(output_path, height = height, width = width)
  on.exit(grDevices::dev.off(), add = TRUE)
  heatmap_object <- ComplexHeatmap::Heatmap(
    plot_matrix,
    col = color_function,
    na_col = "white",
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    row_names_side = "left",
    heatmap_legend_param = list(title = "Score", legend_direction = "vertical"),
    row_names_gp = grid::gpar(fontsize = 10, font = 3),
    column_names_gp = grid::gpar(fontsize = 10, font = 3),
    row_split = factor(human_metadata$tissue, levels = intersect(c("tumor", "normal"), unique(human_metadata$tissue))),
    column_split = factor(mouse_metadata$model, levels = unique(mouse_metadata$model)),
    border = TRUE,
    row_gap = grid::unit(0, "mm"),
    column_gap = grid::unit(0, "mm"),
    show_column_names = FALSE,
    top_annotation = column_annotation,
    left_annotation = row_annotation
  )
  ComplexHeatmap::draw(heatmap_object)
  message("TROMscore_heatmap saved to: ", output_path)

  list(
    heatmap_obj = heatmap_object,
    merge_long = merge_long,
    plot_matrix = plot_matrix
  )
}
