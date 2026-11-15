#' Plot Metabolism Signature Heatmap
#'
#' This function generates a heatmap of metabolism-related signature scores.
#'
#' When current human-reference metabolism signatures are applied to mouse samples, this
#' heatmap is an exploratory, method-specific description. It does not enter cross-species similarity scoring
#' or model-selection decisions and does not establish cross-species validation without
#' mouse-specific signatures or orthogonal metabolic evidence.
#'
#' @param sig_meta Data frame. Output from 'calculate_metabolism_score()'.
#' @param group_as Character vector. Group labels aligned to sample rows in 'sig_meta'.
#' @param output_path Character. File path to save the heatmap (default: "./metabolism/").
#' @param clustering_method Character. Clustering distance (`"pearson"`,
#'   `"manhattan"`, or `"canberra"`). The default `"pearson"` clusters groups
#'   by their pathway patterns rather than absolute score magnitude.
#' @param scale_rows Logical. If `TRUE` (default), standardize each pathway
#'   across the displayed groups before clustering and plotting.
#' @param width Numeric. Width of the output PDF (default: 10).
#' @param height Numeric. Height of the output PDF (default: 22).
#' @param row_name_width Numeric. Maximum width for row names in centimeters (default: 7).
#' @param right_padding Numeric. Right margin in centimeters to avoid label cutoff (default: 6).
#' @param row_height_pt Numeric. Height of each heatmap row in points (default: 5).
#'
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid gpar unit
#'
#' @return Saves the heatmap and returns the heatmap object.
#'
#' @seealso [plot_TME_barplot()], [plot_TME_boxplot()],
#'   [calculate_metabolism_score()]
#'
#' @export
plot_metabolism_heatmap <- function(sig_meta,
                                    group_as,
                                    output_path = "./metabolism",
                                    clustering_method = "pearson",
                                    scale_rows = TRUE,
                                    width = 10,
                                    height = 22,
                                    row_name_width = 7,
                                    right_padding = 6,
                                    row_height_pt = 5) {
  if (!is.data.frame(sig_meta)) stop("Error: 'sig_meta' must be a data frame.")
  if (!is.character(group_as) || length(group_as) != nrow(sig_meta)) {
    stop("Error: 'group_as' must be a character vector with the same length as the number of rows in 'sig_meta'.")
  }
  if (!clustering_method %in% c("pearson", "manhattan", "canberra")) {
    stop("Error: 'clustering_method' must be 'pearson', 'manhattan', or 'canberra'.")
  }
  if (!is.logical(scale_rows) || length(scale_rows) != 1L || is.na(scale_rows)) {
    stop("Error: 'scale_rows' must be TRUE or FALSE.")
  }

  score_data <- sig_meta[, -1, drop = FALSE]
  if (ncol(score_data) == 0L || !all(vapply(score_data, is.numeric, logical(1)))) {
    stop("Error: all metabolism score columns must be numeric.")
  }

  group_levels <- unique(group_as)
  group_matrix <- matrix(
    NA_real_,
    nrow = length(group_levels),
    ncol = ncol(score_data),
    dimnames = list(group_levels, colnames(score_data))
  )
  mean_or_na <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }
  for (i in seq_along(group_levels)) {
    group_scores <- score_data[group_as == group_levels[i], , drop = FALSE]
    group_matrix[i, ] <- vapply(group_scores, mean_or_na, numeric(1))
  }
  df <- t(group_matrix)

  if (scale_rows) {
    scale_one_row <- function(x) {
      finite <- is.finite(x)
      out <- rep(NA_real_, length(x))
      if (!any(finite)) return(out)
      row_sd <- stats::sd(x[finite])
      if (!is.finite(row_sd) || row_sd == 0) {
        out[finite] <- 0
      } else {
        out[finite] <- (x[finite] - mean(x[finite])) / row_sd
      }
      out
    }
    df <- t(apply(df, 1, scale_one_row))
    dimnames(df) <- list(colnames(score_data), group_levels)
  }

  cluster_axis <- function(x) {
    if (nrow(x) < 2L) return(FALSE)
    x_for_distance <- x
    x_for_distance[!is.finite(x_for_distance)] <- 0
    if (clustering_method == "pearson") {
      similarity <- stats::cor(t(x_for_distance))
      similarity[!is.finite(similarity)] <- 0
      diag(similarity) <- 1
      distance_matrix <- 1 - similarity
      distance_matrix[distance_matrix < 0] <- 0
      distance <- stats::as.dist(distance_matrix)
    } else {
      distance <- stats::dist(x_for_distance, method = clustering_method)
    }
    stats::hclust(distance)
  }

  finite_values <- abs(df[is.finite(df)])
  color_limit <- if (length(finite_values) == 0L || max(finite_values) == 0) {
    1
  } else {
    max(finite_values)
  }

  heatmap_obj <- Heatmap(
    df,
    name = if (scale_rows) "row Z-score" else "score",
    col = colorRamp2(
      c(-color_limit, 0, color_limit),
      c("#2166AC", "#F7F7F7", "#B2182B")
    ),
    na_col = "#D9D9D9",
    cluster_columns = cluster_axis(t(df)),
    cluster_rows = cluster_axis(df),
    row_names_gp = gpar(fontsize = 8),
    row_names_max_width = unit(row_name_width, "cm"),
    height = unit(nrow(df) * row_height_pt, "pt")
  )

  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  output_file <- file.path(output_path, "heatmap.pdf")

  grDevices::pdf(output_file, width = width, height = height, useDingbats = FALSE)
  draw(heatmap_obj,
    heatmap_legend_side = "right",
    padding = unit(c(2, 2, 2, right_padding), "cm")
  )
  grDevices::dev.off()

  message("Heatmap saved to: ", output_file)

  return(heatmap_obj)
}
