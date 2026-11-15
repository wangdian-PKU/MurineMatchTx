test_that("prepare_DEG_heatmap_data applies custom TCGA thresholds", {
  tcga_file <- tempfile(fileext = ".tsv")
  mouse_a_file <- tempfile(fileext = ".tsv")
  mouse_b_file <- tempfile(fileext = ".tsv")

  tcga_data <- data.frame(
    logFC = c(0.8, 1.2),
    AveExpr = c(10, 12),
    PValue = c(0.001, 0.001),
    FDR = c(0.01, 0.01),
    row.names = c("TP53", "EGFR")
  )

  mouse_data <- data.frame(
    logFC = c(1.1, 1.3),
    AveExpr = c(8, 9),
    PValue = c(0.001, 0.001),
    FDR = c(0.01, 0.01),
    row.names = c("TP53", "EGFR")
  )

  write.table(tcga_data, tcga_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(mouse_data, mouse_a_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(mouse_data, mouse_b_file, sep = "\t", quote = FALSE, col.names = NA)

  mouse_files <- list(ModelA = mouse_a_file, ModelB = mouse_b_file)

  ortholog_map <- data.frame(
    human_gene_id = c("TP53", "EGFR"),
    mouse_gene_id = c("TP53", "EGFR"),
    human_gene_id_version = c("TP53", "EGFR"),
    mouse_gene_id_version = c("TP53", "EGFR"),
    homology_type = "ortholog_one2one",
    is_high_confidence = TRUE,
    release = 116,
    stringsAsFactors = FALSE
  )

  default_result <- prepare_DEG_heatmap_data(tcga_file, mouse_files, ortholog_map = ortholog_map)
  custom_result <- prepare_DEG_heatmap_data(
    tcga_file,
    mouse_files,
    logFC = 0.5,
    FDR = 0.05,
    ortholog_map = ortholog_map
  )

  expect_equal(nrow(default_result$processed_tcga), 1)
  expect_equal(nrow(custom_result$processed_tcga), 2)
})

test_that("prepare_DEG_heatmap_data validates threshold inputs", {
  tcga_file <- tempfile(fileext = ".tsv")
  mouse_file <- tempfile(fileext = ".tsv")

  deg_data <- data.frame(
    logFC = c(1.2, -1.1),
    AveExpr = c(10, 12),
    PValue = c(0.001, 0.001),
    FDR = c(0.01, 0.02),
    row.names = c("TP53", "EGFR")
  )

  write.table(deg_data, tcga_file, sep = "\t", quote = FALSE, col.names = NA)
  write.table(deg_data, mouse_file, sep = "\t", quote = FALSE, col.names = NA)

  mouse_files <- list(ModelA = mouse_file)

  ortholog_map <- data.frame(
    human_gene_id = c("TP53", "EGFR"),
    mouse_gene_id = c("TP53", "EGFR"),
    human_gene_id_version = c("TP53", "EGFR"),
    mouse_gene_id_version = c("TP53", "EGFR"),
    homology_type = "ortholog_one2one",
    is_high_confidence = TRUE,
    release = 116,
    stringsAsFactors = FALSE
  )

  expect_error(
    prepare_DEG_heatmap_data(tcga_file, mouse_files, logFC = -1, ortholog_map = ortholog_map),
    "non-negative"
  )

  expect_error(
    prepare_DEG_heatmap_data(tcga_file, mouse_files, FDR = 2, ortholog_map = ortholog_map),
    "between 0 and 1"
  )
})
