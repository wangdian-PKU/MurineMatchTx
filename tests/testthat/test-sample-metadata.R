make_expression_fixture <- function() {
  set.seed(42)
  sample_ids <- c(
    "TCGA-AA-0001-01", "TCGA-AA-0002-11",
    "mouse_normal", "mouse_tumor"
  )
  expression <- matrix(
    rnorm(200 * length(sample_ids)),
    nrow = 200,
    dimnames = list(paste0("gene", seq_len(200)), sample_ids)
  )
  metadata <- data.frame(
    sample_id = sample_ids,
    species = c("human", "human", "mouse", "mouse"),
    dataset = c("TCGA", "TCGA", "GSE-test", "GSE-test"),
    model = c("TCGA", "TCGA", "ModelA", "ModelA"),
    tissue = c("tumor", "normal", "normal", "tumor"),
    stringsAsFactors = FALSE
  )
  clinic <- data.frame(
    PATIENT_ID = c("TCGA.AA.0001", "TCGA.AA.0002"),
    stringsAsFactors = FALSE
  )
  list(expression = expression, metadata = metadata, clinic = clinic)
}

test_that("prepare_PCA_data joins shuffled metadata by sample_id", {
  tcga_file <- tempfile(fileext = ".csv")
  mouse_file <- tempfile(fileext = ".csv")

  write.csv(
    data.frame(X = c("g1", "g2"), h_tumor = c(4, 5), h_normal = c(1, 2)),
    tcga_file,
    row.names = FALSE
  )
  write.csv(
    data.frame(X = c("g1", "g2"), m_normal = c(2, 3), m_tumor = c(5, 6)),
    mouse_file,
    row.names = FALSE
  )

  metadata <- data.frame(
    sample_id = c("m_tumor", "h_normal", "m_normal", "h_tumor"),
    species = c("mouse", "human", "mouse", "human"),
    dataset = c("GSE-test", "TCGA", "GSE-test", "TCGA"),
    model = c("ModelA", "TCGA", "ModelA", "TCGA"),
    tissue = c("tumor", "normal", "normal", "tumor"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PCA_data(
    tcga_file,
    list(ModelA = mouse_file),
    sample_metadata = metadata,
    batch_correction = FALSE
  )

  expect_identical(result$sample_metadata$sample_id, colnames(result$merged_data))
  expect_identical(
    unname(result$group_all),
    paste(result$sample_metadata$dataset, result$sample_metadata$tissue, sep = "_")
  )
  expect_identical(unname(result$batch), result$sample_metadata$dataset)
})

test_that("sample metadata errors name missing, extra, and duplicate samples", {
  fixture <- make_expression_fixture()
  metadata <- fixture$metadata
  sample_ids <- colnames(fixture$expression)

  expect_error(
    MurineMatchTx:::.validate_sample_metadata(metadata[-1, ], sample_ids),
    "missing.*TCGA-AA-0001-01"
  )

  extra <- rbind(metadata, transform(metadata[1, ], sample_id = "extra_sample"))
  expect_error(
    MurineMatchTx:::.validate_sample_metadata(extra, sample_ids),
    "extra.*extra_sample"
  )

  duplicate <- rbind(metadata, metadata[1, ])
  expect_error(
    MurineMatchTx:::.validate_sample_metadata(duplicate, sample_ids),
    "duplicate.*TCGA-AA-0001-01"
  )
})

test_that("TROM results are invariant to expression and metadata order", {
  fixture <- make_expression_fixture()

  result_a <- compute_similarity_matrix(
    fixture$expression,
    sample_metadata = fixture$metadata,
    merge_clinic_CNV = fixture$clinic,
    cnv_genes = character(),
    clinical_vars = character()
  )

  shuffled_ids <- c("mouse_tumor", "TCGA-AA-0002-11", "mouse_normal", "TCGA-AA-0001-01")
  result_b <- compute_similarity_matrix(
    fixture$expression[, shuffled_ids],
    sample_metadata = fixture$metadata[c(4, 2, 3, 1), ],
    merge_clinic_CNV = fixture$clinic[2:1, , drop = FALSE],
    cnv_genes = character(),
    clinical_vars = character()
  )

  expect_equal(
    result_a$similarity_matrix,
    result_b$similarity_matrix[
      rownames(result_a$similarity_matrix),
      colnames(result_a$similarity_matrix),
      drop = FALSE
    ]
  )
  expect_true(is.na(result_a$filtered_similarity_matrix["TCGA-AA-0001-01", "mouse_normal"]))
  expect_true(is.na(result_a$filtered_similarity_matrix["TCGA-AA-0002-11", "mouse_tumor"]))
  expect_equal(
    sort(result_a$dm_trom_del_backgroud_clinic$sample_id),
    sort(fixture$metadata$sample_id[fixture$metadata$species == "human"])
  )
})

test_that("zero TCGA normals keeps tumor comparisons and omits empty normal branch", {
  fixture <- make_expression_fixture()
  keep <- c("TCGA-AA-0001-01", "TCGA-AA-0002-11", "mouse_normal", "mouse_tumor")
  expression <- fixture$expression[, keep, drop = FALSE]
  metadata <- fixture$metadata[match(keep, fixture$metadata$sample_id), ]
  metadata$tissue[metadata$species == "human"] <- "tumor"
  clinic <- fixture$clinic

  similarity <- compute_similarity_matrix(
    expression,
    sample_metadata = metadata,
    merge_clinic_CNV = clinic,
    cnv_genes = character(),
    clinical_vars = character()
  )

  expect_true(all(is.na(similarity$filtered_similarity_matrix[, "mouse_normal"])))
  expect_true(all(!is.na(similarity$filtered_similarity_matrix[, "mouse_tumor"])))

  output_file <- tempfile(fileext = ".pdf")
  heatmap <- plot_TROMscore_heatmap(
    dm_trom_del_backgroud_clinic = similarity$dm_trom_del_backgroud_clinic,
    trom_matrix = similarity$similarity_matrix,
    sample_metadata = metadata,
    clinical_vars = character(),
    cnv_genes = character(),
    output_path = output_file
  )

  expect_true(file.exists(output_file))
  expect_identical(colnames(heatmap$plot_matrix), "mouse_tumor")
  expect_identical(heatmap$merge_long$sample_id, "mouse_tumor")
  expect_identical(as.character(heatmap$merge_long$group), "tumor")
})
