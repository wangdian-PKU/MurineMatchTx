test_that("calculate_metabolism_score returns the complete IOBR metabolism panel", {
  eset_tpm <- data.frame(
    s1 = c(1, 2, 1),
    s2 = c(2, 4, 2),
    s3 = c(4, 8, 4),
    row.names = c("HK1", "LDHA", "TBXAS1")
  )

  sig_meta <- calculate_metabolism_score(
    eset_tpm,
    method = "relative_zscore",
    mini_gene_count = 1,
    input_source = "provided_tpm"
  )

  expect_s3_class(sig_meta, "data.frame")
  expect_identical(colnames(sig_meta)[1], "ID")
  expect_length(setdiff(colnames(sig_meta), "ID"), 114L)
  expect_true("Glycolysis" %in% colnames(sig_meta))
  expect_true("Thromboxane_Biosynthesis" %in% colnames(sig_meta))
  expect_false("T cells CD8" %in% colnames(sig_meta))
  expect_false(anyNA(sig_meta$Thromboxane_Biosynthesis))
  expect_equal(mean(sig_meta$Thromboxane_Biosynthesis), 0, tolerance = 1e-12)

  gene_counts <- attr(sig_meta, "MurineMatchTx_pathway_gene_counts")
  expect_identical(unname(gene_counts["Thromboxane_Biosynthesis"]), 1L)
})

test_that("metabolism signature loading suppresses connectivity probe warnings", {
  loader_env <- new.env(parent = environment(MurineMatchTx:::.get_metabolism_signatures))
  loader_env$getNamespaceExports <- function(package) "load_data"
  loader_env$getExportedValue <- function(package, name) {
    function(data_name) {
      warning("URL 'https://www.baidu.com/': Timeout of 5 seconds was reached", call. = FALSE)
      list(Glycolysis = c("HK1", "HK1", NA_character_))
    }
  }
  get_signatures <- MurineMatchTx:::.get_metabolism_signatures
  environment(get_signatures) <- loader_env

  expect_no_warning(get_signatures())
})

test_that("plot_metabolism_heatmap requires explicit grouping and writes output", {
  sig_meta <- data.frame(
    ID = c("s1", "s2", "s3", "s4"),
    Glycolysis = c(1, 2, 3, 4),
    OXPHOS = c(4, 3, 2, 1)
  )

  expect_error(
    plot_metabolism_heatmap(sig_meta, group_as = c("A", "A", "B"), output_path = tempdir()),
    "same length"
  )

  output_dir <- tempfile(pattern = "metabolism-plot-")
  heatmap_obj <- plot_metabolism_heatmap(
    sig_meta,
    group_as = c("A", "A", "B", "B"),
    output_path = output_dir,
    height = 6
  )

  expect_true(file.exists(file.path(output_dir, "heatmap.pdf")))
  expect_s4_class(heatmap_obj, "Heatmap")
  expect_equal(rowMeans(heatmap_obj@matrix), c(Glycolysis = 0, OXPHOS = 0), tolerance = 1e-12)
  expect_equal(apply(heatmap_obj@matrix, 1, stats::sd), c(Glycolysis = 1, OXPHOS = 1), tolerance = 1e-12)
})

test_that("GO symbol map lookup works for supported species", {
  human_map <- MurineMatchTx:::.get_symbol_map("hsa")
  mouse_map <- MurineMatchTx:::.get_symbol_map("mmu")

  expect_s3_class(human_map, "data.frame")
  expect_s3_class(mouse_map, "data.frame")
  expect_true(all(c("gene_id", "symbol") %in% colnames(human_map)))
  expect_true(all(c("gene_id", "symbol") %in% colnames(mouse_map)))
})
