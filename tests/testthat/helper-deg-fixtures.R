write_deg_fixture <- function(data, extension = ".tsv", row_names = TRUE) {
  path <- tempfile(fileext = extension)
  if (extension == ".csv") {
    write.csv(data, path, row.names = row_names)
  } else {
    write.table(
      data,
      path,
      sep = "\t",
      quote = FALSE,
      row.names = row_names,
      col.names = if (row_names) NA else TRUE
    )
  }
  path
}
