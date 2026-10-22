args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
stage <- match.arg(args, c("install", "test", "vignette", "build", "check", "environment"))
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
results <- file.path(dirname(root), "check-results")
dir.create(results, showWarnings = FALSE)
log <- file(file.path(results, paste0(stage, ".log")), open = "wt")
sink(log, split = TRUE)
sink(log, type = "message")
set.seed(42)
cat("Stage:", stage, "\nUTC:", format(Sys.time(), tz = "UTC"), "\n")
cat(paste("Rscript .github/scripts/ci.R", stage), "\n",
    file = file.path(results, "commands.log"), append = TRUE)

library_dir <- file.path(results, "library")
dir.create(library_dir, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
description <- read.dcf(file.path(root, "DESCRIPTION"))
package <- description[1, "Package"]
tarball <- file.path(results, paste0(package, "_", description[1, "Version"], ".tar.gz"))

run_r <- function(arguments, label) {
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
  cat(paste(c(shQuote(executable), shQuote(arguments)), collapse = " "), "\n",
      file = file.path(results, "commands.log"), append = TRUE)
  output <- file.path(results, paste0(label, "-command.log"))
  status <- system2(executable, shQuote(arguments), stdout = output, stderr = output)
  cat(readLines(output, warn = FALSE), sep = "\n")
  if (status != 0L) stop(label, " failed with exit status ", status)
}

tryCatch({
if (stage == "install") {
  run_r(c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", library_dir), root), "install")
  stopifnot(normalizePath(find.package(package), winslash = "/") == paste0(library_dir, "/", package))
} else if (stage == "test") {
  tests <- testthat::test_local(root, reporter = "summary", stop_on_failure = TRUE)
  counts <- colSums(as.data.frame(tests)[, c("failed", "skipped", "error", "warning", "passed")])
  print(counts)
  stopifnot(counts[["failed"]] == 0, counts[["error"]] == 0,
            counts[["warning"]] == 0, counts[["skipped"]] == 0)
} else if (stage == "vignette") {
  stopifnot(normalizePath(find.package(package), winslash = "/") == paste0(library_dir, "/", package))
  rmarkdown::render(file.path(root, "vignettes", "MurineMatchTx_Workflow.Rmd"),
                    output_dir = results, intermediates_dir = results,
                    envir = new.env(parent = globalenv()), quiet = FALSE)
} else if (stage == "build") {
  stopifnot(!file.exists(tarball))
  setwd(results)
  run_r(c("CMD", "build", root), "build")
  stopifnot(file.exists(tarball))
  writeLines(c(basename(tarball), digest::digest(file = tarball, algo = "sha256")),
             file.path(results, "tarball-sha256.txt"))
} else if (stage == "check") {
  stopifnot(file.exists(tarball), !dir.exists(file.path(results, paste0(package, ".Rcheck"))))
  hash <- readLines(file.path(results, "tarball-sha256.txt"))
  stopifnot(identical(hash, c(basename(tarball), digest::digest(file = tarball, algo = "sha256"))))
  setwd(results)
  run_r(c("CMD", "check", "--as-cran", tarball), "check")
  check_log <- readLines(file.path(results, paste0(package, ".Rcheck"), "00check.log"))
  stopifnot(!any(grepl(" (ERROR|WARNING)$", check_log)))
} else if (stage == "environment") {
  info <- capture.output(sessionInfo())
  writeLines(info, file.path(results, "sessionInfo.txt"))
  cat(info, sep = "\n")
  system_info <- capture.output(print(Sys.info()))
  writeLines(system_info, file.path(results, "system-info.txt"))
  cat(system_info, sep = "\n")
  installed <- installed.packages(fields = c("RemoteType", "RemoteSha"))
  write.csv(installed[, c("Package", "Version", "Built", "RemoteType", "RemoteSha")],
            file.path(results, "installed-packages.csv"), row.names = FALSE)
  file.copy(file.path(root, ".github/workflows/R-CMD-check.yaml"), results, overwrite = TRUE)
  file.copy(file.path(root, ".github/scripts/ci.R"), results, overwrite = TRUE)
}
cat("Stage completed:", stage, "\n")
}, finally = {
if (stage != "environment") {
  info <- capture.output(sessionInfo())
  writeLines(info, file.path(results, paste0("sessionInfo-", stage, ".txt")))
  cat(info, sep = "\n")
}
sink(type = "message")
sink()
close(log)
})
