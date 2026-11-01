# Run from the repository root: Rscript --vanilla Code/get_raw_data.R
# Raw data: https://doi.org/10.5281/zenodo.23857927
# Base R only; restore the historical filenames used by the analysis scripts.
get_raw_data <- function(data_dir = "Data") {
  if (!file.exists("Spatial-Multiomics.Rproj"))
    stop("Run this script from the repository root.")
  files <- data.frame(
    remote = c("A26-Q500 Data.xlsx", "A27-Q500 Data.xlsx",
               "5583 FIA Q500 Data.xlsx", "5583 UPLC Q500 Data.xlsx",
               "A26_Flatmap_Grid_20260315.xlsx", "A27_Flatmap_Grid_20260317.xlsx"),
    local = c("A26-Q500 Data.xlsx", "A27-Q500 Data.xlsx",
              "5583 FIA Q500 Data.xlsx", "5583 UPLC Q500 Data.xlsx",
              "A26_Bonnie_Flatmap_Grid_20260315.xlsx",
              "A27_Frizzle_Flatmap_Grid_20260317.xlsx"),
    md5 = c("dd7b693f29c252c50ac348c278ac22c1", "a7c135f4d13cef4c2a544b3240624721",
            "e3a6417954f15432deae312e0baa0401", "a7d5c2e003df225220522853541e64ef",
            "a6235a08737e52774fff5e72b4ffa50e", "b918c58b3661834ea12bd9d05e56e1a6")
  )
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  old_timeout <- getOption("timeout")
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  for (i in seq_len(nrow(files))) {
    dest <- file.path(data_dir, files$local[i])
    if (file.exists(dest)) {
      if (unname(tools::md5sum(dest)) != files$md5[i])
        stop("Existing file differs from Zenodo; move it aside before retrying: ", dest)
      message("Already verified: ", dest)
      next
    }
    tmp <- tempfile(tmpdir = data_dir, fileext = ".xlsx")
    tryCatch({
      url <- paste0("https://zenodo.org/records/23857927/files/",
                    utils::URLencode(files$remote[i], reserved = TRUE), "?download=1")
      utils::download.file(url, tmp, mode = "wb", method = "libcurl")
      if (unname(tools::md5sum(tmp)) != files$md5[i])
        stop("Zenodo checksum mismatch: ", files$remote[i])
      if (!file.rename(tmp, dest)) stop("Could not install: ", dest)
      message("Downloaded and verified: ", dest)
    }, finally = unlink(tmp))
  }
  invisible(file.path(data_dir, files$local))
}
get_raw_data()
