#' Write compressed CSV
#'
#' @inheritParams utils::write.csv
#'
write_compressed_csv <- function(x, file, row.names = FALSE){

  # Remove .gz extensions
  file <- sub(x=file, pattern = "\\.gz", replacement = "")

  utils::write.csv(x, file, row.names = row.names)
  #R.utils::gzip(file, ext="gz", overwrite = TRUE)
  R.utils::bzip2(file, ext="gz", overwrite = TRUE)
}
