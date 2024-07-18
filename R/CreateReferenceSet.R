#' Signature analysis outputs into a reference set
#'
#' Point to a folder full of outputs from \strong{sigminerUtils::sig_analyse_mutations()} and create
#'
#' @return
#' @export
#'
#' @examples
sig_create_reference_set <- function(path_to_signature_directory = "signatures/"){
  dir(path_to_signature_directory)
}


#' Parse Signature Analysis Output Files
#'
#' Parse Signature Analysis Output Files
#'
#' @param x paths to any file produced by sig_analyse_mutations (character vector)
#'
#' @return a data.frame with 1 row per path, that includes filename-encoded metadata and a list-column with dataframes (contents)
#' @export
#'
parse_sig_files <- function(x){
  filenames=basename(x)
  ls_res = strsplit(basename(x), split = "\\.")

  res <- lapply(seq_along(ls_res), \(i) {
    v <- ls_res[[i]]
    sig_class = sub(x= v[1], "_.*", "")
    filename = filenames[i]
    filepath = x[i]

    list(
      "sigtype" = sig_class,
      "sampleId" = v[[2]],
      "refgenome" = v[[3]],
      "filetype" = v[[4]],
      "extension" = paste0(v[5:length(v)], collapse = "."),
      "filename" = filename,
      "filepath" = filepath,
      "contents" = read.csv(filepath, header = TRUE)
      )
    })

  df = as.data.frame(do.call("rbind", res))
  #names(res) <- vapply(res, FUN = \(v){v[names(v) == "sampleId"]}, FUN.VALUE = character(1))

  #names(res) <- sampleIds
  return(df)
}
