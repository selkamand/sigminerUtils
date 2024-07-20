#' Signature analysis outputs into a reference set
#'
#' Point to a folder full of outputs from \strong{sigminerUtils::sig_analyse_mutations()} and create
#'
#' @return this function is run for its side effects. Invisibly returns NULL.
#' @export
#'
#' @examples
sig_create_reference_set <- function(path_to_signature_directory = "signatures/", outfolder = "reference", format=c("parquet", "csv.gz")){
  # ---- Assertions ---- #
  format <- rlang::arg_match(format)
  cli::cli_alert_info
  if(!file.exists(outfolder))  dir.create(outfolder, recursive	= TRUE)
  if(!file.exists(outfolder))  dir.create(outfolder, recursive	= TRUE)

  outfile_tally <- paste0(outfolder, "/refmatrix.tally.", format)
  outfile_exposures <- paste0(outfolder, "/refmatrix.exposures.", format)

  assertions::assert_file_does_not_exist(outfile_tally)
  assertions::assert_file_does_not_exist(outfile_exposures)

  # ---- Catalogue Reference Set ---- #
  cli::cli_alert_info("Reading Catalalogue Tallies")
  path_catalogues <- dir(path_to_signature_directory, pattern = "\\.tally", full.names = TRUE)
  df <- parse_sig_files(path_catalogues)

  df_expanded <- df |>
    dplyr::select(class, sample, contents) |>
    tidyr::unnest(contents)

  cli::cli_alert_info("Writing tally refmatrix to {outfile_tally}")
  if(format == "parquet"){
    df_expanded |>
      dplyr::group_by(class,sample) |>
      arrow::write_dataset(outfile_tally, format = "parquet")
  }
  else{
    df_expanded |>
      readr::write_csv(outfile_tally)
  }
  cli::cli_alert_success("Tally refmatrix succesfully created")

  # Clean up memory
  remove(df)
  remove(df_expanded)


  # ---- Exposures ---- #
  cli::cli_alert_info("Reading Signature Exposure Results")
  path_exposures <- dir(path_to_signature_directory, pattern = "\\.expo\\.", full.names = TRUE)

  if(length(path_exposures) == 0){
    cli::cli_alert_warning("Skipping exposure reference set creation (no expo files describing signature contributions were found)")
    return(invisible(NULL))
  }
  df_exposures_nested <- parse_sig_files(path_exposures)
  df_exposures <- df_exposures_nested |>
    dplyr::select(class, sample, contents) |>
    tidyr::unnest(contents)

  if(format == "parquet"){
    df_exposures |>
      dplyr::group_by(class,sample) |>
      arrow::write_dataset(outfile_exposures, format = "parquet")
  }
  else{
    df_exposures |>
      readr::write_csv(outfile_exposures)
  }

  return(invisible(NULL))
  #return(df_exposures_nested)
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

    tibble::tibble(
      "class" = sig_class,
      "sample" = v[[2]],
      "refgenome" = v[[3]],
      "filetype" = v[[4]],
      "extension" = paste0(v[5:length(v)], collapse = "."),
      "filename" = filename,
      "filepath" = filepath,
      "contents" = list(read.csv(filepath, header = TRUE))
      )

    })

  df = do.call("rbind", res)
  #names(res) <- vapply(res, FUN = \(v){v[names(v) == "sampleId"]}, FUN.VALUE = character(1))

  #names(res) <- sampleIds
  return(df)
}

compute_similarity_against_reference_set <- function(tally_file = "./signatures/SBS96_catalogue.TCGA-A2-A0T5-01.hg19.tally.csv.gz", ref_tallies = "./reference/refmatrix.tally.parquet"){
  assertions::assert_file_exists(tally_file)
  assertions::assert_directory_exists(ref_tallies)

  ref_tallies_parquet <- arrow::open_dataset(ref_tallies)

  df_tally_nested <- parse_sig_files(tally_file)
  target_class <- df_tally_nested$class
  target_sample <- df_tally_nested$sample
  target_tally <- df_tally_nested$contents[[1]]

  df_reference_tallies <- dplyr::filter(ref_tallies_parquet, class == target_class, sample != target_sample) |>
    dplyr::collect() |>
    dplyr::group_by(sample)

  keys = df_reference_tallies |>
    dplyr::group_keys() |>
    unlist()

  ls_reference_tallies <- df_reference_tallies |>
    dplyr::group_split(.keep = FALSE) |>
    setNames(keys)

  similarities <- purrr::map_dbl(
    ls_reference_tallies,
    \(refsig){
      sigstats::sig_cosine_similarity(signature1 = target_tally, signature2 = refsig)
    }
  )

  df_sim <- dplyr::tibble(
    sample = names(similarities),
    cosine_similarity = similarities
  )

  df_sim <- dplyr::arrange(df_sim, dplyr::desc(cosine_similarity))
  names(df_sim)[2] <- paste0("cosine_similarity_", target_class)
  return(df_sim)
}
