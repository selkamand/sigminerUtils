#' Signature analysis outputs into a reference set
#'
#' Point to a folder full of outputs from \strong{sigminerUtils::sig_analyse_mutations()} and create
#'
#' @return this function is run for its side effects. Invisibly returns NULL.
#' @param path_to_signature_directory path to signature directory
#' @param outfolder path to save reference set
#' @param format format to save reference set as (parquet or bgzipped csv)
#' @param umap_n_neighbours This parameter sets the balance between local and global structure in UMAP.
#' The value determines the size of the local neighborhood UMAP which is used to learn the manifold structure of the data.
#' Lower values focus on local details. higher values emphasise more global patterns. Must be less than the number of samples in your dataset.
#' @param create_umaps should umap references be created (flag)
#' @param seed random seed used for umap generation
#' @param verbose verbose (flag)
#' @export
#'
sig_create_reference_set <- function(
    path_to_signature_directory = "signatures/",
    outfolder = "reference",
    format=c("parquet", "csv.gz"),
    create_umaps = TRUE,
    umap_n_neighbours=15,
    seed=111,
    verbose = TRUE
    ){
  # ---- Assertions ---- #
  v <- verbose
  format <- rlang::arg_match(format)
  cli::cli_alert_info
  if(!file.exists(outfolder))  dir.create(outfolder, recursive	= TRUE)
  if(!file.exists(outfolder))  dir.create(outfolder, recursive	= TRUE)
  if(create_umaps) rlang::check_installed("uwot", reason = "Creating UMAPs in reference set requires the 'uwot' pacakge to be installed. Please run install.packages('uwot') and restart R.")

  # Convert outfolder to absolute filepath
  outfolder = normalizePath(outfolder, mustWork = TRUE)
  outfile_tally <- paste0(outfolder, "/refmatrix.tally.", format)
  outfile_exposures <- paste0(outfolder, "/refmatrix.exposures.", format)
  outfile_bootstraps <- paste0(outfolder, "/refmatrix.bootstraps.", format)
  outfile_umaps <- paste0(outfolder, "/refmatrix.umaps")

  assertions::assert_file_does_not_exist(outfile_tally)
  assertions::assert_file_does_not_exist(outfile_exposures)
  assertions::assert_file_does_not_exist(outfile_bootstraps)
  assertions::assert_directory_exists(path_to_signature_directory)


  # Get data.frame describing every catalogue in the reference directory (based on files with .tally extension)
  df_catalogues <- get_catalogue_dataframe(path_to_signature_directory)


  # ---- Catalogue Reference Set ---- #
  if(v) cli::cli_h1("Catalogue/Tally Reference Matrix")
  build_catalogue_reference_set(df_catalogues, outfile_tally, format=format, verbose=verbose)


  # Create UMAP reference ---------------------------------------------------
  if(v) cli::cli_h1("UMAPs (from sample catalogues)")
  if (create_umaps){
    build_umap_reference_set(df_catalogues, outfile_umap_prefix = outfile_umaps, umap_n_neighbours = umap_n_neighbours, seed=seed, verbose=verbose)
  }else
    cli::cli_alert_warning("Skipping UMAP creation since {.arg create_umaps = FALSE}")


  # ---- Exposures ---- #
  df_exposures <- get_exposures_dataframe(path_to_signature_directory)
  if(v) cli::cli_h1("Exposure Reference Matrix")
  build_exposures_reference_set(df_exposures, outfile_exposures, format=format, verbose=verbose)

  # ---- Bootstraps ---- #
  df_bootstraps <- get_boostrap_summary_dataframe(path_to_signature_directory)
  build_bootstrap_reference_set(df_bootstraps, outfile_bootstraps, format=format, verbose=verbose)

  # ---- Instructions: Using Outputs ---- #
  if (v) {
    cli::cli_h1("Usage Instructions")
    cli::cli_text("")

    cli::cli_alert_success(
      "Reference files were all successfully created and can be used in {.code sig_analyse_mutations()} and related functions."
    )
    cli::cli_text("")

    if (create_umaps) {
      cli::cli_alert_info("Project a sample onto the reference umaps include the argument:")
      cli::cli_text("")
      cli::cli_code(glue::glue("ref_umap_prefix = '{outfile_umaps}'"))
      cli::cli_text("")
    }

    cli::cli_alert_info("To compare a sample to this reference set include the argument:")
    cli::cli_text("")
    cli::cli_code(glue::glue("ref_tallies = '{outfile_tally}'"))
    cli::cli_text("")
  }

  return(invisible(NULL))
}


#' Parse Signature Analysis Output Files
#'
#' Parse Signature Analysis Output Files
#'
#' @param x paths to any file produced by sig_analyse_mutations (character vector)
#' @inheritParams  utils::read.csv
#' @return a data.frame with 1 row per path, that includes filename-encoded metadata and a list-column with dataframes (contents)
#' @export
#'
parse_sig_files <- function(x, colClasses = NA){
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
      "contents" = list(utils::read.csv(filepath, header = TRUE, colClasses = colClasses))
      )

    })

  df = do.call("rbind", res)
  #names(res) <- vapply(res, FUN = \(v){v[names(v) == "sampleId"]}, FUN.VALUE = character(1))

  #names(res) <- sampleIds
  return(df)
}

get_catalogue_dataframe <- function(path_to_signature_directory){
  path_catalogues <- dir(path_to_signature_directory, pattern = "\\.tally", full.names = TRUE, recursive = TRUE)

  # Nested dataframe with 1 row per catalogue file (including a list-column with the full catalogue)
  df <- parse_sig_files(path_catalogues)

  df_catalogues <- df |>
    dplyr::select(class, sample, contents) |>
    tidyr::unnest(contents) |>
    # Replace fraction NAs with zeros when total count sums to 0
    # (so that when we pivot to wide format for UMAP creation, any NAs introduced are because that sample wasn't run on that
    # class of signature analysis)
    dplyr::mutate(fraction = ifelse(is.na(fraction), yes = 0, no = fraction))

  return(df_catalogues)
}

get_exposures_dataframe <- function(path_to_signature_directory){
  path_exposures <- dir(path_to_signature_directory, pattern = "\\.expo\\.", full.names = TRUE, recursive = TRUE)

  if(length(path_exposures) == 0){
    cli::cli_alert_warning("Skipping exposure reference set creation (no expo files describing signature contributions were found)")
    return(invisible(NULL))
  }
  df_exposures_nested <- parse_sig_files(path_exposures)
  df_exposures <- df_exposures_nested |>
    dplyr::select(class, sample, contents) |>
    tidyr::unnest(contents)

  return(df_exposures)
}

get_boostrap_summary_dataframe <- function(path_to_signature_directory){
  path <- dir(path_to_signature_directory, pattern = "\\.bootstrap_summary\\.", full.names = TRUE, recursive = TRUE)

  if(length(path) == 0){
    cli::cli_alert_warning("Skipping boostrap summary reference set creation (no bootstrap summary files were found)")
    return(invisible(NULL))
  }

  df_bootstraps_nested <- parse_sig_files(path, colClasses = c(outliers = "character"))
  df_bootstraps <- df_bootstraps_nested |>
    dplyr::select(class, sample, contents) |>
    tidyr::unnest(contents)

  return(df_bootstraps)
}

build_catalogue_reference_set <- function(df_catalogues, outfile_tally, format, verbose=TRUE){
  v <- verbose
  if(v) cli::cli_alert_info("Reading Catalalogue Tallies")
  if(v) cli::cli_alert_info("Writing tally refmatrix to {outfile_tally}")

  if(format == "parquet"){
    df_catalogues |>
      dplyr::group_by(class,sample) |>
      arrow::write_dataset(outfile_tally, format = "parquet")
  }
  else{
    df_catalogues |>
      readr::write_csv(outfile_tally)
  }
  if(v) cli::cli_alert_success("Tally refmatrix succesfully created")
}

build_umap_reference_set <- function(df_catalogues, seed=111, outfile_umap_prefix, umap_n_neighbours, verbose=TRUE){
  v <- verbose

  # Lets turn each feature in our tally into a a column and each sample into a row (values are fraction) - then we can create
  # a UMAP (1 for each sigclass and 1 with all combined)
  if(v) cli::cli_progress_step("Building UMAPs from sample catalogues")
  df_umap_tbl <- df_catalogues |>
    tidyr::pivot_wider(id_cols = c("sample"), names_from = c(class, channel), names_sep = "|", values_from = fraction) |>
    tibble::column_to_rownames(var = "sample")

  all_classes <- unique(df_catalogues$class)

  # Loop through each class of signature, build a umap, and add to the list
  ls_umaps <- list()
  for (curr_class in all_classes){

    if(v) cli::cli_progress_step("Bulding UMAP for class {curr_class}")
    class_specific_table <- df_umap_tbl |>
      dplyr::select(starts_with(paste0(curr_class, "|")))

    # Prepare dataset (with all samples which were run for this type of signature analysis)
    n_samples <- nrow(class_specific_table)
    class_specific_table_no_missing <- class_specific_table |>
      dplyr::filter(dplyr::if_all(dplyr::everything(), ~!is.na(.)))
    n_samples_with_tallies_for_currrent_class <- nrow(class_specific_table_no_missing)
    if(n_samples_with_tallies_for_currrent_class != n_samples) cli::cli_alert_info("{n_samples_with_tallies_for_currrent_class}/{n_samples} have had features extracted and counted as part of {curr_class} signature analysis. Umap will be built using only these {n_samples_with_tallies_for_currrent_class}")

    # Check umap_n_neighbours paramater is appropriate for this dataset
    if(n_samples_with_tallies_for_currrent_class < umap_n_neighbours) {
      cli::cli_abort("{.arg umap_n_neighbours} must be smaller than the number of samples in your dataset (<={n_samples_with_tallies_for_currrent_class})")
    }


    # Perform UMAP
    umap <- uwot::umap(class_specific_table_no_missing, n_neighbors = umap_n_neighbours, ret_model = TRUE, seed = seed, batch=TRUE)

    # Add column order to model (so we can check it when projecting new data)
    umap <- uwot_add_column_order(umap, colnames(class_specific_table_no_missing))

    # Write Umap to compressed Rds file (will be used in predictions)
    if(v) cli::cli_progress_step("Adding {curr_class} umap reference to list")
    ls_umaps[[curr_class]] <- umap
  }

  if(v) cli::cli_h2("UMAP from all signature classes")

  # Do the same UMAP using ALL features across all signature types
  if(v) cli::cli_progress_step("Building a single UMAP using ALL features accross all signature classes")

  # Prepare dataset (with all samples which were run for this type of signature analysis)
  n_samples <- nrow(df_umap_tbl)
  class_specific_table_no_missing <- df_umap_tbl |>
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~!is.na(.)))
  n_samples_with_tallies_for_currrent_class <- nrow(class_specific_table_no_missing)
  if(n_samples_with_tallies_for_currrent_class != n_samples) cli::cli_alert_info("{n_samples_with_tallies_for_currrent_class}/{n_samples} have had features extracted and counted for all signature analyses. Umap will be built using only these {n_samples_with_tallies_for_currrent_class}")

  # Check theres enough samples to perform at all
  if(n_samples_with_tallies_for_currrent_class == 0){
    cli::cli_abort("No samples have been analysed against {.strong all} of the following signature classes [{all_classes}]")
  }

  # Check umap_n_neighbours paramater is appropriate for this dataset
  if(n_samples_with_tallies_for_currrent_class < umap_n_neighbours) {
    cli::cli_abort("{.arg umap_n_neighbours} must be smaller than the number of samples in your dataset (<={n_samples_with_tallies_for_currrent_class})")
  }

  # Perform UMAP
  umap <- uwot::umap(class_specific_table_no_missing, n_neighbors = umap_n_neighbours, ret_model = TRUE, seed = seed, batch=TRUE)

  # Add column order (so we can check it when projecting new data)
  umap <- uwot_add_column_order(umap, colnames(class_specific_table_no_missing))

  # Add to
  ls_umaps[["all_classes"]] <- umap

  # Serialise umaps (need to be loaded with uwot::load_uwot() for projecting new data)
  lapply(names(ls_umaps), \(curr_class){
    curr_umap_outfile <- paste0(outfile_umap_prefix, ".", curr_class)
    if(v) cli::cli_progress_step("Writing {curr_class} umap reference to {.file {curr_umap_outfile}}")

    uwot::save_uwot(
      model = ls_umaps[[curr_class]],
      file = curr_umap_outfile,
      unload = TRUE # drop umap from memory since we won't use it again here
    )
    }
  )
}

# Add column order to a umap model object
uwot_add_column_order <- function(umap, column_order){
  umap$column_order <- column_order
  return(umap)
}

build_exposures_reference_set <- function(df_exposures, outfile, format, verbose=TRUE){
  v <- verbose
  if(v) cli::cli_alert_info("Reading Signature Exposure Results")
  if(format == "parquet"){
    df_exposures |>
      dplyr::group_by(class,sample) |>
      arrow::write_dataset(outfile, format = "parquet")
  }
  else{
    df_exposures |>
      readr::write_csv(outfile)
  }

  if(v) cli::cli_alert_success("Exposure refmatrix successfully created")
}

build_bootstrap_reference_set <- function(df_bootstraps, outfile, format, verbose=TRUE){
  v <- verbose
  if(v) cli::cli_alert_info("Reading Signature Exposure Results")
  if(format == "parquet"){
    df_bootstraps |>
      dplyr::group_by(class,sample) |>
      arrow::write_dataset(outfile, format = "parquet")
  }
  else{
    df_bootstraps |>
      readr::write_csv(outfile)
  }

  if(v) cli::cli_alert_success("Exposure refmatrix successfully created")
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
    stats::setNames(keys)

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
