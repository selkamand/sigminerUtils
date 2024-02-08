
# Fetch Decompositions from db and compute sample similarity
util_compute_cosine_similarity <- function(sample1, sample2, class, df_decomps, assume_sensible_input = TRUE){
  s1_decomp = df_decomps |>
    dplyr::filter(sampleId == sample1, class == !!class) |>
    dplyr::rename(type = class)

  s2_decomp = df_decomps |>
    dplyr::filter(sampleId == sample2, class == !!class) |>
    dplyr::rename(type = class)

  # If any NAs in fraction (indicate zero mutations in the class)
  # Then Cosine Similarity Should be NA (sinc its non-computable
  if (anyNA(c(s1_decomp[['fraction']], s2_decomp[['fraction']])))
    return(NA_real_)

  ## Otherwise Compute Similarity
  sigstats::sig_cosine_similarity(s1_decomp, s2_decomp, assume_sensible_input = assume_sensible_input)
}

extract_sigverse_decomps <- function(sample, class, df_decomps, overwrite){
  df_decomps |>
    dplyr::filter(sampleId == sample, class == !!class) |>
    dplyr::rename(type = class)
}

#' Compute Pairwise Similarity Table from Database
#'
#' @param path_db path to sqlite database
#' @param cores number of cores to use
#' @param overwrite should we only calculate similarity for pairs not already in the pairwiseSimilarity table?
#' @param assume_sensible_input should we speed up cosine sim by ~30x by assuming input we supply will be appropriate (channels in the same order)
#' @param verbose verbose
#'
#' @return a data.frame with pairwise cosine similarities (sample 1 & sample 2 columns based on R 'sort' order of sampleIds & upper triangle only)
#' @export
#'
sig_pairwise_similarity <- function(path_db, cores = future::availableCores(), overwrite = FALSE, assume_sensible_input = TRUE, verbose = TRUE){

  # Connect to SQLite DB
  conn <- DBI::dbConnect(RSQLite::SQLite() ,path_db)

  df_decompositions <- DBI::dbReadTable(conn, "decompositions")
  df_pairwise_sim <- DBI::dbReadTable(conn, "pairwiseSimilarity")


  # Grab all Samples and Classes
  samples <- unique(df_decompositions[["sampleId"]])
  classes <- unique(df_decompositions[["class"]])

  # Create a data.frame representing the unique sample combinations to compute similarity between
  df_rows_to_add <- expand.grid(sample1 = samples, sample2 = samples, class = classes) |>
    dplyr::tibble() |>
    dplyr::rowwise() |>
    dplyr::filter(
      # Choose order of samples based on the order they appear in R base 'sort'
      identical(c(sample1, sample2), sort(c(sample1, sample2))),
      # Exclude self-similarity
      sample1 != sample2
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(pair_id = paste(sample1, sample2, class))

  # Check Which Samples combinations are missing from 'pairwiseSimilarity' table
  if (!overwrite) {
    pairs_already_in_db <- unique(paste(df_pairwise_sim[["sample1"]], df_pairwise_sim[["sample2"]], df_pairwise_sim[["class"]]))

    df_rows_to_add <- df_rows_to_add |>
      dplyr::filter(!pair_id %in% pairs_already_in_db)

    if (verbose) cli::cli_alert('Found {nrow(df_rows_to_add)} unique combinations of sampleId{?s} & classes to add to {?s} pairwiseSimilarity table')
  }
  else{
    if (verbose) cli::cli_alert('Found {nrow(df_rows_to_add)} unique combinations of sampleId{?s} & classes to add to {?s} pairwiseSimilarity table')
  }


  # Create a list of inputs that we'll iterate over using Map/purrr::pmap
  l <- df_rows_to_add |>
    dplyr::select(sample1, sample2, class) |>
    as.list() |>
    lapply(as.character)

  # Compute Cosine Similarity
  df_rows_to_add[['cosine_similarity']] <- parallel::mcMap(
    f = \(sample1, sample2, class) { util_compute_cosine_similarity(sample1, sample2, class, df_decompositions, assume_sensible_input) },
    l[['sample1']],
    l[['sample2']],
    l[['class']],
    mc.cores = cores
  ) |> unlist()

  # Drop Unused Columns
  df_rows_to_add  <- df_rows_to_add |>
    dplyr::select(-pair_id)


  return(df_rows_to_add)
}
