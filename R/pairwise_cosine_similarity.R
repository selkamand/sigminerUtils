sig_update_database_sqlite <- function(path_db){

  # Connect to an sqlite signature database and for any samples lacking pairwise cosine similarity - update it.

}

# Fetch Decompositions from db and compute sample similarity
util_compute_cosine_similarity <- function(sample1, sample2, class, df_decomps){
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
  sigstats::sig_cosine_similarity(s1_decomp, s2_decomp)
}

extract_sigverse_decomps <- function(sample, class, df_decomps, overwrite){
  df_decomps |>
    dplyr::filter(sampleId == sample, class == !!class) |>
    dplyr::rename(type = class)
}

sig_pairwise_similarity <- function(path_db, overwrite = FALSE, verbose = TRUE){

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

    if (verbose) cli::cli_alert('Found {length(df_rows_to_add)} unique combinations of sampleId{?s} & classes to add to {?s} pairwiseSimilarity table')
  }
  else{
    if (verbose) cli::cli_alert('Found {length(df_rows_to_add)} unique combinations of sampleId{?s} & classes to add to {?s} pairwiseSimilarity table')
  }


  # Compute Cosine Similarity
  df_rows_to_add[['cosine_similarity']] <- purrr::pmap_dbl(
    .l = df_rows_to_add |> dplyr::select(sample1, sample2, class),
    .f = \(sample1, sample2, class) { util_compute_cosine_similarity(sample1, sample2, class, df_decompositions) }
  )

  # Drop Unused Columns
  df_rows_to_add  <- df_rows_to_add |>
    dplyr::select(-pair_id)


  return(df_rows_to_add)
}
