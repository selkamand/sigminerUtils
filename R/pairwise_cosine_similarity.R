
util_compute_cosine_similarity_sql <- function(sample1, sample2, class, con, assume_sensible_input = TRUE){

  s1_query <- paste0("SELECT count FROM decompositions WHERE sampleId = '", sample1,"' AND class = '",class,"'")
  s2_query <- paste0("SELECT count FROM decompositions WHERE sampleId = '", sample2,"' AND class = '",class,"'")

  s1_decomp <- DBI::dbGetQuery(conn = con, statement = s1_query)[[1]]
  s2_decomp <- DBI::dbGetQuery(conn = con, statement = s2_query)[[1]]

  ## Otherwise Compute Similarity (Note if either count vector is all zeros, will return NA)
  lsa::cosine(s1_decomp, s2_decomp)
}

#' Compute Pairwise Similarity Table from Database
#'
#' @param path_db path to sqlite database
#' @param cores number of cores to use
#' @param overwrite should we only calculate similarity for pairs not already in the pairwiseSimilarity table?
#' @param assume_sensible_input should we speed up cosine sim by ~30x by assuming input we supply will be appropriate (channels in the same order)
#' @param low_ram run in a low-ram version where instead of loading all decompositions for all samples into memory, we repeatedly query the sqlite db to fetch each sample
#' @param verbose verbose
#'
#' @return a data.frame with pairwise cosine similarities (sample 1 & sample 2 columns based on R 'sort' order of sampleIds & upper triangle only)
#' @export
#'
sig_pairwise_similarity <- function(path_db, cores = future::availableCores(), overwrite = FALSE, assume_sensible_input = TRUE, verbose = TRUE, low_ram = FALSE){

  # Connect to SQLite DB
  conn <- DBI::dbConnect(RSQLite::SQLite() ,path_db)

  db_decompositions <- if (low_ram) { dplyr::tbl(conn, 'decompositions') } else {db_decompositions <- DBI::dbReadTable(conn, "decompositions")}
  df_sample <- DBI::dbReadTable(conn, "sample")

  df_pairwise_sim <- DBI::dbReadTable(conn, "pairwiseSimilarity")

  # Grab all Samples and Classes
  if (verbose) cli::cli_progress_step('Preparing the list of Sample1-Sample2-Class combinations to compute similarity for')
  samples <- unique(df_sample[['sampleId']])
  classes <- unique(db_decompositions |> dplyr::distinct(class) |> dplyr::collect() |> dplyr::pull(class))

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
  if (verbose) cli::cli_progress_step('Computing similarity using {cores} core{?s}')
  df_rows_to_add[['cosine_similarity']] <- parallel::mcMap(
    f = \(sample1, sample2, class) { util_compute_cosine_similarity_sql(sample1, sample2, class, conn, assume_sensible_input) },
    l[['sample1']],
    l[['sample2']],
    l[['class']],
    mc.cores = cores
  ) |> unlist()

  # Drop Unused Columns
  df_rows_to_add  <- df_rows_to_add |>
    dplyr::select(-pair_id)

  # Disconnect from DB
  DBI::dbDisconnect(conn)

  # Return Similarity Table
  return(df_rows_to_add)
}
