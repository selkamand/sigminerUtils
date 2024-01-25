
# Utilities ---------------------------------------------------------------
split_sql_string_by_command <- function(command){
  unlist(strsplit(x=command, split = ';[[:space:]]*'))
}


signature_directory_to_table <- function(signature_directory){
  assertions::assert_directory_exists(x = signature_directory)
  paths = dir(signature_directory, full.names = TRUE)

  df_inputs <- dplyr::tibble(
    filepath = paths,
    filename = basename(paths),
    signature_set = sub(x=filename, pattern = "_.*", ""),
    reference_genome = sub(x=filename, pattern = "^.*\\.(.*?)\\..*$", "\\1"),
    output_type = sub(x=filename, pattern = "^.*\\..*.*\\.(.*).csv$", "\\1"),
  )

  df_inputs <- df_inputs |> dplyr::mutate(
    output_type = ifelse(output_type == "expo","exposures", output_type),
    tablename = paste0(signature_set,"_", output_type,"_",reference_genome)
  )

  return(df_inputs)
}


# Key Commands ------------------------------------------------------------

#' Create a Signature Database
#'
#' @param sqlite_db  path to sql db
#' @param overwrite  overwrite old sqlite db
#'
#' @return invisible(NULL) run for its side effects
#' @export
sig_create_database <- function(sqlite_db, overwrite = TRUE){
  path_sql <- system.file('sql/CreateSignatureDB.sql', package = "sigminerUtils")

  cli::cli_h2('Create Signature Database')

  cli::cli_progress_step("Checking if db exists")
  if(file.exists(sqlite_db)){
    if (!overwrite){
      cli::cli_abort('Signature database already exists. Delete {.path {sqlite_db}} or set {.code overwrite = TRUE} to replace existing database')
    }
    else{
      cli::cli_alert_warning("Signature database already exists ({.path {sqlite_db}})")
      res <- utils::askYesNo(msg = "Are you sure you want to overwrite?")
      if (res) {
        cli::cli_progress_step('Deleting old signature database')
        unlink(sqlite_db)
      }
      else
        cli::cli_abort('Signature database already exists. Delete {.path {sqlite_db}} or set {.code overwrite = TRUE} to replace existing database')
    }
  }

  cli::cli_progress_step("Connecting to the signature database {.path {sqlite_db}}")
  con = RSQLite::dbConnect(drv = RSQLite::SQLite(), sqlite_db)

  cli::cli_progress_step("Reading Database Creation Script {.path {path_sql}}")
  database_creation_script <- readr::read_file(path_sql)
  database_creation_script <- split_sql_string_by_command(database_creation_script)

  cli::cli_progress_step("Running Database Creation Script")
  for (command in database_creation_script){
    DBI::dbExecute(conn = con, statement = command)
  }
}


#' Add data to sqlite DB
#'
#' Load signature data into sqlite DB
#'
#' @param signature_directory path to the directory produced by [sig_analyse_mutations()]
#' @param sqlite_db path to the sqlite database produced by [(sig_create_database)]
#' @param ref reference genome: one of hg19 or hg38 (string)
#' @return invisible(NULL). This function is run for its side effects
#' @export
sig_add_to_database <- function(signature_directory, sqlite_db, ref = c("hg19", "hg38")){
  assertions::assert_directory_exists(signature_directory)
  assertions::assert_file_exists(sqlite_db)
  ref <- rlang::arg_match(ref)

  df_files <- signature_directory_to_table(signature_directory)

  #refs_found <- unique(df_files[['reference_genome']])

  df_files <- df_files |> dplyr::filter(reference_genome == ref)

  assertions::assert_greater_than(nrow(df_files), minimum = 0, msg = "Failed to find any signature files for ref genome {ref} in {.path {signature_directory}}")

  cli::cli_progress_step("Reading data from {.path {signature_directory}}")

  # Exposure
  df_exposures <-  df_files |>
    dplyr::filter(output_type == "exposures") |>
    dplyr::pull(filepath) |>
    readr::read_csv()

  df_exposures <- df_exposures |>
    dplyr::rename(method=Method,
                  sampleId = SampleID,
                  signature=Sig,
                  contribution=Contribution,
                  contributionRelative = ContributionRelative,
                  optimal = IsOptimal,
                  type = Type)


  # Decompositions
  df_decomposition <-  df_files |>
    dplyr::filter(output_type == "decomposition") |>
    dplyr::pull(filepath) |>
    readr::read_csv(id = 'class')

  df_decomposition <- df_decomposition |>
    dplyr::rename(
      sampleId = "SampleID",
      channel = "Context",
      count = "Count",
      fraction = "CountRelative"
      ) |>
    dplyr::mutate(class = sub(x=basename(class), pattern = "_.*",replacement = ""))

  # cosmicErrorAndCosine
  df_error <- df_files |>
    dplyr::filter(output_type == "error") |>
    dplyr::pull(filepath) |>
    readr::read_csv(id = 'class')

  df_cosine <- df_files |>
    dplyr::filter(output_type == "cosine") |>
    dplyr::pull(filepath) |>
    readr::read_csv(id = 'class')

  df_error <- df_error |>
    dplyr::rename(
      sampleId = "SampleID",
      method = "Method",
      type = "Type",
      optimal = "IsOptimal",
      error = "Errors"
      ) |>
    dplyr::mutate(class = sub(x=basename(class), pattern = "_.*",replacement = ""))

  df_cosine <- df_cosine |>
    dplyr::rename(
      sampleId = "SampleID",
      method = "Method",
      type = "Type",
      optimal = "IsOptimal",
      cosine = "Cosine"
    ) |>
    dplyr::mutate(
      class = sub(x=basename(class), pattern = "_.*",replacement = "")
    )

  df_error_and_cosine <- dplyr::left_join(df_error, df_cosine, by = dplyr::join_by(class, sampleId, method, type, optimal), keep = FALSE)

  # pVal
  df_pval <- df_files |>
    dplyr::filter(output_type == "p_val") |>
    dplyr::pull(filepath) |>
    readr::read_csv(id = 'class')

  df_pval <- df_pval |>
    dplyr::rename(
      sampleId = "SampleID",
      method = "Method",
      signature = "Sig",
      threshold = "Threshold",
      p = "Pvalue"
    ) |>
    dplyr::mutate(
      class = sub(x=basename(class), pattern = "_.*",replacement = "")
    )

  cli::cli_progress_step("Connecting to the signature database {.path {sqlite_db}}")
  con <- RSQLite::dbConnect(drv = RSQLite::SQLite(), sqlite_db)

  cli::cli_progress_step("Appending exposure data to database {.path {sqlite_db}}")
  DBI::dbWriteTable(conn = con, name = "cosmicExposures", df_exposures, append = TRUE, row.names = FALSE)

  cli::cli_progress_step("Appending decomposition data to database {.path {sqlite_db}}")
  DBI::dbWriteTable(conn = con, name = "decompositions", df_decomposition, append = TRUE, row.names = FALSE)

  cli::cli_progress_step("Appending error and cosine data to database {.path {sqlite_db}}")
  DBI::dbWriteTable(conn = con, name = "cosmicErrorAndCosine", df_error_and_cosine, append = TRUE, row.names = FALSE)

  cli::cli_progress_step("Appending P-value data to database {.path {sqlite_db}}")
  DBI::dbWriteTable(conn = con, name = "cosmicPvalues", df_pval, append = TRUE, row.names = FALSE)



  return(invisible(NULL))
}
