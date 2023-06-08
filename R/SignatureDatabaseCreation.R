
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
#'
sig_create_database <- function(sqlite_db, overwrite = TRUE){
  path_sql <- here::here('SQL/CreateSignatureDB.sql')

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
#'
#' @return invisible(NULL). This function is run for its side effects
#' @export
#'
sig_add_to_database <- function(signature_directory, sqlite_db){
  assertions::assert_directory_exists(signature_directory)
  assertions::assert_file_exists(sqlite_db)

  df_files <- signature_directory_to_table(signature_directory)

  # Exposure tables
  cli::cli_progress_step("Reading data from {.path {signature_directory}}")

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

  cli::cli_progress_step("Connecting to the signature database {.path {sqlite_db}}")
  con <- RSQLite::dbConnect(drv = RSQLite::SQLite(), sqlite_db)

  cli::cli_progress_step("Appending exposure data to database {.path {sqlite_db}}")
  DBI::dbWriteTable(conn = con, name = "cosmicExposures", df_exposures, append = TRUE, row.names = FALSE)

  return(invisible(NULL))
}
