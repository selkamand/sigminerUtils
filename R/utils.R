#' Parse Signature Database.
#'
#' If db is a string, assumes it is a path to a csv_tidy formatted file and parses it into a sigminer (wide) format.
#' If db is not a string, returns db.
#'
#' @param db a signature collection
#' @param dbtype a string describing the db type (e.g. 'db_sbs', 'db_sv', etc)
#'
#' @return sigminer formatted signature collection if db is a string, else returns db
#'
db_read_if_filepath <- function(db, dbtype = "A"){
  if(is.character(db) & length(db) == 1) {
    cli::cli_alert_info("{dbtype} signature collection was supplied as a string. Attempting to read as a file (assuming csv_tidy format)")
    signature_collection <- sigstash::sig_read_signatures(filepath = db, format = "csv_tidy")
    sigstash::sig_collection_to_sigminer(signature_collection)
  }
  else
    return(db)
}

vcf2versions <- function(path_to_vcf){
  versions <- readLines(path_to_vcf) |> grep(x = _, "^##.*[V]ersion", value = T)
  versions <- gsub(x = versions, "#", "")
  versions_collapsed <- paste0(versions, collapse = "\n")
  return(versions_collapsed)
}


capture_messages <- function(logfile, tee=TRUE, expr) {
  messages <- utils::capture.output(expr, type = "message")
  cat(messages, file = logfile, append = TRUE, sep = "\n")
  if(tee) cat(messages, sep = "\n")
}

silence_messages <- function(verbose, expr){
  if(!verbose)
    suppressMessages(expr)
  else
    expr
}

#' Write compressed CSV
#'
#' @inheritParams utils::write.csv
#'
write_compressed_csv <- function(x, file, row.names = FALSE){

  # Remove .gz extensions
  file <- sub(x=file, pattern = "\\.gz", replacement = "")

  utils::write.csv(x, file, row.names = row.names)
  R.utils::bzip2(file, ext="gz", overwrite = TRUE)
}

# Convert trycatch errors to NULL
try_error_to_null <- function(obj){
  if("try-error" %in% class(obj)) return(NULL)
  else return(obj)
}

#' Tally
#'
#' A wrapper for [sigminer::sig_tally()] that adds 2 features.
#' 1. If supplied with a MAF with none of the type of mutation for the given 'MODE' instead of erroring returns an empty matrix.
#'
#' @inheritParams sigminer::sig_tally
#'
#' @inherit sigminer::sig_tally return
#'
tally <- function(object,
                  mode,
                  ref_genome,
                  cores){

  tryCatch(
    expr = {
      sigminer::sig_tally(
        object = object,
        mode = mode,
        ref_genome = ref_genome,
        keep_only_matrix = FALSE,
        cores = cores
      )
    },
    error = function(err){
      err_string = as.character(err)
      no_mutations <- any(grepl(x = err_string, pattern = "Zero [[:alpha:]]+ to analyze"))
      if (no_mutations) {
        cli::cli_warn("Zero Indels to analyze")
        return(NULL)
      }
      else
        stop(err)
    }
  )
}


#' Post-process a tally matrix.
#'
#' This function performs two operations on the input matrix `mx`:
#'
#' \enumerate{
#'   \item If `mx` is `NULL` (which can occur when the supplied MAF file doesn't contain any mutations
#'   of the relevant class), the function produces a default numeric matrix of size `samples × channels`
#'   where all values are 0.
#'
#'   \item The matrix is transposed to become a `channels × samples` matrix.
#' }
#'
#' @param mx A numeric matrix, or `NULL`. The input tally matrix to be processed.
#' @param class A character vector specifying the class of mutations.
#' @param samples An integer specifying the number of samples.
#'
#' @return A numeric matrix with dimensions `channels × samples`.
#'
postprocess_tally_matrix <- function(mx, class, samples){
  mx <-  replace_null_with_empty_matrix(mx, class, samples)
  mx <- t(mx)
}

replace_null_with_empty_matrix <- function(x, class, samples){
  if(!is.null(x)) return(x)

  sigminer_initiate_empty(class, samples)
}

#' Return an 0-initialised matrix
#'
#' Return a 0-initialised matrix of Sample IDs X Signature Channels
#'
#' @param class class of signatures to return (e.g. one of SBS96, SBS1536, ID83, ID28, DBS78, and DBS1248). Matrix returned will have 1 column per channel of the selected signature class.
#' @param samples a vector of sample IDs will become rows of the returned matrix
sigminer_initiate_empty <- function(class, samples){

  ls_channels <- sigminer_channels()
  classes <- names(ls_channels)
  assertions::assert_subset(class, classes)

  channels <- ls_channels[[class]]

  nsamples <- length(samples)
  nchannels <- length(channels)

  mx <- matrix(
    data = rep(0, times = nsamples * nchannels), nrow = nsamples, dimnames = list(
      samples, channels
    )
  )

  return(mx)
}

sigminer_channels <- function() {
  list(
    "SBS96" = c(
      "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
      "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
      "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
      "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
      "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
      "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
      "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
      "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
      "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
      "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
      "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
      "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
      "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
      "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
      "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
      "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
    ),
    "SBS1536" = c(
      "AA[C>A]AA", "AA[C>A]AC", "AA[C>A]AG", "AA[C>A]AT", "AA[C>A]CA",
      "AA[C>A]CC", "AA[C>A]CG", "AA[C>A]CT", "AA[C>A]GA", "AA[C>A]GC",
      "AA[C>A]GG", "AA[C>A]GT", "AA[C>A]TA", "AA[C>A]TC", "AA[C>A]TG",
      "AA[C>A]TT", "AA[C>G]AA", "AA[C>G]AC", "AA[C>G]AG", "AA[C>G]AT",
      "AA[C>G]CA", "AA[C>G]CC", "AA[C>G]CG", "AA[C>G]CT", "AA[C>G]GA",
      "AA[C>G]GC", "AA[C>G]GG", "AA[C>G]GT", "AA[C>G]TA", "AA[C>G]TC",
      "AA[C>G]TG", "AA[C>G]TT", "AA[C>T]AA", "AA[C>T]AC", "AA[C>T]AG",
      "AA[C>T]AT", "AA[C>T]CA", "AA[C>T]CC", "AA[C>T]CG", "AA[C>T]CT",
      "AA[C>T]GA", "AA[C>T]GC", "AA[C>T]GG", "AA[C>T]GT", "AA[C>T]TA",
      "AA[C>T]TC", "AA[C>T]TG", "AA[C>T]TT", "AA[T>A]AA", "AA[T>A]AC",
      "AA[T>A]AG", "AA[T>A]AT", "AA[T>A]CA", "AA[T>A]CC", "AA[T>A]CG",
      "AA[T>A]CT", "AA[T>A]GA", "AA[T>A]GC", "AA[T>A]GG", "AA[T>A]GT",
      "AA[T>A]TA", "AA[T>A]TC", "AA[T>A]TG", "AA[T>A]TT", "AA[T>C]AA",
      "AA[T>C]AC", "AA[T>C]AG", "AA[T>C]AT", "AA[T>C]CA", "AA[T>C]CC",
      "AA[T>C]CG", "AA[T>C]CT", "AA[T>C]GA", "AA[T>C]GC", "AA[T>C]GG",
      "AA[T>C]GT", "AA[T>C]TA", "AA[T>C]TC", "AA[T>C]TG", "AA[T>C]TT",
      "AA[T>G]AA", "AA[T>G]AC", "AA[T>G]AG", "AA[T>G]AT", "AA[T>G]CA",
      "AA[T>G]CC", "AA[T>G]CG", "AA[T>G]CT", "AA[T>G]GA", "AA[T>G]GC",
      "AA[T>G]GG", "AA[T>G]GT", "AA[T>G]TA", "AA[T>G]TC", "AA[T>G]TG",
      "AA[T>G]TT", "AC[C>A]AA", "AC[C>A]AC", "AC[C>A]AG", "AC[C>A]AT",
      "AC[C>A]CA", "AC[C>A]CC", "AC[C>A]CG", "AC[C>A]CT", "AC[C>A]GA",
      "AC[C>A]GC", "AC[C>A]GG", "AC[C>A]GT", "AC[C>A]TA", "AC[C>A]TC",
      "AC[C>A]TG", "AC[C>A]TT", "AC[C>G]AA", "AC[C>G]AC", "AC[C>G]AG",
      "AC[C>G]AT", "AC[C>G]CA", "AC[C>G]CC", "AC[C>G]CG", "AC[C>G]CT",
      "AC[C>G]GA", "AC[C>G]GC", "AC[C>G]GG", "AC[C>G]GT", "AC[C>G]TA",
      "AC[C>G]TC", "AC[C>G]TG", "AC[C>G]TT", "AC[C>T]AA", "AC[C>T]AC",
      "AC[C>T]AG", "AC[C>T]AT", "AC[C>T]CA", "AC[C>T]CC", "AC[C>T]CG",
      "AC[C>T]CT", "AC[C>T]GA", "AC[C>T]GC", "AC[C>T]GG", "AC[C>T]GT",
      "AC[C>T]TA", "AC[C>T]TC", "AC[C>T]TG", "AC[C>T]TT", "AC[T>A]AA",
      "AC[T>A]AC", "AC[T>A]AG", "AC[T>A]AT", "AC[T>A]CA", "AC[T>A]CC",
      "AC[T>A]CG", "AC[T>A]CT", "AC[T>A]GA", "AC[T>A]GC", "AC[T>A]GG",
      "AC[T>A]GT", "AC[T>A]TA", "AC[T>A]TC", "AC[T>A]TG", "AC[T>A]TT",
      "AC[T>C]AA", "AC[T>C]AC", "AC[T>C]AG", "AC[T>C]AT", "AC[T>C]CA",
      "AC[T>C]CC", "AC[T>C]CG", "AC[T>C]CT", "AC[T>C]GA", "AC[T>C]GC",
      "AC[T>C]GG", "AC[T>C]GT", "AC[T>C]TA", "AC[T>C]TC", "AC[T>C]TG",
      "AC[T>C]TT", "AC[T>G]AA", "AC[T>G]AC", "AC[T>G]AG", "AC[T>G]AT",
      "AC[T>G]CA", "AC[T>G]CC", "AC[T>G]CG", "AC[T>G]CT", "AC[T>G]GA",
      "AC[T>G]GC", "AC[T>G]GG", "AC[T>G]GT", "AC[T>G]TA", "AC[T>G]TC",
      "AC[T>G]TG", "AC[T>G]TT", "AG[C>A]AA", "AG[C>A]AC", "AG[C>A]AG",
      "AG[C>A]AT", "AG[C>A]CA", "AG[C>A]CC", "AG[C>A]CG", "AG[C>A]CT",
      "AG[C>A]GA", "AG[C>A]GC", "AG[C>A]GG", "AG[C>A]GT", "AG[C>A]TA",
      "AG[C>A]TC", "AG[C>A]TG", "AG[C>A]TT", "AG[C>G]AA", "AG[C>G]AC",
      "AG[C>G]AG", "AG[C>G]AT", "AG[C>G]CA", "AG[C>G]CC", "AG[C>G]CG",
      "AG[C>G]CT", "AG[C>G]GA", "AG[C>G]GC", "AG[C>G]GG", "AG[C>G]GT",
      "AG[C>G]TA", "AG[C>G]TC", "AG[C>G]TG", "AG[C>G]TT", "AG[C>T]AA",
      "AG[C>T]AC", "AG[C>T]AG", "AG[C>T]AT", "AG[C>T]CA", "AG[C>T]CC",
      "AG[C>T]CG", "AG[C>T]CT", "AG[C>T]GA", "AG[C>T]GC", "AG[C>T]GG",
      "AG[C>T]GT", "AG[C>T]TA", "AG[C>T]TC", "AG[C>T]TG", "AG[C>T]TT",
      "AG[T>A]AA", "AG[T>A]AC", "AG[T>A]AG", "AG[T>A]AT", "AG[T>A]CA",
      "AG[T>A]CC", "AG[T>A]CG", "AG[T>A]CT", "AG[T>A]GA", "AG[T>A]GC",
      "AG[T>A]GG", "AG[T>A]GT", "AG[T>A]TA", "AG[T>A]TC", "AG[T>A]TG",
      "AG[T>A]TT", "AG[T>C]AA", "AG[T>C]AC", "AG[T>C]AG", "AG[T>C]AT",
      "AG[T>C]CA", "AG[T>C]CC", "AG[T>C]CG", "AG[T>C]CT", "AG[T>C]GA",
      "AG[T>C]GC", "AG[T>C]GG", "AG[T>C]GT", "AG[T>C]TA", "AG[T>C]TC",
      "AG[T>C]TG", "AG[T>C]TT", "AG[T>G]AA", "AG[T>G]AC", "AG[T>G]AG",
      "AG[T>G]AT", "AG[T>G]CA", "AG[T>G]CC", "AG[T>G]CG", "AG[T>G]CT",
      "AG[T>G]GA", "AG[T>G]GC", "AG[T>G]GG", "AG[T>G]GT", "AG[T>G]TA",
      "AG[T>G]TC", "AG[T>G]TG", "AG[T>G]TT", "AT[C>A]AA", "AT[C>A]AC",
      "AT[C>A]AG", "AT[C>A]AT", "AT[C>A]CA", "AT[C>A]CC", "AT[C>A]CG",
      "AT[C>A]CT", "AT[C>A]GA", "AT[C>A]GC", "AT[C>A]GG", "AT[C>A]GT",
      "AT[C>A]TA", "AT[C>A]TC", "AT[C>A]TG", "AT[C>A]TT", "AT[C>G]AA",
      "AT[C>G]AC", "AT[C>G]AG", "AT[C>G]AT", "AT[C>G]CA", "AT[C>G]CC",
      "AT[C>G]CG", "AT[C>G]CT", "AT[C>G]GA", "AT[C>G]GC", "AT[C>G]GG",
      "AT[C>G]GT", "AT[C>G]TA", "AT[C>G]TC", "AT[C>G]TG", "AT[C>G]TT",
      "AT[C>T]AA", "AT[C>T]AC", "AT[C>T]AG", "AT[C>T]AT", "AT[C>T]CA",
      "AT[C>T]CC", "AT[C>T]CG", "AT[C>T]CT", "AT[C>T]GA", "AT[C>T]GC",
      "AT[C>T]GG", "AT[C>T]GT", "AT[C>T]TA", "AT[C>T]TC", "AT[C>T]TG",
      "AT[C>T]TT", "AT[T>A]AA", "AT[T>A]AC", "AT[T>A]AG", "AT[T>A]AT",
      "AT[T>A]CA", "AT[T>A]CC", "AT[T>A]CG", "AT[T>A]CT", "AT[T>A]GA",
      "AT[T>A]GC", "AT[T>A]GG", "AT[T>A]GT", "AT[T>A]TA", "AT[T>A]TC",
      "AT[T>A]TG", "AT[T>A]TT", "AT[T>C]AA", "AT[T>C]AC", "AT[T>C]AG",
      "AT[T>C]AT", "AT[T>C]CA", "AT[T>C]CC", "AT[T>C]CG", "AT[T>C]CT",
      "AT[T>C]GA", "AT[T>C]GC", "AT[T>C]GG", "AT[T>C]GT", "AT[T>C]TA",
      "AT[T>C]TC", "AT[T>C]TG", "AT[T>C]TT", "AT[T>G]AA", "AT[T>G]AC",
      "AT[T>G]AG", "AT[T>G]AT", "AT[T>G]CA", "AT[T>G]CC", "AT[T>G]CG",
      "AT[T>G]CT", "AT[T>G]GA", "AT[T>G]GC", "AT[T>G]GG", "AT[T>G]GT",
      "AT[T>G]TA", "AT[T>G]TC", "AT[T>G]TG", "AT[T>G]TT", "CA[C>A]AA",
      "CA[C>A]AC", "CA[C>A]AG", "CA[C>A]AT", "CA[C>A]CA", "CA[C>A]CC",
      "CA[C>A]CG", "CA[C>A]CT", "CA[C>A]GA", "CA[C>A]GC", "CA[C>A]GG",
      "CA[C>A]GT", "CA[C>A]TA", "CA[C>A]TC", "CA[C>A]TG", "CA[C>A]TT",
      "CA[C>G]AA", "CA[C>G]AC", "CA[C>G]AG", "CA[C>G]AT", "CA[C>G]CA",
      "CA[C>G]CC", "CA[C>G]CG", "CA[C>G]CT", "CA[C>G]GA", "CA[C>G]GC",
      "CA[C>G]GG", "CA[C>G]GT", "CA[C>G]TA", "CA[C>G]TC", "CA[C>G]TG",
      "CA[C>G]TT", "CA[C>T]AA", "CA[C>T]AC", "CA[C>T]AG", "CA[C>T]AT",
      "CA[C>T]CA", "CA[C>T]CC", "CA[C>T]CG", "CA[C>T]CT", "CA[C>T]GA",
      "CA[C>T]GC", "CA[C>T]GG", "CA[C>T]GT", "CA[C>T]TA", "CA[C>T]TC",
      "CA[C>T]TG", "CA[C>T]TT", "CA[T>A]AA", "CA[T>A]AC", "CA[T>A]AG",
      "CA[T>A]AT", "CA[T>A]CA", "CA[T>A]CC", "CA[T>A]CG", "CA[T>A]CT",
      "CA[T>A]GA", "CA[T>A]GC", "CA[T>A]GG", "CA[T>A]GT", "CA[T>A]TA",
      "CA[T>A]TC", "CA[T>A]TG", "CA[T>A]TT", "CA[T>C]AA", "CA[T>C]AC",
      "CA[T>C]AG", "CA[T>C]AT", "CA[T>C]CA", "CA[T>C]CC", "CA[T>C]CG",
      "CA[T>C]CT", "CA[T>C]GA", "CA[T>C]GC", "CA[T>C]GG", "CA[T>C]GT",
      "CA[T>C]TA", "CA[T>C]TC", "CA[T>C]TG", "CA[T>C]TT", "CA[T>G]AA",
      "CA[T>G]AC", "CA[T>G]AG", "CA[T>G]AT", "CA[T>G]CA", "CA[T>G]CC",
      "CA[T>G]CG", "CA[T>G]CT", "CA[T>G]GA", "CA[T>G]GC", "CA[T>G]GG",
      "CA[T>G]GT", "CA[T>G]TA", "CA[T>G]TC", "CA[T>G]TG", "CA[T>G]TT",
      "CC[C>A]AA", "CC[C>A]AC", "CC[C>A]AG", "CC[C>A]AT", "CC[C>A]CA",
      "CC[C>A]CC", "CC[C>A]CG", "CC[C>A]CT", "CC[C>A]GA", "CC[C>A]GC",
      "CC[C>A]GG", "CC[C>A]GT", "CC[C>A]TA", "CC[C>A]TC", "CC[C>A]TG",
      "CC[C>A]TT", "CC[C>G]AA", "CC[C>G]AC", "CC[C>G]AG", "CC[C>G]AT",
      "CC[C>G]CA", "CC[C>G]CC", "CC[C>G]CG", "CC[C>G]CT", "CC[C>G]GA",
      "CC[C>G]GC", "CC[C>G]GG", "CC[C>G]GT", "CC[C>G]TA", "CC[C>G]TC",
      "CC[C>G]TG", "CC[C>G]TT", "CC[C>T]AA", "CC[C>T]AC", "CC[C>T]AG",
      "CC[C>T]AT", "CC[C>T]CA", "CC[C>T]CC", "CC[C>T]CG", "CC[C>T]CT",
      "CC[C>T]GA", "CC[C>T]GC", "CC[C>T]GG", "CC[C>T]GT", "CC[C>T]TA",
      "CC[C>T]TC", "CC[C>T]TG", "CC[C>T]TT", "CC[T>A]AA", "CC[T>A]AC",
      "CC[T>A]AG", "CC[T>A]AT", "CC[T>A]CA", "CC[T>A]CC", "CC[T>A]CG",
      "CC[T>A]CT", "CC[T>A]GA", "CC[T>A]GC", "CC[T>A]GG", "CC[T>A]GT",
      "CC[T>A]TA", "CC[T>A]TC", "CC[T>A]TG", "CC[T>A]TT", "CC[T>C]AA",
      "CC[T>C]AC", "CC[T>C]AG", "CC[T>C]AT", "CC[T>C]CA", "CC[T>C]CC",
      "CC[T>C]CG", "CC[T>C]CT", "CC[T>C]GA", "CC[T>C]GC", "CC[T>C]GG",
      "CC[T>C]GT", "CC[T>C]TA", "CC[T>C]TC", "CC[T>C]TG", "CC[T>C]TT",
      "CC[T>G]AA", "CC[T>G]AC", "CC[T>G]AG", "CC[T>G]AT", "CC[T>G]CA",
      "CC[T>G]CC", "CC[T>G]CG", "CC[T>G]CT", "CC[T>G]GA", "CC[T>G]GC",
      "CC[T>G]GG", "CC[T>G]GT", "CC[T>G]TA", "CC[T>G]TC", "CC[T>G]TG",
      "CC[T>G]TT", "CG[C>A]AA", "CG[C>A]AC", "CG[C>A]AG", "CG[C>A]AT",
      "CG[C>A]CA", "CG[C>A]CC", "CG[C>A]CG", "CG[C>A]CT", "CG[C>A]GA",
      "CG[C>A]GC", "CG[C>A]GG", "CG[C>A]GT", "CG[C>A]TA", "CG[C>A]TC",
      "CG[C>A]TG", "CG[C>A]TT", "CG[C>G]AA", "CG[C>G]AC", "CG[C>G]AG",
      "CG[C>G]AT", "CG[C>G]CA", "CG[C>G]CC", "CG[C>G]CG", "CG[C>G]CT",
      "CG[C>G]GA", "CG[C>G]GC", "CG[C>G]GG", "CG[C>G]GT", "CG[C>G]TA",
      "CG[C>G]TC", "CG[C>G]TG", "CG[C>G]TT", "CG[C>T]AA", "CG[C>T]AC",
      "CG[C>T]AG", "CG[C>T]AT", "CG[C>T]CA", "CG[C>T]CC", "CG[C>T]CG",
      "CG[C>T]CT", "CG[C>T]GA", "CG[C>T]GC", "CG[C>T]GG", "CG[C>T]GT",
      "CG[C>T]TA", "CG[C>T]TC", "CG[C>T]TG", "CG[C>T]TT", "CG[T>A]AA",
      "CG[T>A]AC", "CG[T>A]AG", "CG[T>A]AT", "CG[T>A]CA", "CG[T>A]CC",
      "CG[T>A]CG", "CG[T>A]CT", "CG[T>A]GA", "CG[T>A]GC", "CG[T>A]GG",
      "CG[T>A]GT", "CG[T>A]TA", "CG[T>A]TC", "CG[T>A]TG", "CG[T>A]TT",
      "CG[T>C]AA", "CG[T>C]AC", "CG[T>C]AG", "CG[T>C]AT", "CG[T>C]CA",
      "CG[T>C]CC", "CG[T>C]CG", "CG[T>C]CT", "CG[T>C]GA", "CG[T>C]GC",
      "CG[T>C]GG", "CG[T>C]GT", "CG[T>C]TA", "CG[T>C]TC", "CG[T>C]TG",
      "CG[T>C]TT", "CG[T>G]AA", "CG[T>G]AC", "CG[T>G]AG", "CG[T>G]AT",
      "CG[T>G]CA", "CG[T>G]CC", "CG[T>G]CG", "CG[T>G]CT", "CG[T>G]GA",
      "CG[T>G]GC", "CG[T>G]GG", "CG[T>G]GT", "CG[T>G]TA", "CG[T>G]TC",
      "CG[T>G]TG", "CG[T>G]TT", "CT[C>A]AA", "CT[C>A]AC", "CT[C>A]AG",
      "CT[C>A]AT", "CT[C>A]CA", "CT[C>A]CC", "CT[C>A]CG", "CT[C>A]CT",
      "CT[C>A]GA", "CT[C>A]GC", "CT[C>A]GG", "CT[C>A]GT", "CT[C>A]TA",
      "CT[C>A]TC", "CT[C>A]TG", "CT[C>A]TT", "CT[C>G]AA", "CT[C>G]AC",
      "CT[C>G]AG", "CT[C>G]AT", "CT[C>G]CA", "CT[C>G]CC", "CT[C>G]CG",
      "CT[C>G]CT", "CT[C>G]GA", "CT[C>G]GC", "CT[C>G]GG", "CT[C>G]GT",
      "CT[C>G]TA", "CT[C>G]TC", "CT[C>G]TG", "CT[C>G]TT", "CT[C>T]AA",
      "CT[C>T]AC", "CT[C>T]AG", "CT[C>T]AT", "CT[C>T]CA", "CT[C>T]CC",
      "CT[C>T]CG", "CT[C>T]CT", "CT[C>T]GA", "CT[C>T]GC", "CT[C>T]GG",
      "CT[C>T]GT", "CT[C>T]TA", "CT[C>T]TC", "CT[C>T]TG", "CT[C>T]TT",
      "CT[T>A]AA", "CT[T>A]AC", "CT[T>A]AG", "CT[T>A]AT", "CT[T>A]CA",
      "CT[T>A]CC", "CT[T>A]CG", "CT[T>A]CT", "CT[T>A]GA", "CT[T>A]GC",
      "CT[T>A]GG", "CT[T>A]GT", "CT[T>A]TA", "CT[T>A]TC", "CT[T>A]TG",
      "CT[T>A]TT", "CT[T>C]AA", "CT[T>C]AC", "CT[T>C]AG", "CT[T>C]AT",
      "CT[T>C]CA", "CT[T>C]CC", "CT[T>C]CG", "CT[T>C]CT", "CT[T>C]GA",
      "CT[T>C]GC", "CT[T>C]GG", "CT[T>C]GT", "CT[T>C]TA", "CT[T>C]TC",
      "CT[T>C]TG", "CT[T>C]TT", "CT[T>G]AA", "CT[T>G]AC", "CT[T>G]AG",
      "CT[T>G]AT", "CT[T>G]CA", "CT[T>G]CC", "CT[T>G]CG", "CT[T>G]CT",
      "CT[T>G]GA", "CT[T>G]GC", "CT[T>G]GG", "CT[T>G]GT", "CT[T>G]TA",
      "CT[T>G]TC", "CT[T>G]TG", "CT[T>G]TT", "GA[C>A]AA", "GA[C>A]AC",
      "GA[C>A]AG", "GA[C>A]AT", "GA[C>A]CA", "GA[C>A]CC", "GA[C>A]CG",
      "GA[C>A]CT", "GA[C>A]GA", "GA[C>A]GC", "GA[C>A]GG", "GA[C>A]GT",
      "GA[C>A]TA", "GA[C>A]TC", "GA[C>A]TG", "GA[C>A]TT", "GA[C>G]AA",
      "GA[C>G]AC", "GA[C>G]AG", "GA[C>G]AT", "GA[C>G]CA", "GA[C>G]CC",
      "GA[C>G]CG", "GA[C>G]CT", "GA[C>G]GA", "GA[C>G]GC", "GA[C>G]GG",
      "GA[C>G]GT", "GA[C>G]TA", "GA[C>G]TC", "GA[C>G]TG", "GA[C>G]TT",
      "GA[C>T]AA", "GA[C>T]AC", "GA[C>T]AG", "GA[C>T]AT", "GA[C>T]CA",
      "GA[C>T]CC", "GA[C>T]CG", "GA[C>T]CT", "GA[C>T]GA", "GA[C>T]GC",
      "GA[C>T]GG", "GA[C>T]GT", "GA[C>T]TA", "GA[C>T]TC", "GA[C>T]TG",
      "GA[C>T]TT", "GA[T>A]AA", "GA[T>A]AC", "GA[T>A]AG", "GA[T>A]AT",
      "GA[T>A]CA", "GA[T>A]CC", "GA[T>A]CG", "GA[T>A]CT", "GA[T>A]GA",
      "GA[T>A]GC", "GA[T>A]GG", "GA[T>A]GT", "GA[T>A]TA", "GA[T>A]TC",
      "GA[T>A]TG", "GA[T>A]TT", "GA[T>C]AA", "GA[T>C]AC", "GA[T>C]AG",
      "GA[T>C]AT", "GA[T>C]CA", "GA[T>C]CC", "GA[T>C]CG", "GA[T>C]CT",
      "GA[T>C]GA", "GA[T>C]GC", "GA[T>C]GG", "GA[T>C]GT", "GA[T>C]TA",
      "GA[T>C]TC", "GA[T>C]TG", "GA[T>C]TT", "GA[T>G]AA", "GA[T>G]AC",
      "GA[T>G]AG", "GA[T>G]AT", "GA[T>G]CA", "GA[T>G]CC", "GA[T>G]CG",
      "GA[T>G]CT", "GA[T>G]GA", "GA[T>G]GC", "GA[T>G]GG", "GA[T>G]GT",
      "GA[T>G]TA", "GA[T>G]TC", "GA[T>G]TG", "GA[T>G]TT", "GC[C>A]AA",
      "GC[C>A]AC", "GC[C>A]AG", "GC[C>A]AT", "GC[C>A]CA", "GC[C>A]CC",
      "GC[C>A]CG", "GC[C>A]CT", "GC[C>A]GA", "GC[C>A]GC", "GC[C>A]GG",
      "GC[C>A]GT", "GC[C>A]TA", "GC[C>A]TC", "GC[C>A]TG", "GC[C>A]TT",
      "GC[C>G]AA", "GC[C>G]AC", "GC[C>G]AG", "GC[C>G]AT", "GC[C>G]CA",
      "GC[C>G]CC", "GC[C>G]CG", "GC[C>G]CT", "GC[C>G]GA", "GC[C>G]GC",
      "GC[C>G]GG", "GC[C>G]GT", "GC[C>G]TA", "GC[C>G]TC", "GC[C>G]TG",
      "GC[C>G]TT", "GC[C>T]AA", "GC[C>T]AC", "GC[C>T]AG", "GC[C>T]AT",
      "GC[C>T]CA", "GC[C>T]CC", "GC[C>T]CG", "GC[C>T]CT", "GC[C>T]GA",
      "GC[C>T]GC", "GC[C>T]GG", "GC[C>T]GT", "GC[C>T]TA", "GC[C>T]TC",
      "GC[C>T]TG", "GC[C>T]TT", "GC[T>A]AA", "GC[T>A]AC", "GC[T>A]AG",
      "GC[T>A]AT", "GC[T>A]CA", "GC[T>A]CC", "GC[T>A]CG", "GC[T>A]CT",
      "GC[T>A]GA", "GC[T>A]GC", "GC[T>A]GG", "GC[T>A]GT", "GC[T>A]TA",
      "GC[T>A]TC", "GC[T>A]TG", "GC[T>A]TT", "GC[T>C]AA", "GC[T>C]AC",
      "GC[T>C]AG", "GC[T>C]AT", "GC[T>C]CA", "GC[T>C]CC", "GC[T>C]CG",
      "GC[T>C]CT", "GC[T>C]GA", "GC[T>C]GC", "GC[T>C]GG", "GC[T>C]GT",
      "GC[T>C]TA", "GC[T>C]TC", "GC[T>C]TG", "GC[T>C]TT", "GC[T>G]AA",
      "GC[T>G]AC", "GC[T>G]AG", "GC[T>G]AT", "GC[T>G]CA", "GC[T>G]CC",
      "GC[T>G]CG", "GC[T>G]CT", "GC[T>G]GA", "GC[T>G]GC", "GC[T>G]GG",
      "GC[T>G]GT", "GC[T>G]TA", "GC[T>G]TC", "GC[T>G]TG", "GC[T>G]TT",
      "GG[C>A]AA", "GG[C>A]AC", "GG[C>A]AG", "GG[C>A]AT", "GG[C>A]CA",
      "GG[C>A]CC", "GG[C>A]CG", "GG[C>A]CT", "GG[C>A]GA", "GG[C>A]GC",
      "GG[C>A]GG", "GG[C>A]GT", "GG[C>A]TA", "GG[C>A]TC", "GG[C>A]TG",
      "GG[C>A]TT", "GG[C>G]AA", "GG[C>G]AC", "GG[C>G]AG", "GG[C>G]AT",
      "GG[C>G]CA", "GG[C>G]CC", "GG[C>G]CG", "GG[C>G]CT", "GG[C>G]GA",
      "GG[C>G]GC", "GG[C>G]GG", "GG[C>G]GT", "GG[C>G]TA", "GG[C>G]TC",
      "GG[C>G]TG", "GG[C>G]TT", "GG[C>T]AA", "GG[C>T]AC", "GG[C>T]AG",
      "GG[C>T]AT", "GG[C>T]CA", "GG[C>T]CC", "GG[C>T]CG", "GG[C>T]CT",
      "GG[C>T]GA", "GG[C>T]GC", "GG[C>T]GG", "GG[C>T]GT", "GG[C>T]TA",
      "GG[C>T]TC", "GG[C>T]TG", "GG[C>T]TT", "GG[T>A]AA", "GG[T>A]AC",
      "GG[T>A]AG", "GG[T>A]AT", "GG[T>A]CA", "GG[T>A]CC", "GG[T>A]CG",
      "GG[T>A]CT", "GG[T>A]GA", "GG[T>A]GC", "GG[T>A]GG", "GG[T>A]GT",
      "GG[T>A]TA", "GG[T>A]TC", "GG[T>A]TG", "GG[T>A]TT", "GG[T>C]AA",
      "GG[T>C]AC", "GG[T>C]AG", "GG[T>C]AT", "GG[T>C]CA", "GG[T>C]CC",
      "GG[T>C]CG", "GG[T>C]CT", "GG[T>C]GA", "GG[T>C]GC", "GG[T>C]GG",
      "GG[T>C]GT", "GG[T>C]TA", "GG[T>C]TC", "GG[T>C]TG", "GG[T>C]TT",
      "GG[T>G]AA", "GG[T>G]AC", "GG[T>G]AG", "GG[T>G]AT", "GG[T>G]CA",
      "GG[T>G]CC", "GG[T>G]CG", "GG[T>G]CT", "GG[T>G]GA", "GG[T>G]GC",
      "GG[T>G]GG", "GG[T>G]GT", "GG[T>G]TA", "GG[T>G]TC", "GG[T>G]TG",
      "GG[T>G]TT", "GT[C>A]AA", "GT[C>A]AC", "GT[C>A]AG", "GT[C>A]AT",
      "GT[C>A]CA", "GT[C>A]CC", "GT[C>A]CG", "GT[C>A]CT", "GT[C>A]GA",
      "GT[C>A]GC", "GT[C>A]GG", "GT[C>A]GT", "GT[C>A]TA", "GT[C>A]TC",
      "GT[C>A]TG", "GT[C>A]TT", "GT[C>G]AA", "GT[C>G]AC", "GT[C>G]AG",
      "GT[C>G]AT", "GT[C>G]CA", "GT[C>G]CC", "GT[C>G]CG", "GT[C>G]CT",
      "GT[C>G]GA", "GT[C>G]GC", "GT[C>G]GG", "GT[C>G]GT", "GT[C>G]TA",
      "GT[C>G]TC", "GT[C>G]TG", "GT[C>G]TT", "GT[C>T]AA", "GT[C>T]AC",
      "GT[C>T]AG", "GT[C>T]AT", "GT[C>T]CA", "GT[C>T]CC", "GT[C>T]CG",
      "GT[C>T]CT", "GT[C>T]GA", "GT[C>T]GC", "GT[C>T]GG", "GT[C>T]GT",
      "GT[C>T]TA", "GT[C>T]TC", "GT[C>T]TG", "GT[C>T]TT", "GT[T>A]AA",
      "GT[T>A]AC", "GT[T>A]AG", "GT[T>A]AT", "GT[T>A]CA", "GT[T>A]CC",
      "GT[T>A]CG", "GT[T>A]CT", "GT[T>A]GA", "GT[T>A]GC", "GT[T>A]GG",
      "GT[T>A]GT", "GT[T>A]TA", "GT[T>A]TC", "GT[T>A]TG", "GT[T>A]TT",
      "GT[T>C]AA", "GT[T>C]AC", "GT[T>C]AG", "GT[T>C]AT", "GT[T>C]CA",
      "GT[T>C]CC", "GT[T>C]CG", "GT[T>C]CT", "GT[T>C]GA", "GT[T>C]GC",
      "GT[T>C]GG", "GT[T>C]GT", "GT[T>C]TA", "GT[T>C]TC", "GT[T>C]TG",
      "GT[T>C]TT", "GT[T>G]AA", "GT[T>G]AC", "GT[T>G]AG", "GT[T>G]AT",
      "GT[T>G]CA", "GT[T>G]CC", "GT[T>G]CG", "GT[T>G]CT", "GT[T>G]GA",
      "GT[T>G]GC", "GT[T>G]GG", "GT[T>G]GT", "GT[T>G]TA", "GT[T>G]TC",
      "GT[T>G]TG", "GT[T>G]TT", "TA[C>A]AA", "TA[C>A]AC", "TA[C>A]AG",
      "TA[C>A]AT", "TA[C>A]CA", "TA[C>A]CC", "TA[C>A]CG", "TA[C>A]CT",
      "TA[C>A]GA", "TA[C>A]GC", "TA[C>A]GG", "TA[C>A]GT", "TA[C>A]TA",
      "TA[C>A]TC", "TA[C>A]TG", "TA[C>A]TT", "TA[C>G]AA", "TA[C>G]AC",
      "TA[C>G]AG", "TA[C>G]AT", "TA[C>G]CA", "TA[C>G]CC", "TA[C>G]CG",
      "TA[C>G]CT", "TA[C>G]GA", "TA[C>G]GC", "TA[C>G]GG", "TA[C>G]GT",
      "TA[C>G]TA", "TA[C>G]TC", "TA[C>G]TG", "TA[C>G]TT", "TA[C>T]AA",
      "TA[C>T]AC", "TA[C>T]AG", "TA[C>T]AT", "TA[C>T]CA", "TA[C>T]CC",
      "TA[C>T]CG", "TA[C>T]CT", "TA[C>T]GA", "TA[C>T]GC", "TA[C>T]GG",
      "TA[C>T]GT", "TA[C>T]TA", "TA[C>T]TC", "TA[C>T]TG", "TA[C>T]TT",
      "TA[T>A]AA", "TA[T>A]AC", "TA[T>A]AG", "TA[T>A]AT", "TA[T>A]CA",
      "TA[T>A]CC", "TA[T>A]CG", "TA[T>A]CT", "TA[T>A]GA", "TA[T>A]GC",
      "TA[T>A]GG", "TA[T>A]GT", "TA[T>A]TA", "TA[T>A]TC", "TA[T>A]TG",
      "TA[T>A]TT", "TA[T>C]AA", "TA[T>C]AC", "TA[T>C]AG", "TA[T>C]AT",
      "TA[T>C]CA", "TA[T>C]CC", "TA[T>C]CG", "TA[T>C]CT", "TA[T>C]GA",
      "TA[T>C]GC", "TA[T>C]GG", "TA[T>C]GT", "TA[T>C]TA", "TA[T>C]TC",
      "TA[T>C]TG", "TA[T>C]TT", "TA[T>G]AA", "TA[T>G]AC", "TA[T>G]AG",
      "TA[T>G]AT", "TA[T>G]CA", "TA[T>G]CC", "TA[T>G]CG", "TA[T>G]CT",
      "TA[T>G]GA", "TA[T>G]GC", "TA[T>G]GG", "TA[T>G]GT", "TA[T>G]TA",
      "TA[T>G]TC", "TA[T>G]TG", "TA[T>G]TT", "TC[C>A]AA", "TC[C>A]AC",
      "TC[C>A]AG", "TC[C>A]AT", "TC[C>A]CA", "TC[C>A]CC", "TC[C>A]CG",
      "TC[C>A]CT", "TC[C>A]GA", "TC[C>A]GC", "TC[C>A]GG", "TC[C>A]GT",
      "TC[C>A]TA", "TC[C>A]TC", "TC[C>A]TG", "TC[C>A]TT", "TC[C>G]AA",
      "TC[C>G]AC", "TC[C>G]AG", "TC[C>G]AT", "TC[C>G]CA", "TC[C>G]CC",
      "TC[C>G]CG", "TC[C>G]CT", "TC[C>G]GA", "TC[C>G]GC", "TC[C>G]GG",
      "TC[C>G]GT", "TC[C>G]TA", "TC[C>G]TC", "TC[C>G]TG", "TC[C>G]TT",
      "TC[C>T]AA", "TC[C>T]AC", "TC[C>T]AG", "TC[C>T]AT", "TC[C>T]CA",
      "TC[C>T]CC", "TC[C>T]CG", "TC[C>T]CT", "TC[C>T]GA", "TC[C>T]GC",
      "TC[C>T]GG", "TC[C>T]GT", "TC[C>T]TA", "TC[C>T]TC", "TC[C>T]TG",
      "TC[C>T]TT", "TC[T>A]AA", "TC[T>A]AC", "TC[T>A]AG", "TC[T>A]AT",
      "TC[T>A]CA", "TC[T>A]CC", "TC[T>A]CG", "TC[T>A]CT", "TC[T>A]GA",
      "TC[T>A]GC", "TC[T>A]GG", "TC[T>A]GT", "TC[T>A]TA", "TC[T>A]TC",
      "TC[T>A]TG", "TC[T>A]TT", "TC[T>C]AA", "TC[T>C]AC", "TC[T>C]AG",
      "TC[T>C]AT", "TC[T>C]CA", "TC[T>C]CC", "TC[T>C]CG", "TC[T>C]CT",
      "TC[T>C]GA", "TC[T>C]GC", "TC[T>C]GG", "TC[T>C]GT", "TC[T>C]TA",
      "TC[T>C]TC", "TC[T>C]TG", "TC[T>C]TT", "TC[T>G]AA", "TC[T>G]AC",
      "TC[T>G]AG", "TC[T>G]AT", "TC[T>G]CA", "TC[T>G]CC", "TC[T>G]CG",
      "TC[T>G]CT", "TC[T>G]GA", "TC[T>G]GC", "TC[T>G]GG", "TC[T>G]GT",
      "TC[T>G]TA", "TC[T>G]TC", "TC[T>G]TG", "TC[T>G]TT", "TG[C>A]AA",
      "TG[C>A]AC", "TG[C>A]AG", "TG[C>A]AT", "TG[C>A]CA", "TG[C>A]CC",
      "TG[C>A]CG", "TG[C>A]CT", "TG[C>A]GA", "TG[C>A]GC", "TG[C>A]GG",
      "TG[C>A]GT", "TG[C>A]TA", "TG[C>A]TC", "TG[C>A]TG", "TG[C>A]TT",
      "TG[C>G]AA", "TG[C>G]AC", "TG[C>G]AG", "TG[C>G]AT", "TG[C>G]CA",
      "TG[C>G]CC", "TG[C>G]CG", "TG[C>G]CT", "TG[C>G]GA", "TG[C>G]GC",
      "TG[C>G]GG", "TG[C>G]GT", "TG[C>G]TA", "TG[C>G]TC", "TG[C>G]TG",
      "TG[C>G]TT", "TG[C>T]AA", "TG[C>T]AC", "TG[C>T]AG", "TG[C>T]AT",
      "TG[C>T]CA", "TG[C>T]CC", "TG[C>T]CG", "TG[C>T]CT", "TG[C>T]GA",
      "TG[C>T]GC", "TG[C>T]GG", "TG[C>T]GT", "TG[C>T]TA", "TG[C>T]TC",
      "TG[C>T]TG", "TG[C>T]TT", "TG[T>A]AA", "TG[T>A]AC", "TG[T>A]AG",
      "TG[T>A]AT", "TG[T>A]CA", "TG[T>A]CC", "TG[T>A]CG", "TG[T>A]CT",
      "TG[T>A]GA", "TG[T>A]GC", "TG[T>A]GG", "TG[T>A]GT", "TG[T>A]TA",
      "TG[T>A]TC", "TG[T>A]TG", "TG[T>A]TT", "TG[T>C]AA", "TG[T>C]AC",
      "TG[T>C]AG", "TG[T>C]AT", "TG[T>C]CA", "TG[T>C]CC", "TG[T>C]CG",
      "TG[T>C]CT", "TG[T>C]GA", "TG[T>C]GC", "TG[T>C]GG", "TG[T>C]GT",
      "TG[T>C]TA", "TG[T>C]TC", "TG[T>C]TG", "TG[T>C]TT", "TG[T>G]AA",
      "TG[T>G]AC", "TG[T>G]AG", "TG[T>G]AT", "TG[T>G]CA", "TG[T>G]CC",
      "TG[T>G]CG", "TG[T>G]CT", "TG[T>G]GA", "TG[T>G]GC", "TG[T>G]GG",
      "TG[T>G]GT", "TG[T>G]TA", "TG[T>G]TC", "TG[T>G]TG", "TG[T>G]TT",
      "TT[C>A]AA", "TT[C>A]AC", "TT[C>A]AG", "TT[C>A]AT", "TT[C>A]CA",
      "TT[C>A]CC", "TT[C>A]CG", "TT[C>A]CT", "TT[C>A]GA", "TT[C>A]GC",
      "TT[C>A]GG", "TT[C>A]GT", "TT[C>A]TA", "TT[C>A]TC", "TT[C>A]TG",
      "TT[C>A]TT", "TT[C>G]AA", "TT[C>G]AC", "TT[C>G]AG", "TT[C>G]AT",
      "TT[C>G]CA", "TT[C>G]CC", "TT[C>G]CG", "TT[C>G]CT", "TT[C>G]GA",
      "TT[C>G]GC", "TT[C>G]GG", "TT[C>G]GT", "TT[C>G]TA", "TT[C>G]TC",
      "TT[C>G]TG", "TT[C>G]TT", "TT[C>T]AA", "TT[C>T]AC", "TT[C>T]AG",
      "TT[C>T]AT", "TT[C>T]CA", "TT[C>T]CC", "TT[C>T]CG", "TT[C>T]CT",
      "TT[C>T]GA", "TT[C>T]GC", "TT[C>T]GG", "TT[C>T]GT", "TT[C>T]TA",
      "TT[C>T]TC", "TT[C>T]TG", "TT[C>T]TT", "TT[T>A]AA", "TT[T>A]AC",
      "TT[T>A]AG", "TT[T>A]AT", "TT[T>A]CA", "TT[T>A]CC", "TT[T>A]CG",
      "TT[T>A]CT", "TT[T>A]GA", "TT[T>A]GC", "TT[T>A]GG", "TT[T>A]GT",
      "TT[T>A]TA", "TT[T>A]TC", "TT[T>A]TG", "TT[T>A]TT", "TT[T>C]AA",
      "TT[T>C]AC", "TT[T>C]AG", "TT[T>C]AT", "TT[T>C]CA", "TT[T>C]CC",
      "TT[T>C]CG", "TT[T>C]CT", "TT[T>C]GA", "TT[T>C]GC", "TT[T>C]GG",
      "TT[T>C]GT", "TT[T>C]TA", "TT[T>C]TC", "TT[T>C]TG", "TT[T>C]TT",
      "TT[T>G]AA", "TT[T>G]AC", "TT[T>G]AG", "TT[T>G]AT", "TT[T>G]CA",
      "TT[T>G]CC", "TT[T>G]CG", "TT[T>G]CT", "TT[T>G]GA", "TT[T>G]GC",
      "TT[T>G]GG", "TT[T>G]GT", "TT[T>G]TA", "TT[T>G]TC", "TT[T>G]TG",
      "TT[T>G]TT"
    ),
    "ID83" = c(
      "1:Del:C:0", "1:Del:C:1", "1:Del:C:2", "1:Del:C:3", "1:Del:C:4",
      "1:Del:C:5", "1:Del:T:0", "1:Del:T:1", "1:Del:T:2", "1:Del:T:3",
      "1:Del:T:4", "1:Del:T:5", "1:Ins:C:0", "1:Ins:C:1", "1:Ins:C:2",
      "1:Ins:C:3", "1:Ins:C:4", "1:Ins:C:5", "1:Ins:T:0", "1:Ins:T:1",
      "1:Ins:T:2", "1:Ins:T:3", "1:Ins:T:4", "1:Ins:T:5", "2:Del:M:1",
      "2:Del:R:0", "2:Del:R:1", "2:Del:R:2", "2:Del:R:3", "2:Del:R:4",
      "2:Del:R:5", "2:Ins:R:0", "2:Ins:R:1", "2:Ins:R:2", "2:Ins:R:3",
      "2:Ins:R:4", "2:Ins:R:5", "3:Del:M:1", "3:Del:M:2", "3:Del:R:0",
      "3:Del:R:1", "3:Del:R:2", "3:Del:R:3", "3:Del:R:4", "3:Del:R:5",
      "3:Ins:R:0", "3:Ins:R:1", "3:Ins:R:2", "3:Ins:R:3", "3:Ins:R:4",
      "3:Ins:R:5", "4:Del:M:1", "4:Del:M:2", "4:Del:M:3", "4:Del:R:0",
      "4:Del:R:1", "4:Del:R:2", "4:Del:R:3", "4:Del:R:4", "4:Del:R:5",
      "4:Ins:R:0", "4:Ins:R:1", "4:Ins:R:2", "4:Ins:R:3", "4:Ins:R:4",
      "4:Ins:R:5", "5:Del:M:1", "5:Del:M:2", "5:Del:M:3", "5:Del:M:4",
      "5:Del:M:5", "5:Del:R:0", "5:Del:R:1", "5:Del:R:2", "5:Del:R:3",
      "5:Del:R:4", "5:Del:R:5", "5:Ins:R:0", "5:Ins:R:1", "5:Ins:R:2",
      "5:Ins:R:3", "5:Ins:R:4", "5:Ins:R:5"
    ),
    "ID28" = c(
      "1:Del:C:0", "1:Del:C:1", "1:Del:C:2", "1:Del:C:3", "1:Del:C:4",
      "1:Del:C:5", "1:Del:T:0", "1:Del:T:1", "1:Del:T:2", "1:Del:T:3",
      "1:Del:T:4", "1:Del:T:5", "1:Ins:C:0", "1:Ins:C:1", "1:Ins:C:2",
      "1:Ins:C:3", "1:Ins:C:4", "1:Ins:C:5", "1:Ins:T:0", "1:Ins:T:1",
      "1:Ins:T:2", "1:Ins:T:3", "1:Ins:T:4", "1:Ins:T:5", "complex",
      "long_Del", "long_Ins", "MH"
    ),
    "DBS78" = c(
      "AC>CA", "AC>CG", "AC>CT", "AC>GA", "AC>GG", "AC>GT", "AC>TA",
      "AC>TG", "AC>TT", "AT>CA", "AT>CC", "AT>CG", "AT>GA", "AT>GC",
      "AT>TA", "CC>AA", "CC>AG", "CC>AT", "CC>GA", "CC>GG", "CC>GT",
      "CC>TA", "CC>TG", "CC>TT", "CG>AT", "CG>GC", "CG>GT", "CG>TA",
      "CG>TC", "CG>TT", "CT>AA", "CT>AC", "CT>AG", "CT>GA", "CT>GC",
      "CT>GG", "CT>TA", "CT>TC", "CT>TG", "GC>AA", "GC>AG", "GC>AT",
      "GC>CA", "GC>CG", "GC>TA", "TA>AT", "TA>CG", "TA>CT", "TA>GC",
      "TA>GG", "TA>GT", "TC>AA", "TC>AG", "TC>AT", "TC>CA", "TC>CG",
      "TC>CT", "TC>GA", "TC>GG", "TC>GT", "TG>AA", "TG>AC", "TG>AT",
      "TG>CA", "TG>CC", "TG>CT", "TG>GA", "TG>GC", "TG>GT", "TT>AA",
      "TT>AC", "TT>AG", "TT>CA", "TT>CC", "TT>CG", "TT>GA", "TT>GC",
      "TT>GG"
    ),
    "DBS1248" = c(
      "A[AC>CA]A", "A[AC>CA]C", "A[AC>CA]G", "A[AC>CA]T", "A[AC>CG]A",
      "A[AC>CG]C", "A[AC>CG]G", "A[AC>CG]T", "A[AC>CT]A", "A[AC>CT]C",
      "A[AC>CT]G", "A[AC>CT]T", "A[AC>GA]A", "A[AC>GA]C", "A[AC>GA]G",
      "A[AC>GA]T", "A[AC>GG]A", "A[AC>GG]C", "A[AC>GG]G", "A[AC>GG]T",
      "A[AC>GT]A", "A[AC>GT]C", "A[AC>GT]G", "A[AC>GT]T", "A[AC>TA]A",
      "A[AC>TA]C", "A[AC>TA]G", "A[AC>TA]T", "A[AC>TG]A", "A[AC>TG]C",
      "A[AC>TG]G", "A[AC>TG]T", "A[AC>TT]A", "A[AC>TT]C", "A[AC>TT]G",
      "A[AC>TT]T", "A[AT>CA]A", "A[AT>CA]C", "A[AT>CA]G", "A[AT>CA]T",
      "A[AT>CC]A", "A[AT>CC]C", "A[AT>CC]G", "A[AT>CC]T", "A[AT>CG]A",
      "A[AT>CG]C", "A[AT>CG]G", "A[AT>CG]T", "A[AT>GA]A", "A[AT>GA]C",
      "A[AT>GA]G", "A[AT>GA]T", "A[AT>GC]A", "A[AT>GC]C", "A[AT>GC]G",
      "A[AT>GC]T", "A[AT>TA]A", "A[AT>TA]C", "A[AT>TA]G", "A[AT>TA]T",
      "A[CC>AA]A", "A[CC>AA]C", "A[CC>AA]G", "A[CC>AA]T", "A[CC>AG]A",
      "A[CC>AG]C", "A[CC>AG]G", "A[CC>AG]T", "A[CC>AT]A", "A[CC>AT]C",
      "A[CC>AT]G", "A[CC>AT]T", "A[CC>GA]A", "A[CC>GA]C", "A[CC>GA]G",
      "A[CC>GA]T", "A[CC>GG]A", "A[CC>GG]C", "A[CC>GG]G", "A[CC>GG]T",
      "A[CC>GT]A", "A[CC>GT]C", "A[CC>GT]G", "A[CC>GT]T", "A[CC>TA]A",
      "A[CC>TA]C", "A[CC>TA]G", "A[CC>TA]T", "A[CC>TG]A", "A[CC>TG]C",
      "A[CC>TG]G", "A[CC>TG]T", "A[CC>TT]A", "A[CC>TT]C", "A[CC>TT]G",
      "A[CC>TT]T", "A[CG>AT]A", "A[CG>AT]C", "A[CG>AT]G", "A[CG>AT]T",
      "A[CG>GC]A", "A[CG>GC]C", "A[CG>GC]G", "A[CG>GC]T", "A[CG>GT]A",
      "A[CG>GT]C", "A[CG>GT]G", "A[CG>GT]T", "A[CG>TA]A", "A[CG>TA]C",
      "A[CG>TA]G", "A[CG>TA]T", "A[CG>TC]A", "A[CG>TC]C", "A[CG>TC]G",
      "A[CG>TC]T", "A[CG>TT]A", "A[CG>TT]C", "A[CG>TT]G", "A[CG>TT]T",
      "A[CT>AA]A", "A[CT>AA]C", "A[CT>AA]G", "A[CT>AA]T", "A[CT>AC]A",
      "A[CT>AC]C", "A[CT>AC]G", "A[CT>AC]T", "A[CT>AG]A", "A[CT>AG]C",
      "A[CT>AG]G", "A[CT>AG]T", "A[CT>GA]A", "A[CT>GA]C", "A[CT>GA]G",
      "A[CT>GA]T", "A[CT>GC]A", "A[CT>GC]C", "A[CT>GC]G", "A[CT>GC]T",
      "A[CT>GG]A", "A[CT>GG]C", "A[CT>GG]G", "A[CT>GG]T", "A[CT>TA]A",
      "A[CT>TA]C", "A[CT>TA]G", "A[CT>TA]T", "A[CT>TC]A", "A[CT>TC]C",
      "A[CT>TC]G", "A[CT>TC]T", "A[CT>TG]A", "A[CT>TG]C", "A[CT>TG]G",
      "A[CT>TG]T", "A[GC>AA]A", "A[GC>AA]C", "A[GC>AA]G", "A[GC>AA]T",
      "A[GC>AG]A", "A[GC>AG]C", "A[GC>AG]G", "A[GC>AG]T", "A[GC>AT]A",
      "A[GC>AT]C", "A[GC>AT]G", "A[GC>AT]T", "A[GC>CA]A", "A[GC>CA]C",
      "A[GC>CA]G", "A[GC>CA]T", "A[GC>CG]A", "A[GC>CG]C", "A[GC>CG]G",
      "A[GC>CG]T", "A[GC>TA]A", "A[GC>TA]C", "A[GC>TA]G", "A[GC>TA]T",
      "A[TA>AT]A", "A[TA>AT]C", "A[TA>AT]G", "A[TA>AT]T", "A[TA>CG]A",
      "A[TA>CG]C", "A[TA>CG]G", "A[TA>CG]T", "A[TA>CT]A", "A[TA>CT]C",
      "A[TA>CT]G", "A[TA>CT]T", "A[TA>GC]A", "A[TA>GC]C", "A[TA>GC]G",
      "A[TA>GC]T", "A[TA>GG]A", "A[TA>GG]C", "A[TA>GG]G", "A[TA>GG]T",
      "A[TA>GT]A", "A[TA>GT]C", "A[TA>GT]G", "A[TA>GT]T", "A[TC>AA]A",
      "A[TC>AA]C", "A[TC>AA]G", "A[TC>AA]T", "A[TC>AG]A", "A[TC>AG]C",
      "A[TC>AG]G", "A[TC>AG]T", "A[TC>AT]A", "A[TC>AT]C", "A[TC>AT]G",
      "A[TC>AT]T", "A[TC>CA]A", "A[TC>CA]C", "A[TC>CA]G", "A[TC>CA]T",
      "A[TC>CG]A", "A[TC>CG]C", "A[TC>CG]G", "A[TC>CG]T", "A[TC>CT]A",
      "A[TC>CT]C", "A[TC>CT]G", "A[TC>CT]T", "A[TC>GA]A", "A[TC>GA]C",
      "A[TC>GA]G", "A[TC>GA]T", "A[TC>GG]A", "A[TC>GG]C", "A[TC>GG]G",
      "A[TC>GG]T", "A[TC>GT]A", "A[TC>GT]C", "A[TC>GT]G", "A[TC>GT]T",
      "A[TG>AA]A", "A[TG>AA]C", "A[TG>AA]G", "A[TG>AA]T", "A[TG>AC]A",
      "A[TG>AC]C", "A[TG>AC]G", "A[TG>AC]T", "A[TG>AT]A", "A[TG>AT]C",
      "A[TG>AT]G", "A[TG>AT]T", "A[TG>CA]A", "A[TG>CA]C", "A[TG>CA]G",
      "A[TG>CA]T", "A[TG>CC]A", "A[TG>CC]C", "A[TG>CC]G", "A[TG>CC]T",
      "A[TG>CT]A", "A[TG>CT]C", "A[TG>CT]G", "A[TG>CT]T", "A[TG>GA]A",
      "A[TG>GA]C", "A[TG>GA]G", "A[TG>GA]T", "A[TG>GC]A", "A[TG>GC]C",
      "A[TG>GC]G", "A[TG>GC]T", "A[TG>GT]A", "A[TG>GT]C", "A[TG>GT]G",
      "A[TG>GT]T", "A[TT>AA]A", "A[TT>AA]C", "A[TT>AA]G", "A[TT>AA]T",
      "A[TT>AC]A", "A[TT>AC]C", "A[TT>AC]G", "A[TT>AC]T", "A[TT>AG]A",
      "A[TT>AG]C", "A[TT>AG]G", "A[TT>AG]T", "A[TT>CA]A", "A[TT>CA]C",
      "A[TT>CA]G", "A[TT>CA]T", "A[TT>CC]A", "A[TT>CC]C", "A[TT>CC]G",
      "A[TT>CC]T", "A[TT>CG]A", "A[TT>CG]C", "A[TT>CG]G", "A[TT>CG]T",
      "A[TT>GA]A", "A[TT>GA]C", "A[TT>GA]G", "A[TT>GA]T", "A[TT>GC]A",
      "A[TT>GC]C", "A[TT>GC]G", "A[TT>GC]T", "A[TT>GG]A", "A[TT>GG]C",
      "A[TT>GG]G", "A[TT>GG]T", "C[AC>CA]A", "C[AC>CA]C", "C[AC>CA]G",
      "C[AC>CA]T", "C[AC>CG]A", "C[AC>CG]C", "C[AC>CG]G", "C[AC>CG]T",
      "C[AC>CT]A", "C[AC>CT]C", "C[AC>CT]G", "C[AC>CT]T", "C[AC>GA]A",
      "C[AC>GA]C", "C[AC>GA]G", "C[AC>GA]T", "C[AC>GG]A", "C[AC>GG]C",
      "C[AC>GG]G", "C[AC>GG]T", "C[AC>GT]A", "C[AC>GT]C", "C[AC>GT]G",
      "C[AC>GT]T", "C[AC>TA]A", "C[AC>TA]C", "C[AC>TA]G", "C[AC>TA]T",
      "C[AC>TG]A", "C[AC>TG]C", "C[AC>TG]G", "C[AC>TG]T", "C[AC>TT]A",
      "C[AC>TT]C", "C[AC>TT]G", "C[AC>TT]T", "C[AT>CA]A", "C[AT>CA]C",
      "C[AT>CA]G", "C[AT>CA]T", "C[AT>CC]A", "C[AT>CC]C", "C[AT>CC]G",
      "C[AT>CC]T", "C[AT>CG]A", "C[AT>CG]C", "C[AT>CG]G", "C[AT>CG]T",
      "C[AT>GA]A", "C[AT>GA]C", "C[AT>GA]G", "C[AT>GA]T", "C[AT>GC]A",
      "C[AT>GC]C", "C[AT>GC]G", "C[AT>GC]T", "C[AT>TA]A", "C[AT>TA]C",
      "C[AT>TA]G", "C[AT>TA]T", "C[CC>AA]A", "C[CC>AA]C", "C[CC>AA]G",
      "C[CC>AA]T", "C[CC>AG]A", "C[CC>AG]C", "C[CC>AG]G", "C[CC>AG]T",
      "C[CC>AT]A", "C[CC>AT]C", "C[CC>AT]G", "C[CC>AT]T", "C[CC>GA]A",
      "C[CC>GA]C", "C[CC>GA]G", "C[CC>GA]T", "C[CC>GG]A", "C[CC>GG]C",
      "C[CC>GG]G", "C[CC>GG]T", "C[CC>GT]A", "C[CC>GT]C", "C[CC>GT]G",
      "C[CC>GT]T", "C[CC>TA]A", "C[CC>TA]C", "C[CC>TA]G", "C[CC>TA]T",
      "C[CC>TG]A", "C[CC>TG]C", "C[CC>TG]G", "C[CC>TG]T", "C[CC>TT]A",
      "C[CC>TT]C", "C[CC>TT]G", "C[CC>TT]T", "C[CG>AT]A", "C[CG>AT]C",
      "C[CG>AT]G", "C[CG>AT]T", "C[CG>GC]A", "C[CG>GC]C", "C[CG>GC]G",
      "C[CG>GC]T", "C[CG>GT]A", "C[CG>GT]C", "C[CG>GT]G", "C[CG>GT]T",
      "C[CG>TA]A", "C[CG>TA]C", "C[CG>TA]G", "C[CG>TA]T", "C[CG>TC]A",
      "C[CG>TC]C", "C[CG>TC]G", "C[CG>TC]T", "C[CG>TT]A", "C[CG>TT]C",
      "C[CG>TT]G", "C[CG>TT]T", "C[CT>AA]A", "C[CT>AA]C", "C[CT>AA]G",
      "C[CT>AA]T", "C[CT>AC]A", "C[CT>AC]C", "C[CT>AC]G", "C[CT>AC]T",
      "C[CT>AG]A", "C[CT>AG]C", "C[CT>AG]G", "C[CT>AG]T", "C[CT>GA]A",
      "C[CT>GA]C", "C[CT>GA]G", "C[CT>GA]T", "C[CT>GC]A", "C[CT>GC]C",
      "C[CT>GC]G", "C[CT>GC]T", "C[CT>GG]A", "C[CT>GG]C", "C[CT>GG]G",
      "C[CT>GG]T", "C[CT>TA]A", "C[CT>TA]C", "C[CT>TA]G", "C[CT>TA]T",
      "C[CT>TC]A", "C[CT>TC]C", "C[CT>TC]G", "C[CT>TC]T", "C[CT>TG]A",
      "C[CT>TG]C", "C[CT>TG]G", "C[CT>TG]T", "C[GC>AA]A", "C[GC>AA]C",
      "C[GC>AA]G", "C[GC>AA]T", "C[GC>AG]A", "C[GC>AG]C", "C[GC>AG]G",
      "C[GC>AG]T", "C[GC>AT]A", "C[GC>AT]C", "C[GC>AT]G", "C[GC>AT]T",
      "C[GC>CA]A", "C[GC>CA]C", "C[GC>CA]G", "C[GC>CA]T", "C[GC>CG]A",
      "C[GC>CG]C", "C[GC>CG]G", "C[GC>CG]T", "C[GC>TA]A", "C[GC>TA]C",
      "C[GC>TA]G", "C[GC>TA]T", "C[TA>AT]A", "C[TA>AT]C", "C[TA>AT]G",
      "C[TA>AT]T", "C[TA>CG]A", "C[TA>CG]C", "C[TA>CG]G", "C[TA>CG]T",
      "C[TA>CT]A", "C[TA>CT]C", "C[TA>CT]G", "C[TA>CT]T", "C[TA>GC]A",
      "C[TA>GC]C", "C[TA>GC]G", "C[TA>GC]T", "C[TA>GG]A", "C[TA>GG]C",
      "C[TA>GG]G", "C[TA>GG]T", "C[TA>GT]A", "C[TA>GT]C", "C[TA>GT]G",
      "C[TA>GT]T", "C[TC>AA]A", "C[TC>AA]C", "C[TC>AA]G", "C[TC>AA]T",
      "C[TC>AG]A", "C[TC>AG]C", "C[TC>AG]G", "C[TC>AG]T", "C[TC>AT]A",
      "C[TC>AT]C", "C[TC>AT]G", "C[TC>AT]T", "C[TC>CA]A", "C[TC>CA]C",
      "C[TC>CA]G", "C[TC>CA]T", "C[TC>CG]A", "C[TC>CG]C", "C[TC>CG]G",
      "C[TC>CG]T", "C[TC>CT]A", "C[TC>CT]C", "C[TC>CT]G", "C[TC>CT]T",
      "C[TC>GA]A", "C[TC>GA]C", "C[TC>GA]G", "C[TC>GA]T", "C[TC>GG]A",
      "C[TC>GG]C", "C[TC>GG]G", "C[TC>GG]T", "C[TC>GT]A", "C[TC>GT]C",
      "C[TC>GT]G", "C[TC>GT]T", "C[TG>AA]A", "C[TG>AA]C", "C[TG>AA]G",
      "C[TG>AA]T", "C[TG>AC]A", "C[TG>AC]C", "C[TG>AC]G", "C[TG>AC]T",
      "C[TG>AT]A", "C[TG>AT]C", "C[TG>AT]G", "C[TG>AT]T", "C[TG>CA]A",
      "C[TG>CA]C", "C[TG>CA]G", "C[TG>CA]T", "C[TG>CC]A", "C[TG>CC]C",
      "C[TG>CC]G", "C[TG>CC]T", "C[TG>CT]A", "C[TG>CT]C", "C[TG>CT]G",
      "C[TG>CT]T", "C[TG>GA]A", "C[TG>GA]C", "C[TG>GA]G", "C[TG>GA]T",
      "C[TG>GC]A", "C[TG>GC]C", "C[TG>GC]G", "C[TG>GC]T", "C[TG>GT]A",
      "C[TG>GT]C", "C[TG>GT]G", "C[TG>GT]T", "C[TT>AA]A", "C[TT>AA]C",
      "C[TT>AA]G", "C[TT>AA]T", "C[TT>AC]A", "C[TT>AC]C", "C[TT>AC]G",
      "C[TT>AC]T", "C[TT>AG]A", "C[TT>AG]C", "C[TT>AG]G", "C[TT>AG]T",
      "C[TT>CA]A", "C[TT>CA]C", "C[TT>CA]G", "C[TT>CA]T", "C[TT>CC]A",
      "C[TT>CC]C", "C[TT>CC]G", "C[TT>CC]T", "C[TT>CG]A", "C[TT>CG]C",
      "C[TT>CG]G", "C[TT>CG]T", "C[TT>GA]A", "C[TT>GA]C", "C[TT>GA]G",
      "C[TT>GA]T", "C[TT>GC]A", "C[TT>GC]C", "C[TT>GC]G", "C[TT>GC]T",
      "C[TT>GG]A", "C[TT>GG]C", "C[TT>GG]G", "C[TT>GG]T", "G[AC>CA]A",
      "G[AC>CA]C", "G[AC>CA]G", "G[AC>CA]T", "G[AC>CG]A", "G[AC>CG]C",
      "G[AC>CG]G", "G[AC>CG]T", "G[AC>CT]A", "G[AC>CT]C", "G[AC>CT]G",
      "G[AC>CT]T", "G[AC>GA]A", "G[AC>GA]C", "G[AC>GA]G", "G[AC>GA]T",
      "G[AC>GG]A", "G[AC>GG]C", "G[AC>GG]G", "G[AC>GG]T", "G[AC>GT]A",
      "G[AC>GT]C", "G[AC>GT]G", "G[AC>GT]T", "G[AC>TA]A", "G[AC>TA]C",
      "G[AC>TA]G", "G[AC>TA]T", "G[AC>TG]A", "G[AC>TG]C", "G[AC>TG]G",
      "G[AC>TG]T", "G[AC>TT]A", "G[AC>TT]C", "G[AC>TT]G", "G[AC>TT]T",
      "G[AT>CA]A", "G[AT>CA]C", "G[AT>CA]G", "G[AT>CA]T", "G[AT>CC]A",
      "G[AT>CC]C", "G[AT>CC]G", "G[AT>CC]T", "G[AT>CG]A", "G[AT>CG]C",
      "G[AT>CG]G", "G[AT>CG]T", "G[AT>GA]A", "G[AT>GA]C", "G[AT>GA]G",
      "G[AT>GA]T", "G[AT>GC]A", "G[AT>GC]C", "G[AT>GC]G", "G[AT>GC]T",
      "G[AT>TA]A", "G[AT>TA]C", "G[AT>TA]G", "G[AT>TA]T", "G[CC>AA]A",
      "G[CC>AA]C", "G[CC>AA]G", "G[CC>AA]T", "G[CC>AG]A", "G[CC>AG]C",
      "G[CC>AG]G", "G[CC>AG]T", "G[CC>AT]A", "G[CC>AT]C", "G[CC>AT]G",
      "G[CC>AT]T", "G[CC>GA]A", "G[CC>GA]C", "G[CC>GA]G", "G[CC>GA]T",
      "G[CC>GG]A", "G[CC>GG]C", "G[CC>GG]G", "G[CC>GG]T", "G[CC>GT]A",
      "G[CC>GT]C", "G[CC>GT]G", "G[CC>GT]T", "G[CC>TA]A", "G[CC>TA]C",
      "G[CC>TA]G", "G[CC>TA]T", "G[CC>TG]A", "G[CC>TG]C", "G[CC>TG]G",
      "G[CC>TG]T", "G[CC>TT]A", "G[CC>TT]C", "G[CC>TT]G", "G[CC>TT]T",
      "G[CG>AT]A", "G[CG>AT]C", "G[CG>AT]G", "G[CG>AT]T", "G[CG>GC]A",
      "G[CG>GC]C", "G[CG>GC]G", "G[CG>GC]T", "G[CG>GT]A", "G[CG>GT]C",
      "G[CG>GT]G", "G[CG>GT]T", "G[CG>TA]A", "G[CG>TA]C", "G[CG>TA]G",
      "G[CG>TA]T", "G[CG>TC]A", "G[CG>TC]C", "G[CG>TC]G", "G[CG>TC]T",
      "G[CG>TT]A", "G[CG>TT]C", "G[CG>TT]G", "G[CG>TT]T", "G[CT>AA]A",
      "G[CT>AA]C", "G[CT>AA]G", "G[CT>AA]T", "G[CT>AC]A", "G[CT>AC]C",
      "G[CT>AC]G", "G[CT>AC]T", "G[CT>AG]A", "G[CT>AG]C", "G[CT>AG]G",
      "G[CT>AG]T", "G[CT>GA]A", "G[CT>GA]C", "G[CT>GA]G", "G[CT>GA]T",
      "G[CT>GC]A", "G[CT>GC]C", "G[CT>GC]G", "G[CT>GC]T", "G[CT>GG]A",
      "G[CT>GG]C", "G[CT>GG]G", "G[CT>GG]T", "G[CT>TA]A", "G[CT>TA]C",
      "G[CT>TA]G", "G[CT>TA]T", "G[CT>TC]A", "G[CT>TC]C", "G[CT>TC]G",
      "G[CT>TC]T", "G[CT>TG]A", "G[CT>TG]C", "G[CT>TG]G", "G[CT>TG]T",
      "G[GC>AA]A", "G[GC>AA]C", "G[GC>AA]G", "G[GC>AA]T", "G[GC>AG]A",
      "G[GC>AG]C", "G[GC>AG]G", "G[GC>AG]T", "G[GC>AT]A", "G[GC>AT]C",
      "G[GC>AT]G", "G[GC>AT]T", "G[GC>CA]A", "G[GC>CA]C", "G[GC>CA]G",
      "G[GC>CA]T", "G[GC>CG]A", "G[GC>CG]C", "G[GC>CG]G", "G[GC>CG]T",
      "G[GC>TA]A", "G[GC>TA]C", "G[GC>TA]G", "G[GC>TA]T", "G[TA>AT]A",
      "G[TA>AT]C", "G[TA>AT]G", "G[TA>AT]T", "G[TA>CG]A", "G[TA>CG]C",
      "G[TA>CG]G", "G[TA>CG]T", "G[TA>CT]A", "G[TA>CT]C", "G[TA>CT]G",
      "G[TA>CT]T", "G[TA>GC]A", "G[TA>GC]C", "G[TA>GC]G", "G[TA>GC]T",
      "G[TA>GG]A", "G[TA>GG]C", "G[TA>GG]G", "G[TA>GG]T", "G[TA>GT]A",
      "G[TA>GT]C", "G[TA>GT]G", "G[TA>GT]T", "G[TC>AA]A", "G[TC>AA]C",
      "G[TC>AA]G", "G[TC>AA]T", "G[TC>AG]A", "G[TC>AG]C", "G[TC>AG]G",
      "G[TC>AG]T", "G[TC>AT]A", "G[TC>AT]C", "G[TC>AT]G", "G[TC>AT]T",
      "G[TC>CA]A", "G[TC>CA]C", "G[TC>CA]G", "G[TC>CA]T", "G[TC>CG]A",
      "G[TC>CG]C", "G[TC>CG]G", "G[TC>CG]T", "G[TC>CT]A", "G[TC>CT]C",
      "G[TC>CT]G", "G[TC>CT]T", "G[TC>GA]A", "G[TC>GA]C", "G[TC>GA]G",
      "G[TC>GA]T", "G[TC>GG]A", "G[TC>GG]C", "G[TC>GG]G", "G[TC>GG]T",
      "G[TC>GT]A", "G[TC>GT]C", "G[TC>GT]G", "G[TC>GT]T", "G[TG>AA]A",
      "G[TG>AA]C", "G[TG>AA]G", "G[TG>AA]T", "G[TG>AC]A", "G[TG>AC]C",
      "G[TG>AC]G", "G[TG>AC]T", "G[TG>AT]A", "G[TG>AT]C", "G[TG>AT]G",
      "G[TG>AT]T", "G[TG>CA]A", "G[TG>CA]C", "G[TG>CA]G", "G[TG>CA]T",
      "G[TG>CC]A", "G[TG>CC]C", "G[TG>CC]G", "G[TG>CC]T", "G[TG>CT]A",
      "G[TG>CT]C", "G[TG>CT]G", "G[TG>CT]T", "G[TG>GA]A", "G[TG>GA]C",
      "G[TG>GA]G", "G[TG>GA]T", "G[TG>GC]A", "G[TG>GC]C", "G[TG>GC]G",
      "G[TG>GC]T", "G[TG>GT]A", "G[TG>GT]C", "G[TG>GT]G", "G[TG>GT]T",
      "G[TT>AA]A", "G[TT>AA]C", "G[TT>AA]G", "G[TT>AA]T", "G[TT>AC]A",
      "G[TT>AC]C", "G[TT>AC]G", "G[TT>AC]T", "G[TT>AG]A", "G[TT>AG]C",
      "G[TT>AG]G", "G[TT>AG]T", "G[TT>CA]A", "G[TT>CA]C", "G[TT>CA]G",
      "G[TT>CA]T", "G[TT>CC]A", "G[TT>CC]C", "G[TT>CC]G", "G[TT>CC]T",
      "G[TT>CG]A", "G[TT>CG]C", "G[TT>CG]G", "G[TT>CG]T", "G[TT>GA]A",
      "G[TT>GA]C", "G[TT>GA]G", "G[TT>GA]T", "G[TT>GC]A", "G[TT>GC]C",
      "G[TT>GC]G", "G[TT>GC]T", "G[TT>GG]A", "G[TT>GG]C", "G[TT>GG]G",
      "G[TT>GG]T", "T[AC>CA]A", "T[AC>CA]C", "T[AC>CA]G", "T[AC>CA]T",
      "T[AC>CG]A", "T[AC>CG]C", "T[AC>CG]G", "T[AC>CG]T", "T[AC>CT]A",
      "T[AC>CT]C", "T[AC>CT]G", "T[AC>CT]T", "T[AC>GA]A", "T[AC>GA]C",
      "T[AC>GA]G", "T[AC>GA]T", "T[AC>GG]A", "T[AC>GG]C", "T[AC>GG]G",
      "T[AC>GG]T", "T[AC>GT]A", "T[AC>GT]C", "T[AC>GT]G", "T[AC>GT]T",
      "T[AC>TA]A", "T[AC>TA]C", "T[AC>TA]G", "T[AC>TA]T", "T[AC>TG]A",
      "T[AC>TG]C", "T[AC>TG]G", "T[AC>TG]T", "T[AC>TT]A", "T[AC>TT]C",
      "T[AC>TT]G", "T[AC>TT]T", "T[AT>CA]A", "T[AT>CA]C", "T[AT>CA]G",
      "T[AT>CA]T", "T[AT>CC]A", "T[AT>CC]C", "T[AT>CC]G", "T[AT>CC]T",
      "T[AT>CG]A", "T[AT>CG]C", "T[AT>CG]G", "T[AT>CG]T", "T[AT>GA]A",
      "T[AT>GA]C", "T[AT>GA]G", "T[AT>GA]T", "T[AT>GC]A", "T[AT>GC]C",
      "T[AT>GC]G", "T[AT>GC]T", "T[AT>TA]A", "T[AT>TA]C", "T[AT>TA]G",
      "T[AT>TA]T", "T[CC>AA]A", "T[CC>AA]C", "T[CC>AA]G", "T[CC>AA]T",
      "T[CC>AG]A", "T[CC>AG]C", "T[CC>AG]G", "T[CC>AG]T", "T[CC>AT]A",
      "T[CC>AT]C", "T[CC>AT]G", "T[CC>AT]T", "T[CC>GA]A", "T[CC>GA]C",
      "T[CC>GA]G", "T[CC>GA]T", "T[CC>GG]A", "T[CC>GG]C", "T[CC>GG]G",
      "T[CC>GG]T", "T[CC>GT]A", "T[CC>GT]C", "T[CC>GT]G", "T[CC>GT]T",
      "T[CC>TA]A", "T[CC>TA]C", "T[CC>TA]G", "T[CC>TA]T", "T[CC>TG]A",
      "T[CC>TG]C", "T[CC>TG]G", "T[CC>TG]T", "T[CC>TT]A", "T[CC>TT]C",
      "T[CC>TT]G", "T[CC>TT]T", "T[CG>AT]A", "T[CG>AT]C", "T[CG>AT]G",
      "T[CG>AT]T", "T[CG>GC]A", "T[CG>GC]C", "T[CG>GC]G", "T[CG>GC]T",
      "T[CG>GT]A", "T[CG>GT]C", "T[CG>GT]G", "T[CG>GT]T", "T[CG>TA]A",
      "T[CG>TA]C", "T[CG>TA]G", "T[CG>TA]T", "T[CG>TC]A", "T[CG>TC]C",
      "T[CG>TC]G", "T[CG>TC]T", "T[CG>TT]A", "T[CG>TT]C", "T[CG>TT]G",
      "T[CG>TT]T", "T[CT>AA]A", "T[CT>AA]C", "T[CT>AA]G", "T[CT>AA]T",
      "T[CT>AC]A", "T[CT>AC]C", "T[CT>AC]G", "T[CT>AC]T", "T[CT>AG]A",
      "T[CT>AG]C", "T[CT>AG]G", "T[CT>AG]T", "T[CT>GA]A", "T[CT>GA]C",
      "T[CT>GA]G", "T[CT>GA]T", "T[CT>GC]A", "T[CT>GC]C", "T[CT>GC]G",
      "T[CT>GC]T", "T[CT>GG]A", "T[CT>GG]C", "T[CT>GG]G", "T[CT>GG]T",
      "T[CT>TA]A", "T[CT>TA]C", "T[CT>TA]G", "T[CT>TA]T", "T[CT>TC]A",
      "T[CT>TC]C", "T[CT>TC]G", "T[CT>TC]T", "T[CT>TG]A", "T[CT>TG]C",
      "T[CT>TG]G", "T[CT>TG]T", "T[GC>AA]A", "T[GC>AA]C", "T[GC>AA]G",
      "T[GC>AA]T", "T[GC>AG]A", "T[GC>AG]C", "T[GC>AG]G", "T[GC>AG]T",
      "T[GC>AT]A", "T[GC>AT]C", "T[GC>AT]G", "T[GC>AT]T", "T[GC>CA]A",
      "T[GC>CA]C", "T[GC>CA]G", "T[GC>CA]T", "T[GC>CG]A", "T[GC>CG]C",
      "T[GC>CG]G", "T[GC>CG]T", "T[GC>TA]A", "T[GC>TA]C", "T[GC>TA]G",
      "T[GC>TA]T", "T[TA>AT]A", "T[TA>AT]C", "T[TA>AT]G", "T[TA>AT]T",
      "T[TA>CG]A", "T[TA>CG]C", "T[TA>CG]G", "T[TA>CG]T", "T[TA>CT]A",
      "T[TA>CT]C", "T[TA>CT]G", "T[TA>CT]T", "T[TA>GC]A", "T[TA>GC]C",
      "T[TA>GC]G", "T[TA>GC]T", "T[TA>GG]A", "T[TA>GG]C", "T[TA>GG]G",
      "T[TA>GG]T", "T[TA>GT]A", "T[TA>GT]C", "T[TA>GT]G", "T[TA>GT]T",
      "T[TC>AA]A", "T[TC>AA]C", "T[TC>AA]G", "T[TC>AA]T", "T[TC>AG]A",
      "T[TC>AG]C", "T[TC>AG]G", "T[TC>AG]T", "T[TC>AT]A", "T[TC>AT]C",
      "T[TC>AT]G", "T[TC>AT]T", "T[TC>CA]A", "T[TC>CA]C", "T[TC>CA]G",
      "T[TC>CA]T", "T[TC>CG]A", "T[TC>CG]C", "T[TC>CG]G", "T[TC>CG]T",
      "T[TC>CT]A", "T[TC>CT]C", "T[TC>CT]G", "T[TC>CT]T", "T[TC>GA]A",
      "T[TC>GA]C", "T[TC>GA]G", "T[TC>GA]T", "T[TC>GG]A", "T[TC>GG]C",
      "T[TC>GG]G", "T[TC>GG]T", "T[TC>GT]A", "T[TC>GT]C", "T[TC>GT]G",
      "T[TC>GT]T", "T[TG>AA]A", "T[TG>AA]C", "T[TG>AA]G", "T[TG>AA]T",
      "T[TG>AC]A", "T[TG>AC]C", "T[TG>AC]G", "T[TG>AC]T", "T[TG>AT]A",
      "T[TG>AT]C", "T[TG>AT]G", "T[TG>AT]T", "T[TG>CA]A", "T[TG>CA]C",
      "T[TG>CA]G", "T[TG>CA]T", "T[TG>CC]A", "T[TG>CC]C", "T[TG>CC]G",
      "T[TG>CC]T", "T[TG>CT]A", "T[TG>CT]C", "T[TG>CT]G", "T[TG>CT]T",
      "T[TG>GA]A", "T[TG>GA]C", "T[TG>GA]G", "T[TG>GA]T", "T[TG>GC]A",
      "T[TG>GC]C", "T[TG>GC]G", "T[TG>GC]T", "T[TG>GT]A", "T[TG>GT]C",
      "T[TG>GT]G", "T[TG>GT]T", "T[TT>AA]A", "T[TT>AA]C", "T[TT>AA]G",
      "T[TT>AA]T", "T[TT>AC]A", "T[TT>AC]C", "T[TT>AC]G", "T[TT>AC]T",
      "T[TT>AG]A", "T[TT>AG]C", "T[TT>AG]G", "T[TT>AG]T", "T[TT>CA]A",
      "T[TT>CA]C", "T[TT>CA]G", "T[TT>CA]T", "T[TT>CC]A", "T[TT>CC]C",
      "T[TT>CC]G", "T[TT>CC]T", "T[TT>CG]A", "T[TT>CG]C", "T[TT>CG]G",
      "T[TT>CG]T", "T[TT>GA]A", "T[TT>GA]C", "T[TT>GA]G", "T[TT>GA]T",
      "T[TT>GC]A", "T[TT>GC]C", "T[TT>GC]G", "T[TT>GC]T", "T[TT>GG]A",
      "T[TT>GG]C", "T[TT>GG]G", "T[TT>GG]T"
    )
  )
}


catalogue_to_wide <- function(df_catalogue, class, col_sample = NULL){
  tidyr::pivot_wider(df_catalogue, names_from = "channel", values_from = fraction, id_cols = col_sample, names_prefix = paste0(class, "|"))
}
