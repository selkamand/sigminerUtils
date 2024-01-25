# Utilities --------------------------------------------------------

prepare_matrix <- function(decomp, sample_id = NULL){
  if(is.vector(decomp)){
    res <- as.matrix(decomp)
    colnames(res) <- sample_id
    return(res)
  }
  else
    return(t(decomp))

}

# Mutational Signature Analysis -------------------------------------------

#' Mutational Signature Analysis
#'
#' Run all signature mutation analyses possible from MAF inputs on
#'
#' @param maf The input MAF file.
#' @param somatic_ids A vector of sample IDs for somatic mutations.
#' @param ref A character vector specifying the reference genome. One of 'hg38' or 'hg19'.
#' @param output_dir The output directory for storing results. Default is "./signatures".
#' @param exposure_type The type of exposure. Can be "absolute" or "relative". One of "absolute" or "relative"
#' @param n_bootstraps The number of bootstrap iterations for fitting signatures. Default is 100.
#' @param temp_dir The temporary directory for storing intermediate files. Default is tempdir().
#' @param db_sbs a signature collection data.frame where rows are channels and columns are signatures.
#' @param db_indel a signature collection data.frame where rows are channels and columns are signatures.
#' @param db_dbs a signature collection data.frame where rows are channels and columns are signatures.
#'
#' @return None.
#' @export
#' @importFrom rlang `%||%`
#'
#' @details
#'
#'
sig_analyse_mutations <- function(maf, somatic_ids, db_sbs = NULL, db_indel = NULL, db_dbs = NULL, ref = c('hg38', 'hg19'), output_dir = "./signatures", exposure_type = c("absolute", "relative"), n_bootstraps = 100, temp_dir = tempdir()){

  cli::cli_h1("Mutational Signature Analysis")
  cli::cli_h2("Checking arguments")


  ref <- rlang::arg_match(ref)
  exposure_type <- rlang::arg_match(exposure_type)


  # Define default signature collections based on reference genome
  if(ref == "hg38"){
    default_sbs =  sigstash::sig_load("COSMIC_v3.3.1_SBS_GRCh38", format = "sigminer")
    default_indel = sigstash::sig_load("COSMIC_v3.3_ID_GRCh37", format = "sigminer") # no hg38 renormalised data is available in cosmic
    default_dbs = sigstash::sig_load("COSMIC_v3.3_DBS_GRCh37", format = "sigminer") #no hg38 renormalised data is available in cosmic
  }
  else if (ref == "hg19"){
    default_sbs = sigstash::sig_load("COSMIC_v3.3.1_SBS_GRCh37", format = "sigminer")
    default_indel = sigstash::sig_load("COSMIC_v3.3_ID_GRCh37", format = "sigminer")
    default_dbs = sigstash::sig_load("COSMIC_v3.3_DBS_GRCh37", format = "sigminer")
  }
  else
    stop('Unexpected value of ref: ', ref)

  # Set signaure collection to defaults if null.
  db_sbs <- db_sbs %||% default_sbs
  db_indel <- db_indel %||% default_indel
  db_dbs <- db_dbs %||% default_dbs

  # Pick appropriate reference gene
  if(ref == "hg19"){
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  }else if (ref == "hg38"){
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  }else
    stop("Unknown reference genome", ref)


  if(!file.exists(output_dir)){
    cli::cli_alert_info("Creating Output Directory at: {.path {output_dir}}")
    dir.create(output_dir)
  }
  cli::cli_alert_info("Output Directory: {.path {output_dir}}")
  cli::cli_alert_info("Reference Genome: {.strong {ref_genome}}")

  cli::cli_h2("Decomposition")
  decompositions <- sigminer::sig_tally(
    object = maf,
    mode = "ALL",
    ref_genome = ref_genome,
    keep_only_matrix = FALSE
  )

  cli::cli_h2("Fitting")
  sbs_96_matrices <- prepare_matrix(decompositions$SBS_96, sample_id = somatic_ids)
  id_83_matrices <- prepare_matrix(decompositions$ID_83, sample_id = somatic_ids)
  dbs_78_matrices <- prepare_matrix(decompositions$DBS_78, sample_id = somatic_ids)

  # Sort Signature databases so that rows match out sample catalogues
  db_sbs <- sort_so_rownames_match(db_sbs, rowname_desired_order = rownames(sbs_96_matrices))
  db_indel <- sort_so_rownames_match(db_indel, rowname_desired_order = rownames(id_83_matrices))
  db_dbs <- sort_so_rownames_match(db_dbs, rowname_desired_order = rownames(dbs_78_matrices))

  assertions::assert_identical(rownames(sbs_96_matrices), rownames(db_sbs))
  assertions::assert_identical(rownames(id_83_matrices), rownames(db_indel))
  assertions::assert_identical(rownames(dbs_78_matrices), rownames(db_dbs))


  cli::cli_h3("Single Base Substitutions (SBS96)")
  # Single Base Substitution
  sbs96_fit <- sigminer::sig_fit_bootstrap_batch(
    catalog = sbs_96_matrices,
    sig = db_sbs,
    #sig_db = "SBS",  # Use 'legacy' for V2 or 'SBS' for V3
    #sig_index= "ALL",
    n = n_bootstraps,
    method = "QP",
    min_count = 1L,
    p_val_thresholds = c(0.05),
    use_parallel = TRUE,
    seed = 123456L,
    job_id = NULL,
    result_dir = temp_dir,
    #mode = "SBS",
    type = exposure_type # could be 'relative' or 'absolute'
  )

  cli::cli_h3("INDELS (ID83)")
  ## Indel
  id83_fit <- sigminer::sig_fit_bootstrap_batch(
    catalog = id_83_matrices,
    sig = db_indel,
    #sig_db = "ID",  # Use 'legacy' for V2
    sig_index= NULL,
    n = n_bootstraps,
    method = "QP",
    min_count = 1L,
    p_val_thresholds = c(0.05),
    use_parallel = TRUE,
    seed = 123456L,
    job_id = NULL,
    result_dir = temp_dir,
    mode = "ID",
    type = exposure_type # could be 'relative' or 'absolute'
  )

  cli::cli_h3("Doublet Mutations (DBS78)")
  dbs78_fit <- sigminer::sig_fit_bootstrap_batch(
    catalog = dbs_78_matrices,
    sig = db_dbs,
    #sig_db = "DBS",  # Use 'legacy' for V2
    #sig_index= "ALL",
    n = n_bootstraps,
    method = "QP",
    min_count = 1L,
    p_val_thresholds = c(0.05),
    use_parallel = TRUE,
    seed = 123456L,
    job_id = NULL,
    result_dir = temp_dir,
    mode = "DBS",
    type = exposure_type # could be 'relative' or 'absolute'
  )



  cli::cli_h2("Write Output")
  cli::cli_h3("Raw counts (decompositions)")
  outfile_sbs96_decompositions <- glue::glue("{output_dir}/SBS96_catalogue.{ref}.decomposition.csv")
  outfile_id83_decompositions <- glue::glue("{output_dir}/ID83_catalogue.{ref}.decomposition.csv")
  outfile_dbs78_decompositions <- glue::glue("{output_dir}/DBS78_catalogue.{ref}.decomposition.csv")


  # Longify Decomposition Matrices (result of prepare_matrix)
  fix_decomposition <- function(matrix){
    matrix |>
      as.data.frame() |>
      tibble::rownames_to_column(var="Context") |>
      tidyr::pivot_longer(cols = -1, names_to = "SampleID", values_to = "Count") |>
      dplyr::mutate(CountRelative = Count / sum(Count), .by = SampleID) |>
      dplyr::relocate(SampleID)
  }

  # Add class column
  fix_decomposition(sbs_96_matrices) |>
    utils::write.csv(file = outfile_sbs96_decompositions, row.names = FALSE)
  cli::cli_alert_success("SBS96 decomposition written to csv: {.path {outfile_sbs96_decompositions}}")

  fix_decomposition(id_83_matrices) |>
    utils::write.csv(file = outfile_id83_decompositions, row.names = FALSE)
  cli::cli_alert_success("ID83 decomposition written to csv: {.path {outfile_id83_decompositions}}")

  fix_decomposition(dbs_78_matrices) |>
    utils::write.csv(file = outfile_dbs78_decompositions, row.names = FALSE)
  cli::cli_alert_success("DBS_78 decomposition written to csv: {.path {outfile_dbs78_decompositions}}")


  cli::cli_h3("Fit (Exposures)")

  # Functions
  bootstrap_pluck_expo <- function(fit){
    fit$expo |>
      dplyr::rename(SampleID = sample, Contribution = exposure, Sig=sig, Method = method, Type = type)  |>
      dplyr::mutate(ContributionRelative = Contribution / sum(Contribution, na.rm = TRUE), .by = c(SampleID, Type)) |>
      dplyr::mutate(IsOptimal = Type == "optimal") |>
      dplyr::relocate(.after = dplyr::everything(), c(Contribution, ContributionRelative))
  }

  bootstrap_pluck_error <- function(fit){
    fit$error |>
      dplyr::rename(Errors = errors, SampleID = sample, Type = type, Method = method) |>
      dplyr::mutate(IsOptimal = Type == "optimal") |>
      dplyr::relocate(Method, SampleID, Type, IsOptimal)
  }

  bootstrap_pluck_cosine<- function(fit){
    fit$cosine |>
      dplyr::rename(Cosine = cosine, SampleID = sample, Type = type, Method = method) |>
      dplyr::mutate(IsOptimal = Type == "optimal") |>
      dplyr::relocate(Method, SampleID, Type, IsOptimal)
  }

  bootstrap_pluck_pval <- function(fit){
    fit$p_val |>
      dplyr::rename(Threshold = threshold, Sig = sig, SampleID = sample, Method = method, Pvalue = p_value) |>
      dplyr::relocate(Method, SampleID)
  }


  write_model_outputs <- function(output_dir, fit, fit_type = "SBS96", ref){

    for (fit_metric in c("expo", "error", "cosine", "p_val")){
      tmp_outfile=glue::glue("{output_dir}/{fit_type}_fit.{ref}.{fit_metric}.csv")

      if(fit_metric == "expo")
        res <- bootstrap_pluck_expo(fit)
      else if(fit_metric == "error")
        res <- bootstrap_pluck_error(fit)
      else if(fit_metric == "cosine")
        res <- bootstrap_pluck_cosine(fit)
      else if(fit_metric == "p_val")
        res <- bootstrap_pluck_pval(fit)

      utils::write.csv(res, tmp_outfile, row.names = FALSE)
      cli::cli_alert_success("{fit_type} model fit [{fit_metric}] has been written to csv: {.path {tmp_outfile}}")
    }
  }

  write_model_outputs(fit = sbs96_fit, fit_type = "SBS96", output_dir = output_dir, ref = ref)
  write_model_outputs(fit = dbs78_fit, fit_type = "DBS78", output_dir = output_dir, ref = ref)
  write_model_outputs(fit = id83_fit, fit_type = "ID83", output_dir = output_dir, ref = ref)
}


sort_so_rownames_match <- function(data, rowname_desired_order){
  assertions::assert(nrow(data) == length(rowname_desired_order), msg = "Number of sample catalog rows must be equal to rows of signature collection dataframe")
  rnames = rownames(data)
  indexes = match(rowname_desired_order, rnames)
  assertions::assert_no_missing(indexes, msg = "Sample Catalog Matrix contains rows not present in the signature collection dataframe")

  return(data[indexes,])
}
