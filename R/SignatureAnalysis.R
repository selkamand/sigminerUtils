# Mutational Signature Analysis -------------------------------------------

#' Mutational Signature Analysis
#'
#' Run all signature mutation analyses possible from MAF inputs on
#'
#' @param maf The input MAF file. Can be a maf object or the path to a MAF file
#' @param copynumber The input copynumber data.frame. See [sigminer::read_copynumber()] for details.
#' @param structuralvariant The input structural variant data.frame. See [sigminer::read_sv_as_rs()] for details
#' @param ref A character vector specifying the reference genome. One of 'hg38' or 'hg19'.
#' @param output_dir The output directory for storing results. Default is "./signatures".
#' @param exposure_type The type of exposure. Can be "absolute" or "relative". One of "absolute" or "relative"
#' @param n_bootstraps The number of bootstrap iterations for fitting signatures. Default is 100.
#' @param temp_dir The temporary directory for storing intermediate files. Default is tempdir().
#' @param db_sbs,db_indel,db_dbs,db_cn a signature collection data.frame where rows are channels and columns are signatures.
#' @param ref_tallies path to a parquet file describing catalogues of a reference database. Can be produced from a folder full of sigminerUtils signature outputs using [sig_create_reference_set()].
#' If building yourself, it must contain columns class,sample,channel,type,fraction,count. If building your own, we recommend partitioning on class then sample.
#' @param cores Number of cores to use.
#' @return None.
#' @export
#' @importFrom rlang `%||%`
#'
sig_analyse_mutations <- function(
    maf, copynumber = NULL, structuralvariant = NULL,
    db_sbs = NULL, db_indel = NULL, db_dbs = NULL, db_cn = NULL, db_sv = NULL,
    ref_tallies = NULL,
    ref = c('hg38', 'hg19'), output_dir = "./signatures", exposure_type = c("absolute", "relative"),
    n_bootstraps = 100, temp_dir = tempdir(), cores = future::availableCores()){

  cli::cli_h1("Mutational Signature Analysis")
  cli::cli_h2("Checking arguments")
  ref <- rlang::arg_match(ref)
  exposure_type <- rlang::arg_match(exposure_type)
  if(!is.null(copynumber)) { assertions::assert_dataframe(copynumber); cn=TRUE} else cn = FALSE
  if(!is.null(structuralvariant)) { assertions::assert_dataframe(structuralvariant); sv = TRUE} else sv = FALSE


  # Define default signature collections based on reference genome
  if(ref == "hg38"){
    default_sbs =  sigstash::sig_load("COSMIC_v3.4_SBS_GRCh38", format = "sigminer")
    default_dbs = sigstash::sig_load("COSMIC_v3.4_DBS_GRCh38", format = "sigminer")
    default_indel = sigstash::sig_load("COSMIC_v3.4_ID_GRCh37", format = "sigminer") # no hg38 renormalised data is available in cosmic
    default_cn = sigstash::sig_load("COSMIC_v3.4_CN_GRCh37", format = "sigminer") #no hg38 renormalised data is available in cosmic
    default_sv = sigstash::sig_load("COSMIC_v3.4_SV_GRCh38", format = "sigminer")
  }
  else if (ref == "hg19"){
    default_sbs = sigstash::sig_load("COSMIC_v3.4_SBS_GRCh37", format = "sigminer")
    default_dbs = sigstash::sig_load("COSMIC_v3.4_DBS_GRCh37", format = "sigminer")
    default_indel = sigstash::sig_load("COSMIC_v3.4_ID_GRCh37", format = "sigminer")
    default_cn = sigstash::sig_load("COSMIC_v3.4_CN_GRCh37", format = "sigminer")
    default_sv = sigstash::sig_load("COSMIC_v3.4_SV_GRCh38", format = "sigminer") #no hg37 renormalised data is available in cosmic
  }
  else
    stop('Unexpected value of ref: ', ref)

  # Set signature collection to defaults if null.
  db_sbs <- db_sbs %||% default_sbs
  db_indel <- db_indel %||% default_indel
  db_dbs <- db_dbs %||% default_dbs
  db_cn <- db_cn %||% default_cn
  db_sv <- db_sv %||% default_sv

  # Pick appropriate reference gene
  if(ref == "hg19"){
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  }else if (ref == "hg38"){
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  }else
    stop("Unknown reference genome", ref)

  # Create Output directory
  if(!file.exists(output_dir)){
    cli::cli_alert_info("Creating Output Directory at: {.path {output_dir}}")
    dir.create(output_dir)
  }
  cli::cli_alert_info("Output Directory: {.path {output_dir}}")
  cli::cli_alert_info("Reference Genome: {.strong {ref_genome}}")

  # Read MAF file if supplied as filepath
  if(is.character(maf))
   maf <- readr::read_tsv(maf, show_col_types = FALSE)

  # If maf is a data.frame, read it into a maf object
  if(is.data.frame(maf)){
    maf <- sigminer::read_maf_minimal(maf)
  }

  # By this point the maf variable must contain a MAF object
  assertions::assert_class(maf, class = "MAF", msg = "maf input in an unexpected format. Please supply maf argument as either a path to a MAF file, a data.frame with MAF columns, or a MAF object from maftools. If you're input format is acceptable please check that your file/object conforms to the MAF specification")

  cli::cli_h2("Tally (Small Variants)")
  tally <- sigminer::sig_tally(
    object = maf,
    mode = "ALL",
    ref_genome = ref_genome,
    keep_only_matrix = FALSE,
    cores = cores
  )

  if(cn){
    cli::cli_h2("Tally (CopyNumber)")

   cn_object <- sigminer::read_copynumber(
     input = copynumber,
     loh_min_len = 10000,
     loh_min_frac = 0.05,
     join_adj_seg = FALSE, # Must be set this way for steele tally method
     genome_measure = "wg",
     genome_build = ref,
     samp_col = "sample",
     complement = TRUE,
     add_loh = TRUE
     )

   tally_copynumber_steele <- sigminer::sig_tally(cn_object, method = "S", ref_genome=ref_genome, cores = cores)
   tally_copynumber_wang <- sigminer::sig_tally(cn_object, method = "W", ref_genome=ref_genome, cores = cores)
   tally_copynumber_tao <- sigminer::sig_tally(cn_object, method = "X", ref_genome=ref_genome, cores = cores)
  }

  if(sv){
    cli::cli_h2("Tally (SV)")
    sv_object <- sigminer::read_sv_as_rs(structuralvariant)
    tally_sv <- sigminer::sig_tally(sv_object, ref_genome=ref_genome, cores = cores)
  }

  cli::cli_h2("Fitting")
  sbs_96_matrices <- t(tally$SBS_96)
  id_83_matrices <- t(tally$ID_83)
  dbs_78_matrices <- t(tally$DBS_78)

  # Matrices without sig databases
  sbs_1536_matrices <- t(tally$SBS_1536)

  if(cn) {
    cn_48_matrices <- t(tally_copynumber_steele$all_matrices$CN_48)

    # Matrices without sig databases (might need to replace with simpler versions
    cn_80_wang_matrices <- t(tally_copynumber_wang$nmf_matrix)
    cn_176_tao_matrices <- t(tally_copynumber_tao$all_matrices$standard_matrix)
  }

  if(sv){
    sv_32_matrices <- t(tally_sv$all_matrices$RS_32)

    # Matrix without a sig database
    sv_38_matrices <- t(tally_sv$all_matrices$RS_38)
  }

  # Sort Signature databases so that rows match out sample catalogues
  db_sbs <- sort_so_rownames_match(db_sbs, rowname_desired_order = rownames(sbs_96_matrices))
  db_indel <- sort_so_rownames_match(db_indel, rowname_desired_order = rownames(id_83_matrices))
  db_dbs <- sort_so_rownames_match(db_dbs, rowname_desired_order = rownames(dbs_78_matrices))
  if(cn) db_cn <- sort_so_rownames_match(db_cn, rowname_desired_order = rownames(cn_48_matrices))
  if(sv) db_sv <- sort_so_rownames_match(db_sv, rowname_desired_order = rownames(sv_32_matrices))

  assertions::assert_identical(rownames(sbs_96_matrices), rownames(db_sbs))
  assertions::assert_identical(rownames(id_83_matrices), rownames(db_indel))
  assertions::assert_identical(rownames(dbs_78_matrices), rownames(db_dbs))
  if(cn) assertions::assert_identical(rownames(cn_48_matrices), rownames(db_cn))
  if(sv) assertions::assert_identical(rownames(sv_32_matrices), rownames(db_sv))


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
    use_parallel = cores,
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

  if(cn){
    cli::cli_h3("Copy Number Alterations (CN48)")
    cn48_fit <- sigminer::sig_fit_bootstrap_batch(
      catalog = cn_48_matrices,
      sig = db_cn,
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
      #mode = "copynumber",
      type = exposure_type # could be 'relative' or 'absolute'
    )
  }

  if(sv){
    #browser()
    cli::cli_h3("Structural Variants (SV32)")
    sv32_fit <- sigminer::sig_fit_bootstrap_batch(
      catalog = sv_32_matrices,
      sig = db_sv,
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
      type = exposure_type # could be 'relative' or 'absolute'
    )
  }

  cli::cli_h2("Write Output")
  cli::cli_h3("Raw counts (tally)")

  # Longify tally Matrices (result of prepare_matrix)
  fix_tally <- function(matrix){
    matrix |>
      as.data.frame() |>
      tibble::rownames_to_column(var="Context") |>
      tidyr::pivot_longer(cols = -1, names_to = "SampleID", values_to = "Count") |>
      dplyr::mutate(CountRelative = Count / sum(Count), .by = SampleID) |>
      dplyr::relocate(SampleID)
  }

  write_tally_matrix <- function(matrix, class, output_dir, ref){
    samples = colnames(tally_ls[[1]])

    df_longform_tally <- fix_tally(matrix)
    ls_longform <- split(df_longform_tally, f = df_longform_tally[["SampleID"]])
    for (sample in names(ls_longform)){
      tmp_tally_outfile = glue::glue("{output_dir}/{class}_catalogue.{sample}.{ref}.tally.csv.gz")
      df_tally <- dplyr::select(ls_longform[[sample]], -SampleID)
      df_tally <- dplyr::rename(df_tally, channel=Context, count=Count, fraction=CountRelative)
      #if(class == "SV38") browser()
      df_tally[["channel"]] <- sigstash::sig_convert_channel_name(df_tally[["channel"]], from = "sigminer", to = "cosmic")
      df_tally[["type"]] <- sigstash::sig_convert_channel2type(df_tally[["channel"]], sigclass = class)
      df_tally[["fraction"]] <- ifelse(is.na(df_tally[["fraction"]]), yes = 0, no = df_tally[["fraction"]])
      df_tally <- df_tally[c("channel", "type", "fraction", "count")]

      write_compressed_csv(
        x = df_tally,
        file = tmp_tally_outfile
      )
      cli::cli_alert_success("{class} tally written to csv: {.path {tmp_tally_outfile}}")

      # Compare to Reference Set
      if(!is.null(ref_tallies)){
        tally_similarity_outfile = glue::glue("{output_dir}/{class}_comparison.{sample}.{ref}.similarity.csv.gz")
        cli::cli_alert_info("Computing tally similarity to reference dataset: {.file {ref_tallies}}")

        #if(sample == "TCGA-CA-6717-01" & class == "ID83") browser()

        df_similarity <- compute_similarity_against_reference_set(tally_file = tmp_tally_outfile, ref_tallies = ref_tallies)
        write_compressed_csv(
          x = df_similarity,
          file = tally_similarity_outfile
        )
        cli::cli_alert_success("Similarities written to csv : {.path {tmp_tally_outfile}.gz}")
      }
    }
  }


  # Create Tally List
  tally_ls <- list(
    "SBS96" = sbs_96_matrices,
    "SBS1536" = sbs_1536_matrices,
    "ID83" = id_83_matrices,
    "DBS78" = dbs_78_matrices
  )
  if(cn) {
    tally_ls[["CN48"]] <- cn_48_matrices
    tally_ls[["CN80"]] <- cn_80_wang_matrices
    tally_ls[["CN176"]] <- cn_176_tao_matrices
  }
  if(sv){
   tally_ls[["SV32"]] <- sv_32_matrices
   tally_ls[["SV38"]] <- sv_38_matrices
  }

  # Write each matrix
  for (class in names(tally_ls)){
    write_tally_matrix(matrix = tally_ls[[class]], class = class, output_dir=output_dir,ref=ref)
  }

  cli::cli_h3("Fit (Exposures)")

  # Functions to pluck different types
  bootstrap_pluck_expo <- function(fit){
    fit$expo |>
      dplyr::rename(SampleID = sample, Contribution = exposure, Sig=sig, Method = method, Type = type)  |>
      dplyr::filter(Type == "optimal") |>
      dplyr::mutate(ContributionRelative = Contribution / sum(Contribution, na.rm = TRUE), .by = c(SampleID, Type)) |>
      #dplyr::mutate(IsOptimal = Type == "optimal") |>
      dplyr::relocate(.after = dplyr::everything(), c(Contribution, ContributionRelative))
  }

  bootstrap_pluck_expo_bootstraps <- function(fit){
    fit$expo |>
      dplyr::rename(SampleID = sample, Contribution = exposure, Sig=sig, Method = method, Type = type)  |>
      dplyr::filter(Type != "optimal") |>
      dplyr::mutate(ContributionRelative = Contribution / sum(Contribution, na.rm = TRUE), .by = c(SampleID, Type)) |>
      #dplyr::mutate(IsOptimal = Type == "optimal") |>
      dplyr::relocate(.after = dplyr::everything(), c(Contribution, ContributionRelative))

  }

  bootstrap_pluck_error_and_cosine <- function(fit){
    df_error <- fit$error |>
      dplyr::rename(Errors = errors, SampleID = sample, Type = type, Method = method) |>
      dplyr::filter(Type == "optimal") |>
      #dplyr::mutate(IsOptimal = Type == "optimal") |>
      dplyr::relocate(Method, SampleID, Type)

    df_cosine <- fit$cosine |>
      dplyr::rename(Cosine = cosine, SampleID = sample, Type = type, Method = method) |>
      dplyr::filter(Type == "optimal") |>
      dplyr::relocate(Method, SampleID, Type)

    df_error_and_cosine <- dplyr::left_join(x = df_error, y = df_cosine, by = c("SampleID", "Method", "Type"))
    return(df_error_and_cosine)
  }

  bootstrap_pluck_error_and_cosine_bootstraps <- function(fit){
    df_error <- fit$error |>
      dplyr::rename(Errors = errors, SampleID = sample, Type = type, Method = method) |>
      dplyr::filter(Type != "optimal") |>
      #dplyr::mutate(IsOptimal = Type == "optimal") |>
      dplyr::relocate(Method, SampleID, Type)

    df_cosine <- fit$cosine |>
      dplyr::rename(Cosine = cosine, SampleID = sample, Type = type, Method = method) |>
      dplyr::filter(Type != "optimal") |>
      dplyr::relocate(Method, SampleID, Type)

    df_error_and_cosine <- dplyr::left_join(x = df_error, y = df_cosine, by = c("SampleID", "Method", "Type"))
    return(df_error_and_cosine)
  }

  bootstrap_pluck_pval <- function(fit){
    fit$p_val |>
      dplyr::rename(Threshold = threshold, Sig = sig, SampleID = sample, Method = method, Pvalue = p_value) |>
      dplyr::relocate(Method, SampleID)
  }

  bootstrap_pluck_summary <- function(fit){
    # Summarise the bootstraps by signature with all the information required to filter sigs.
    df_bootstraps <- bootstrap_pluck_expo_bootstraps(fit)
    df_summary <- df_bootstraps |>
      dplyr::summarise(
        boxplotstats::calculate_boxplot_stats(ContributionRelative, outliers_as_strings = TRUE),
        experimental_pval = sigstats::sig_compute_experimental_p_value(ContributionRelative, threshold = 0.05),
        .by = c(SampleID, Sig)
        )
    df_summary
    return(df_summary)
  }


  write_model_outputs <- function(output_dir, fit, fit_type = "SBS96", ref){

    for (fit_metric in c(
      "expo", "expo_bootstraps", "bootstrap_summary" , "error_and_cosine", "error_and_cosine_bootstraps",  "p_val")
      ){

      if(fit_metric == "expo")
        res <- bootstrap_pluck_expo(fit)
      else if(fit_metric == "expo_bootstraps")
        res <- bootstrap_pluck_expo_bootstraps(fit)
      else if(fit_metric == "bootstrap_summary")
        res <- bootstrap_pluck_summary(fit)
      else if(fit_metric == "error_and_cosine")
        res <- bootstrap_pluck_error_and_cosine(fit)
      else if(fit_metric == "error_and_cosine_bootstraps")
        res <- bootstrap_pluck_error_and_cosine_bootstraps(fit)
      else if(fit_metric == "p_val")
        res <- bootstrap_pluck_pval(fit)

      ls_res <- split(res, f = res[["SampleID"]])
      for (sample in names(ls_res)) {
        tmp_outfile=glue::glue("{output_dir}/{fit_type}_fit.{sample}.{ref}.{fit_metric}.csv")
        write_compressed_csv(
          x=dplyr::select(.data = ls_res[[sample]], -SampleID),
          tmp_outfile
          )
        cli::cli_alert_success("{fit_type} model fit [{fit_metric}] has been written to csv: {.path {tmp_outfile}.gz}")
      }
    }
  }

  write_model_outputs(fit = sbs96_fit, fit_type = "SBS96", output_dir = output_dir, ref = ref)
  write_model_outputs(fit = dbs78_fit, fit_type = "DBS78", output_dir = output_dir, ref = ref)
  write_model_outputs(fit = id83_fit, fit_type = "ID83", output_dir = output_dir, ref = ref)
  if(cn) write_model_outputs(fit = cn48_fit, fit_type = "CN48", output_dir = output_dir, ref = ref)
  if(sv) write_model_outputs(fit = cn48_fit, fit_type = "SV32", output_dir = output_dir, ref = ref)

}


sort_so_rownames_match <- function(data, rowname_desired_order){
  assertions::assert(nrow(data) == length(rowname_desired_order), msg = "Number of sample catalog rows must be equal to rows of signature collection dataframe")
  rnames = rownames(data)
  indexes = match(rowname_desired_order, rnames)
  assertions::assert_no_missing(indexes, msg = "Sample Catalog Matrix contains rows not present in the signature collection dataframe")

  return(data[indexes,])
}
