# Mutational Signature Analysis -------------------------------------------

#' Mutational Signature Analysis
#'
#' Screen datasets for known mutational signatures by
#' supplying both sample mutation data and collections of known
#' signatures.
#'
#' @param maf The input MAF file. Can be a maf object or the path to a MAF file
#' @param copynumber The input copynumber data.frame. See [sigminer::read_copynumber()] for details.
#' @param structuralvariant The input structural variant data.frame. See [sigminer::read_sv_as_rs()] for details
#' @param ref A character vector specifying the reference genome. One of 'hg38' or 'hg19'.
#' @param output_dir The output directory for storing results. Default is "./signatures".
#' @param exposure_type The type of exposure. Can be "absolute" or "relative". One of "absolute" or "relative"
#' @param n_bootstraps The number of bootstrap iterations for fitting signatures. Default is 100.
#' @param temp_dir The temporary directory for storing intermediate files. Default is tempdir().
#' @param min_contribution_threshold The minimum contribution threshold a signature must surpass for it to count as 'present' in a bootstrap (typically 0.05). This is sometimes referred to as the sparsity threshold. See [sigstats::sig_compute_experimental_p_value()] for details.
#' @param db_sbs,db_indel,db_dbs,db_cn,db_sv a signature collection data.frame where rows are channels and columns are signatures. Row names should be signature channels.
#' See \code{sigstash::sig_load("COSMIC_v3.4_SBS_GRCh38", format = "sigminer")}  for an example.
#' Alternatively, you can supply a path to signature collections in tidy_csv format (see [sigstash::sig_read_signatures()] for details)
#' @param db_sbs_name,db_indel_name,db_dbs_name,db_cn_name,db_sv_name names of signature databases (used for logs). If NULL, will try and lookup using [sigstash::sig_identify_collection()].
#' @param ref_tallies path to a parquet file describing catalogues of a reference database. Can be produced from a folder full of sigminerUtils signature outputs using [sig_create_reference_set()].
#' If building yourself, it must contain columns class,sample,channel,type,fraction,count. If building your own, we recommend partitioning on class then sample.
#' @param cores Number of cores to use.
#' @return None.
#' @export
#' @importFrom rlang `%||%`
#'
#' @details
#' Produces the following files:
#'
#' Per Run:
#' signature_collections.csv: Describes the signature collection used COSMIC_v3.4_SBS_GRCh38
#'  thresholds.csv: Describes any the thresholds that were used
#'
#' For each signature class:
#' <signatureclass>_catalogue.<sampleid>.<refgenome>.tally.csv.gz: A tally of mutations in the sample by channel.
#' <signatureclass>_fit.<sampleid>.<refgenome>.expo.csv.gz: The signature model produced by sigminer (% contribution of each signature). Annotated with p-val from bootstrap_summary.csv.gz
#' <signatureclass>_fit.<sampleid>.<refgenome>.bootstrap_summary.csv.gz: A signature-level summary of stability across bootstraps. Includes experimentally defined P-values.
#' <signatureclass>_fit.<sampleid>.<refgenome>.expo_bootstraps.csv.gz: The optimal signature model produced in each bootstrap.
#' <signatureclass>_fit.<sampleid>.<refgenome>.error_and_cosine.csv.gz: Error and cosine similarity of signature model VS original data (unfiltered model, so may include unstable signatures)
#' <signatureclass>_fit.<sampleid>.<refgenome>.error_and_cosine_bootstraps.csv.gz: Error and cosine similarity of each bootstraps optimal signature model vs original data (unfiltered model, so may include unstable signatures).
#' <signatureclass>_fit.<sampleid>.<refgenome>.p_val.csv.gz: Per-signature p value computed by sigminer. Note this is different to the p-value in bootstrap_summary and expo datasets (we recommend using the latter, which has been added for greater consistency with other signature pipelines).
#'
sig_analyse_mutations <- function(
    maf, copynumber = NULL, structuralvariant = NULL,
    db_sbs = NULL, db_indel = NULL, db_dbs = NULL, db_cn = NULL, db_sv = NULL,
    db_sbs_name = NULL, db_indel_name = NULL, db_dbs_name = NULL, db_cn_name = NULL, db_sv_name = NULL,
    ref_tallies = NULL,
    min_contribution_threshold = 0.05,
    ref = c('hg38', 'hg19'), output_dir = "./signatures", exposure_type = c("absolute", "relative"),
    n_bootstraps = 100, temp_dir = tempdir(), cores = future::availableCores()){

  # TODO: REMOVE locale setting once sigstash issue https://github.com/selkamand/sigstash/issues/43 is resolved
  Sys.setlocale("LC_COLLATE", "C")

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

  # If user supplies string as a signature db, assume its a file (in csv_tidy format) and parse it to a sigminer-compatible format
  db_sbs <- db_read_if_filepath(db_sbs, dbtype="db_sbs")
  db_indel <- db_read_if_filepath(db_indel, dbtype="db_indel")
  db_dbs <- db_read_if_filepath(db_dbs, dbtype="db_dbs")
  db_cn <- db_read_if_filepath(db_cn, dbtype="db_cn")
  db_sv <- db_read_if_filepath(db_sv, dbtype="db_sv")

  # If signature collection names are NULL, see if we can identify based on their md5
  if(!is.null(db_sbs)) db_sbs_name <- db_sbs_name %||% sigstash::sig_identify_collection(db_sbs, return = 'name')
  if(!is.null(db_indel)) db_indel_name <- db_indel_name %||% sigstash::sig_identify_collection(db_indel, return = 'name')
  if(!is.null(db_dbs)) db_dbs_name <- db_dbs_name %||% sigstash::sig_identify_collection(db_dbs, return = 'name')
  if(!is.null(db_cn)) db_cn_name <- db_cn_name %||% sigstash::sig_identify_collection(db_cn, return = 'name')
  if(!is.null(db_sv)) db_sv_name <- db_sv_name %||% sigstash::sig_identify_collection(db_sv, return = 'name')


  # If signature collection names couldn't be identified, throw an error
  assertions::assert(is.null(db_sbs) | db_sbs_name != "Uncertain", msg = "For logging purposes, signature collections must be named. Either supply name to {.arg db_sbs_name} argument or use a sigstash collection whose name can be identified using sigstash::sig_identify_collection()")
  assertions::assert(is.null(db_indel) | db_indel_name != "Uncertain", msg = "For logging purposes, signature collections must be named. Either supply name to {.arg db_indel_name} argument or use a sigstash collection whose name can be identified using sigstash::sig_identify_collection()")
  assertions::assert(is.null(db_dbs) | db_dbs_name != "Uncertain", msg = "For logging purposes, signature collections must be named. Either supply name to {.arg db_dbs_name} argument or use a sigstash collection whose name can be identified using sigstash::sig_identify_collection()")
  assertions::assert(is.null(db_cn) | db_cn_name != "Uncertain", msg = "For logging purposes, signature collections must be named. Either supply name to {.arg db_cn_name} argument or use a sigstash collection whose name can be identified using sigstash::sig_identify_collection()")
  assertions::assert(is.null(db_sv) | db_sv_name != "Uncertain", msg = "For logging purposes, signature collections must be named. Either supply name to {.arg db_sv_name} argument or use a sigstash collection whose name can be identified using sigstash::sig_identify_collection()")


  # Pick appropriate reference genome
  if(ref == "hg19"){
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  }else if (ref == "hg38"){
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  }else
    stop("Unknown reference genome", ref)

  # Create Output directory
  if(!file.exists(output_dir)){
    cli::cli_alert_info("Creating Output Directory at: {.path {output_dir}}")
    dir.create(output_dir, recursive = TRUE)
  }
  cli::cli_alert_info("Output Directory: {.path {output_dir}}")
  cli::cli_alert_info("Reference Genome: {.strong {ref_genome}}")

  # Log which signature databases were used
  sigdbs_log <- glue::glue_safe("{output_dir}/signature_collections.csv")
  file.create(sigdbs_log)
  dbname_map = c(
    "SBS96" = db_sbs_name,
    "ID83" = db_indel_name,
    "DBS78" = db_dbs_name,
    "CN48" = db_cn_name,
    "SV32" = db_sv_name
  )
  db_names <- data.frame(collection_type = names(dbname_map), dataset = unname(dbname_map))
  write.csv(db_names, file = sigdbs_log, row.names = FALSE)

  # Log which thresholds were used
  thresholds_log <- glue::glue_safe("{output_dir}/thresholds.csv")
  file.create(thresholds_log)
  thresholds_map = c(
    "min_contribution_threshold" = min_contribution_threshold
  )
  df_thresholds <- data.frame(threshold = names(thresholds_map), value = unname(thresholds_map))
  write.csv(df_thresholds, file = thresholds_log, row.names = FALSE)


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
  tally_sbs <- tally(
    object = maf,
    mode = "SBS",
    ref_genome = ref_genome,
    cores = cores
  )

  tally_dbs <- tally(
    object = maf,
    mode = "DBS",
    ref_genome = ref_genome,
    cores = cores
  )

  tally_id <- tally(
    object = maf,
    mode = "ID",
    ref_genome = ref_genome,
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

  samples <- unique(maf@data$Tumor_Sample_Barcode)
  sbs_96_matrices <- postprocess_tally_matrix(tally_sbs$all_matrices$SBS_96, samples = samples, class = "SBS96")
  id_83_matrices <- postprocess_tally_matrix(tally_id$all_matrices$ID_83, samples = samples, class = "ID83")
  dbs_78_matrices <- postprocess_tally_matrix(tally_dbs$all_matrices$DBS_78, samples = samples, class = "DBS78")

  # Matrices without sig databases
  sbs_1536_matrices <-  postprocess_tally_matrix(tally_sbs$all_matrices$SBS_1536, samples = samples, class = "SBS1536")
  dbs_1248_matrices <- postprocess_tally_matrix(tally_dbs$all_matrices$DBS_1248, samples = samples, class = "DBS1248")


  if(cn) {
    cn_48_matrices <- t(tally_copynumber_steele$all_matrices$CN_48)

    # Matrices without sig databases (might need to replace with simpler versions
    cn_176_tao_matrices <- t(tally_copynumber_tao$all_matrices$standard_matrix)
    cn_80_wang_matrices <- t(tally_copynumber_wang$nmf_matrix)

    #Hotfix for cn80_wang_matrices which lack a sample names.
    # Can Remove once issue https://github.com/ShixiangWang/sigminer/issues/465 is resolved.
    # Can test by checking if there CN80 tally includes sample 'V1')
    colnames(cn_80_wang_matrices) <- colnames(cn_48_matrices)
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
    p_val_thresholds = c(0.05), # This value does inform how sigminer 'p-values' in the `p_val.csv` file are computed,
    # we actually don't use these p-values for any downstream filtering/calculation.
    # Rather we use sigstats::sig_compute_experimental_p_value() based on min_contribution_threshold
    # Which can be seen in the bootstrap_summary.csv
    use_parallel = cores,
    seed = 123456L,
    job_id = NULL,
    result_dir = temp_dir,
    #mode = "SBS",
    type = exposure_type # could be 'relative' or 'absolute'
  ) |> try() |> try_error_to_null()

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
  ) |> try() |> try_error_to_null()

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
  ) |> try() |> try_error_to_null()

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
    ) |> try() |> try_error_to_null()
  }

  if(sv){
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
    ) |> try() |> try_error_to_null()
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

    # Pluck the summary so we can annotate signatures with their bootstrap-computed p-value
    df_bootstrap_summary <- dplyr::select(.data = bootstrap_pluck_summary(fit), SampleID, Sig, experimental_pval)

    fit$expo |>
      dplyr::rename(SampleID = sample, Contribution = exposure, Sig=sig, Method = method, Type = type)  |>
      dplyr::filter(Type == "optimal") |>
      dplyr::select(-Type, -Method) |>
      dplyr::mutate(ContributionRelative = Contribution / sum(Contribution, na.rm = TRUE), .by = c(SampleID, Type)) |>
      dplyr::left_join(df_bootstrap_summary, by = c("SampleID", "Sig")) |>
      dplyr::relocate(.after = dplyr::everything(), c(Contribution, ContributionRelative))
  }

  bootstrap_pluck_expo_bootstraps <- function(fit){
    fit$expo |>
      dplyr::rename(SampleID = sample, Contribution = exposure, Sig=sig, Method = method, Type = type)  |>
      dplyr::filter(Type != "optimal") |>
      dplyr::mutate(ContributionRelative = Contribution / sum(Contribution, na.rm = TRUE), .by = c(SampleID, Type)) |>
      dplyr::relocate(.after = dplyr::everything(), c(Contribution, ContributionRelative))

  }

  bootstrap_pluck_error_and_cosine <- function(fit){
    df_error <- fit$error |>
      dplyr::rename(Errors = errors, SampleID = sample, Type = type, Method = method) |>
      dplyr::filter(Type == "optimal") |>
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
        experimental_pval = sigstats::sig_compute_experimental_p_value(ContributionRelative, threshold = min_contribution_threshold),
        .by = c(SampleID, Sig)
        )
    df_summary
    return(df_summary)
  }


  write_model_outputs <- function(output_dir, fit, fit_type = "SBS96", ref){

    # Skip if fit is NULL (means there were no mutations of fit_type)
    if(is.null(fit)) {
      cli::cli_alert_warning("No model outputs produced for class {fit_type}")
      return(NULL)
    }

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
  if(sv) write_model_outputs(fit = sv32_fit, fit_type = "SV32", output_dir = output_dir, ref = ref)
}


sort_so_rownames_match <- function(data, rowname_desired_order){
  assertions::assert(nrow(data) == length(rowname_desired_order), msg = "Number of sample catalog rows must be equal to rows of signature collection dataframe")
  rnames = rownames(data)
  indexes = match(rowname_desired_order, rnames)
  assertions::assert_no_missing(indexes, msg = "Sample Catalog Matrix contains rows not present in the signature collection dataframe")

  return(data[indexes,])
}

#' Mutational Signature Analysis
#'
#' Run all signature mutation analyses possible from file inputs.
#'
#' @inheritParams sig_analyse_mutations
#' @inheritParams sigstart::parse_purple_cnv_to_sigminer
#' @inheritParams sigstart::parse_purple_sv_vcf_to_sigminer
#' @inheritParams sigstart::parse_vcf_to_sigminer_maf
#' @param sample_id string representing the tumour sample identifier (in your VCFs and other files).
#' @return None.
#' @export
#'
#' @examples
#' \dontrun{
#' path_snvs <- system.file(
#'   "colo829_testfiles/COLO829v003T.purple.somatic.vcf.gz",
#'   package = "sigminerUtils"
#' )
#' path_cnvs <- system.file(
#'   "colo829_testfiles/COLO829v003T.purple.cnv.somatic.tsv",
#'   package = "sigminerUtils"
#' )
#' path_svs <- system.file(
#'   "colo829_testfiles/COLO829v003T.purple.sv.vcf.gz",
#'   package = "sigminerUtils"
#' )
#'
#' sig_analyse_mutations_single_sample_from_files(
#'   sample_id = "COLO829v003T",
#'   vcf_snv = path_snvs,
#'   segment = path_cnvs,
#'   vcf_sv = path_svs,
#'   pass_only = TRUE,
#'   ref = "hg38",
#'   output_dir = "colo829_signature_results"
#' )
#' }
sig_analyse_mutations_single_sample_from_files <- function(
    sample_id,
    vcf_snv = NULL,
    segment = NULL,
    vcf_sv = NULL,
    pass_only = TRUE,
    exclude_sex_chromosomes = TRUE,
    allow_multisample = TRUE,
    db_sbs = NULL, db_indel = NULL, db_dbs = NULL, db_cn = NULL, db_sv = NULL,
    ref_tallies = NULL,
    ref = c('hg38', 'hg19'),
    output_dir = "./signatures",
    exposure_type = c("absolute", "relative"),
    n_bootstraps = 100,
    temp_dir = tempdir(),
    cores = future::availableCores())
  {

    # Check files exist
    if(!is.null(vcf_snv)) assertions::assert_file_exists(vcf_snv)
    if(!is.null(segment)) assertions::assert_file_exists(segment)
    if(!is.null(vcf_sv)) assertions::assert_file_exists(vcf_sv)


    # Create Output Directory
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)
    }

    # Create a directory for this sample
    sample_dir <- paste0(output_dir, "/", sample_id)
    if (file.exists(sample_dir)) {
      error <- glue::glue_safe("Sample folder already exists. To overwrite please manually delete [{sample_dir}]")
      cli::cli_abort(error)
    }
    cli::cli_progress_step("Creating Output Directory at {.file {sample_dir}}")
    dir.create(sample_dir, recursive = TRUE, showWarnings = TRUE)


    # Parse the variant files into sigminer-compatible formats
    small_variants <- if(!is.null(vcf_snv)) sigstart::parse_vcf_to_sigminer_maf(vcf_snv = vcf_snv, sample_id = sample_id, pass_only = pass_only, allow_multisample = allow_multisample) else NULL
    cnvs <- if(!is.null(segment)) sigstart::parse_purple_cnv_to_sigminer(segment = segment, sample_id = sample_id, exclude_sex_chromosomes = exclude_sex_chromosomes) else NULL
    svs <- if(!is.null(vcf_sv)) sigstart::parse_purple_sv_vcf_to_sigminer(vcf_sv = vcf_sv, sample_id = sample_id, pass_only = pass_only) else NULL


    # Create a log with just the most important info about each run
    sample_log <- glue::glue_safe("{sample_dir}/{sample_id}.summary.log")
    file.create(sample_log)

    # Mine versions from VCFs
    snv_versions <- if(!is.null(vcf_snv)) vcf2versions(vcf_snv) else "No SNV VCF supplied"
    sv_versions <- if(!is.null(vcf_sv))  vcf2versions(vcf_sv) else "No SV VCF supplied"

    # Write VCF versions to logfile
    write(glue::glue("\n\nSNV VCF versions:\n{snv_versions}"), sample_log, append = TRUE)
    write(glue::glue("\n\nSV VCF versions:\n{sv_versions}"), sample_log, append = TRUE)


    # Create an additional logfile for sigminerUtils output
    sigminer_log <- glue::glue_safe("{sample_dir}/{sample_id}.sigminer.log")
    file.create(sigminer_log)

    # Sink all sigminer messages to sigminer.log
    # zz <- file(sigminer_log, open = "wt")
    # sink(type = "message", file = zz, split = TRUE)


    cli::cli_h2("Running Signature Analysis. This will take some time")

    try_output <- try({ # Try so we can explicitly log failure
      capture_messages(logfile = sigminer_log, expr = {
          sig_analyse_mutations(
          maf = small_variants,
          copynumber = cnvs,
          structuralvariant = svs,
          db_sbs = db_sbs, db_indel = db_indel, db_dbs = db_dbs, db_cn = db_cn, db_sv = db_sv,
          ref_tallies=ref_tallies,
          ref = ref,
          output_dir = sample_dir,
          exposure_type = exposure_type,
          n_bootstraps = n_bootstraps,
          temp_dir = temp_dir,
          cores = cores
      )})
    })

    # Figure out if sigminer run failed
    sample_failed <- "try-error" %in% class(try_output)

    # Finished Successfully / Failed (update sample_log with result)
    if(sample_failed){
      cli::cli_alert_warning("Sample {sample_id} Failed")
      write(glue::glue_safe("\n\nSample Failed. See {sigminer_log} for details"), sample_log, append = TRUE)
      stop("Signature analysis failed for: ", sample_id, ". Please see ", sigminer_log, " for details")
    }
    else{
      cli::cli_progress_step("Finished successfully")
      write("\n\nFinished successfully!", sample_log, append = TRUE)
    }
}
