# Functions to pluck different types
bootstrap_pluck_expo <- function(fit, min_contribution_threshold){

  # Pluck the summary so we can annotate signatures with their bootstrap-computed p-value
  df_bootstrap_summary <- dplyr::select(.data = bootstrap_pluck_summary(fit, min_contribution_threshold = min_contribution_threshold), sample, signature, experimental_pval)

  fit$expo |>
    dplyr::rename(sample = sample, contribution_absolute = exposure, signature=sig, Method = method, Type = type)  |>
    dplyr::filter(Type == "optimal") |>
    dplyr::mutate(contribution = contribution_absolute / sum(contribution_absolute, na.rm = TRUE), .by = c(sample, Type)) |>
    dplyr::left_join(df_bootstrap_summary, by = c("sample", "signature")) |>
    dplyr::select(-Type, -Method) |>
    dplyr::relocate(.after = dplyr::everything(), c(contribution_absolute, contribution))
}

bootstrap_pluck_expo_bootstraps <- function(fit){
  fit$expo |>
    dplyr::rename(sample = sample, contribution_absolute = exposure, signature=sig, Method = method, Type = type)  |>
    dplyr::filter(Type != "optimal") |>
    dplyr::mutate(contribution = contribution_absolute / sum(contribution_absolute, na.rm = TRUE), .by = c(sample, Type)) |>
    dplyr::relocate(.after = dplyr::everything(), c(contribution_absolute, contribution))

}

bootstrap_pluck_error_and_cosine <- function(fit){
  df_error <- fit$error |>
    dplyr::rename(Errors = errors, sample = sample, Type = type, Method = method) |>
    dplyr::filter(Type == "optimal") |>
    dplyr::relocate(Method, sample, Type)

  df_cosine <- fit$cosine |>
    dplyr::rename(Cosine = cosine, sample = sample, Type = type, Method = method) |>
    dplyr::filter(Type == "optimal") |>
    dplyr::relocate(Method, sample, Type)

  df_error_and_cosine <- dplyr::left_join(x = df_error, y = df_cosine, by = c("sample", "Method", "Type"))
  return(df_error_and_cosine)
}

bootstrap_pluck_error_and_cosine_bootstraps <- function(fit){
  df_error <- fit$error |>
    dplyr::rename(Errors = errors, sample = sample, Type = type, Method = method) |>
    dplyr::filter(Type != "optimal") |>
    #dplyr::mutate(IsOptimal = Type == "optimal") |>
    dplyr::relocate(Method, sample, Type)

  df_cosine <- fit$cosine |>
    dplyr::rename(Cosine = cosine, sample = sample, Type = type, Method = method) |>
    dplyr::filter(Type != "optimal") |>
    dplyr::relocate(Method, sample, Type)

  df_error_and_cosine <- dplyr::left_join(x = df_error, y = df_cosine, by = c("sample", "Method", "Type"))
  return(df_error_and_cosine)
}

bootstrap_pluck_pval <- function(fit){
  fit$p_val |>
    dplyr::rename(Threshold = threshold, signature = sig, sample = sample, Method = method, Pvalue = p_value) |>
    dplyr::relocate(Method, sample)
}

bootstrap_pluck_summary <- function(fit, min_contribution_threshold){
  # Summarise the bootstraps by signature with all the information required to filter sigs.
  df_bootstraps <- bootstrap_pluck_expo_bootstraps(fit)
  df_summary <- df_bootstraps |>
    dplyr::summarise(
      boxplotstats::calculate_boxplot_stats(contribution, outliers_as_strings = TRUE),
      experimental_pval = sigstats::sig_compute_experimental_p_value(contribution, threshold = min_contribution_threshold),
      .by = c(sample, signature)
    )
  return(df_summary)
}


write_model_outputs <- function(output_dir, fit, fit_type = "SBS96", ref, min_contribution_threshold){

  # Skip if fit is NULL (means there were no mutations of fit_type)
  if(is.null(fit)) {
    cli::cli_alert_warning("No model outputs produced for class {fit_type}")
    return(NULL)
  }

  for (fit_metric in c(
    "expo", "expo_bootstraps", "bootstrap_summary" , "error_and_cosine", "error_and_cosine_bootstraps",  "p_val")
  ){
    if(fit_metric == "expo")
      res <- bootstrap_pluck_expo(fit, min_contribution_threshold = min_contribution_threshold)
    else if(fit_metric == "expo_bootstraps")
      res <- bootstrap_pluck_expo_bootstraps(fit)
    else if(fit_metric == "bootstrap_summary")
      res <- bootstrap_pluck_summary(fit, min_contribution_threshold = min_contribution_threshold)
    else if(fit_metric == "error_and_cosine")
      res <- bootstrap_pluck_error_and_cosine(fit)
    else if(fit_metric == "error_and_cosine_bootstraps")
      res <- bootstrap_pluck_error_and_cosine_bootstraps(fit)
    else if(fit_metric == "p_val")
      res <- bootstrap_pluck_pval(fit)

    ls_res <- split(res, f = res[["sample"]])
    for (sample in names(ls_res)) {
      tmp_outfile=glue::glue("{output_dir}/{fit_type}_fit.{sample}.{ref}.{fit_metric}.csv")
      write_compressed_csv(
        x=dplyr::select(.data = ls_res[[sample]], -sample),
        tmp_outfile
      )
      cli::cli_alert_success("{fit_type} model fit [{fit_metric}] has been written to csv: {.path {tmp_outfile}.gz}")
    }
  }
}

extract_model_info <-  function(fit, ref, min_contribution_threshold){

  # Skip if fit is NULL (means there were no mutations of fit_type)
  if(is.null(fit)) {
    cli::cli_alert_warning("No model outputs produced for class {fit_type}")
    return(NULL)
  }

  fit_metrics <- c("expo", "expo_bootstraps", "bootstrap_summary", "error_and_cosine", "error_and_cosine_bootstraps",  "p_val")

  ls_model_outputs <- purrr::map(fit_metrics, function(fit_metric){
    if(fit_metric == "expo")
      res <- bootstrap_pluck_expo(fit, min_contribution_threshold = min_contribution_threshold)
    else if(fit_metric == "expo_bootstraps")
      res <- bootstrap_pluck_expo_bootstraps(fit)
    else if(fit_metric == "bootstrap_summary")
      res <- bootstrap_pluck_summary(fit, min_contribution_threshold = min_contribution_threshold)
    else if(fit_metric == "error_and_cosine")
      res <- bootstrap_pluck_error_and_cosine(fit)
    else if(fit_metric == "error_and_cosine_bootstraps")
      res <- bootstrap_pluck_error_and_cosine_bootstraps(fit)
    else if(fit_metric == "p_val")
      res <- bootstrap_pluck_pval(fit)

      ls_res <- split(res, f = res[["sample"]])
      return(ls_res)
    })

  names(ls_model_outputs) <- fit_metrics

  # Flip Inside Out so sample is the outer layer
  samples <- unique(fit$expo$sample)

  ls_sample_model_outputs <- lapply(samples, function(sample){
    ls_model_outputs <- lapply(fit_metrics, function(fit_metric){
      ls_model_outputs[[fit_metric]][[sample]]
      })

    names(ls_model_outputs) <- fit_metrics
    return(ls_model_outputs)
  })

  names(ls_sample_model_outputs) <- samples

  return(ls_sample_model_outputs)
}

# Add model_unfiltered, etc
enrich_model_info <- function(ls_model_info){

  lapply(ls_model_info, function(ls_sample_model_info){

    # Add unfiltered model in sigverse format (describing all non-zero contributions)
    expo <- ls_sample_model_info[["expo"]]
    expo_sorted <- dplyr::arrange(expo, dplyr::desc(contribution))
    model_unfiltered <- expo_sorted$contribution
    names(model_unfiltered) <- expo_sorted$signature
    ls_sample_model_info$model_unfiltered <- model_unfiltered[model_unfiltered > 0]


    return(ls_sample_model_info)
  })
}

#' Sigminer tally matrix -> a sigverse collection
#'
#' @param matrix sigminer tally matrix (must be transposed)
#' @param class signature class. Must be one of [sigstash::sig_valid_sigclass()]
#' @return sigverse catalogue collection. See [sigshared::example_catalogue_collection()]
#'
sigminer_tally_to_sigverse_catalogue_collection <- function(matrix, class){
  samples = colnames(matrix)
  df_longform_tally <- fix_tally(matrix)
  ls_longform <- split(df_longform_tally, f = df_longform_tally[["SampleID"]])

  ls_catalogue_collection <- lapply(ls_longform, function(df_tally){
    df_tally <- dplyr::select(df_tally, -SampleID)
    df_tally <- dplyr::rename(df_tally, channel=Context, count=Count, fraction=CountRelative)
    df_tally[["channel"]] <- sigstash::sig_convert_channel_name(df_tally[["channel"]], from = "sigminer", to = "cosmic")
    df_tally[["type"]] <- sigstash::sig_convert_channel2type(df_tally[["channel"]], sigclass = class)
    df_tally[["fraction"]] <- ifelse(is.na(df_tally[["fraction"]]), yes = 0, no = df_tally[["fraction"]])
    df_tally <- df_tally[c("channel", "type", "fraction", "count")]
    return(df_tally)
  })

  names(ls_catalogue_collection) <- names(ls_longform)

  return(ls_catalogue_collection)
}


write_tally_matrices <- function(catalogues, class, output_dir, ref){

  for (sample in names(catalogues)){
    outfile <- glue::glue("{output_dir}/{class}_catalogue.{sample}.{ref}.tally.csv.gz")
    write_compressed_csv(
      x = catalogues[[sample]],
      file = outfile
    )
  }
}

write_tally_matrix <- function(matrix, class, output_dir, ref){
  samples = colnames(matrix)

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

    # Project to existing UMAP
    if(is.null(ref_umaps_prefix)){
      cli::cli_alert_info("Skipping projection onto {class} reference umaps because {.arg ref_umaps_prefix} argument was not supplied")
      next
    }
    path_umap = paste0(ref_umaps_prefix, '.', class)
    if(!file.exists(path_umap)){
      cli::cli_alert_info("Skipping projection onto {class} reference umap since could not find file {path_umap}")
      next
    }

    # Read the umap reference
    rlang::check_installed("uwot", reason = "Creating UMAPs requires the 'uwot' pacakge to be installed. Please run install.packages('uwot') and restart R.")
    umap_model <- uwot::load_uwot(path_umap)

    # Convert catalogue to the right form
    df_tally_wide <- catalogue_to_wide(df_tally, class=class)

    # Check column order matches expected from umap
    assertions::assert_identical(
      colnames(df_tally_wide), umap_model$column_order,
      msg = "Failed to project {class} features onto umap for sample {sample} because channel order
        does not match what was used build the original UMAP"
    )

    # Project onto existing umap
    coords <- uwot::umap_transform(df_tally_wide, model=umap_model, seed = seed, batch = TRUE)

    # Prepare ref matrix dataframe
    df_coords_refmatrix <- as.data.frame(umap_model$embedding)
    df_coords_refmatrix[["sample"]] <- rownames(df_coords_refmatrix)
    df_coords_refmatrix <- df_coords_refmatrix[!df_coords_refmatrix$sample %in% sample,]
    rownames(df_coords_refmatrix) <- NULL

    # Prepare UMAP coordinate data.frame for sample of interest
    df_coords_sample <- as.data.frame(coords)
    df_coords_sample[["sample"]] <- sample

    # Combine the two
    df_coords <- rbind(df_coords_sample, df_coords_refmatrix)
    colnames(df_coords) <- c("dim1" , "dim2", "sample")

    # Add sample metadata
    # TODO: add sample_metadata arg (first add to reference matrix creation function)

    # Write the resulting UMAP dataset
    path_umap_output <- glue::glue("{output_dir}/{class}_umap.{sample}.{ref}.umap.csv.gz")
    write_compressed_csv(x = df_coords, file = path_umap_output)

    cli::cli_alert_success("UMAP written to csv : {.path {path_umap_output}.gz}")
  }
}
