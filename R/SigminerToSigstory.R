
#' Create Sigstory Inputs from Sigminer Outputs
#'
#'  This function takes the outputs of sigminerUtils and parses them into sigstory ready visualisations & analysies
#'
#' @param signature_folder folder containing signature data
#' @param rds_outfile file to store serialised signature results that can serve as sigstory input.
#' @param sparsity_pvalue  max p-value, below which we consider the signatures to be stable across bootstraps. See [sigstats::sig_compute_experimental_p_value()] for pvalue computation.
#' @param ref_tallies = path to refmatrix parquet files
#' @param ref_exposures = path to refmatrix parquet files
#' @param ref_metadata = path to refmatrix parquet files
#'
#' @return A list with all the information required to build a signature report
#' @export
#'
sigminer2sigstory <- function(signature_folder = "colo829_signature_results_with_refset/COLO829v003T",
                              rds_outfile = "sigstory.rds",
                              ref_tallies = NULL, #ref_tallies= "pcawg_reference_set/refmatrix.tally.parquet",
                              ref_exposures = NULL, #ref_exposures="pcawg_reference_set/refmatrix.exposures.parquet",
                              ref_metadata= NULL, #ref_metadata="pcawg_reference_set/",
                              sparsity_pvalue = 0.05){

  # Check Signature Folder Exists
  assertions::assert_directory_exists(signature_folder)

  # Get all files in the signature directory
  df_files = parse_sigminer_utils_outputs(signature_folder)

  sample = unique(na.omit(df_files[["sample"]]))
  nsamples = length(sample)
  assertions::assert(nsamples == 1, msg = "Found {nsamples} samples in the signature folder [{sample}]. Expected only 1. Please ensure there is only only 1 samples results in {.path {signature_folder}}")
  assertions::assert_no_missing(sample, msg = "empty sample names are not supported")


  sigclasses <- unique(na.omit(df_files[["sigclass"]]))
  sigclasses_with_expo <- df_files[df_files$type == "expo",][["sigclass"]] |>
    unlist() |>
    na.omit() |>
    unique()

  # Get Signature Collection Names
  signature_collection_path <- subset(df_files, filename == "signature_collections.csv", select=filepath, drop=TRUE)[1]
  assertions::assert(length(signature_collection_path) > 0, msg = "Failed to find signature_collections.csv file for sample: [{sample}]")
  ls_sig_collection_names <- read_sig_collection_names(signature_collection_path)



  # Get Thresholds Used
  thresholds_path <- subset(df_files, filename == "thresholds.csv", select=filepath, drop=TRUE)[1]
  assertions::assert(length(thresholds_path) > 0, msg = "Failed to find thresholds.csv file for sample: [{sample}]")
  ls_thresholds <- read_sig_thresholds(thresholds_path)

  result_tree <- list(
    sample = sample,
    collection_names = ls_sig_collection_names
    )

  for (sigclass in sigclasses_with_expo){
  # for (sigclass in "SBS96"){
    collection_name=ls_sig_collection_names[[sigclass]]
    tally = as.list(df_files[df_files$sigclass == sigclass & df_files$type == "tally", .drop = FALSE][1,])
    expo = as.list(df_files[df_files$sigclass == sigclass & df_files$type == "expo", .drop = FALSE][1,])
    similarity = as.list(df_files[df_files$sigclass == sigclass & df_files$type == "similarity", .drop = FALSE][1,])
    umap = as.list(df_files[df_files$sigclass == sigclass & df_files$type == "umap", .drop = FALSE][1,])
    bootstrap_summary = as.list(df_files[df_files$sigclass == sigclass & df_files$type == "bootstrap_summary", .drop = FALSE][1,])
    bootstraps = as.list(df_files[df_files$sigclass == sigclass & df_files$type == "expo_bootstraps", .drop = FALSE][1,])

    assertions::assert(length(tally) > 0, msg = "Failed to find tally file for sigclass [{sigclass}] (sample: {sample})")
    assertions::assert(length(expo) > 0, msg = "Failed to find expo file for sigclass [{sigclass}] (sample: {sample})")
    assertions::assert(length(bootstrap_summary) > 0, msg = "Failed to find bootstrap_summary file for sigclass [{sigclass}] (sample: {sample})")
    assertions::assert(length(bootstraps) > 0, msg = "Failed to find expo_bootstraps file for sigclass [{sigclass}] (sample: {sample})")

    df_tally <- read_tally(tally$filepath, convert_channels_to_cosmic = TRUE)
    df_exposures <- read_expo(expo$filepath)
    df_bootstraps <- read_bootstraps(bootstraps$filepath)
    df_bootstrap_summary <- read_bootstrap_summary(bootstrap_summary$filepath)
    df_similarity <- read_similarity(similarity$filepath)
    n_comparison_samples <- if(is.null(df_similarity)) 0 else nrow(df_similarity)
    df_umap <- read_umap(umap$filepath)

    # browser()

    # df_bootstrap_summary <- dplyr::rename(df_bootstrap_summary, signature=Sig)
    # df_bootstraps <- dplyr::rename(df_bootstraps, signature=Sig)
    # df_exposures <- dplyr::rename(df_exposures, signature=Sig)

    valid_sigs <- df_bootstrap_summary |>
      subset(experimental_pval < sparsity_pvalue, select=signature, drop = TRUE)
      # dplyr::filter(experimental_pval < sparsity_pvalue) |>
      # # dplyr::select(signature) |>
      # dplyr::pull(Sig)

    model <- sigminerUtils_expo_to_model(df_exposures, valid_sigs)

    total_mutations <- sum(df_tally$count)
    df_exposures_valid <- subset(df_exposures, signature %in% valid_sigs)
    explained_mutations <- sum(df_exposures_valid[["contribution_absolute"]])
    unexplained_mutations <- total_mutations - explained_mutations

    # Plot Catalogues
    gg_tally <- sigvis::sig_visualise(df_tally, title = paste0(sample, "(", sigclass ,")"), options =  sigvis::vis_options(fontsize_title = 11))


    # Plot Reconstructed Vs Observed
    df_reconstructed <- sigstats::sig_combine(
      signatures = sigstash::sig_load(collection_name, format = "sigstash"),
      model = model,
      format = "signature"
    )

    cosine_reconstructed_vs_tally <- sigstats::sig_cosine_similarity(signature1 = df_reconstructed, df_tally)

    # Create Signature Analysis Result (Lets just move this to sig_analyse_mutations and output RDS)
    # sig_analysis_result <- sigshared::signature_analysis_result(
    #   sample = sample,
    #   sigclass = sigclass, model = model,
    #   number_of_mutations = total_mutations,
    #   unexplained_mutations = unexplained_mutations,
    #   model_fit = cosine_reconstructed_vs_tally,
    #   catalogue = df_tally,
    #   # signatures = actual signature collections,
    #   signature_annotations,
    #   analysis_details = "sigminerUtils",
    #   )

    gg_reconstructed_vs_observed <-  sigvis::sig_visualise_compare_reconstructed_to_observed(
      signature = df_reconstructed,
      catalogue = df_tally,
      title = sigstats::sig_model_to_string(model, pair_sep = "  "),
      options = sigvis::vis_options(fontsize_title = 11),
    )

    # Bootstraps
    # Convert df_bootstraps to sigverse bootstrap format
    df_bootstraps <- sigshared::bselect(
      df_bootstraps,
      c(
        "signature",
        bootstrap = "Type",
        contribution_absolute = "contribution_absolute",
        contribution = "contribution"
      )
    )
    # Plot bootstraps
    gg_signature_stability <- sigvis::sig_visualise_bootstraps(
      bootstraps = df_bootstraps,
      min_contribution_threshold = ls_thresholds$min_contribution_threshold,
      pvalue = sparsity_pvalue,
      horizontal = TRUE
    )

    # Plot Exposure Dotplots of exposure comparisons
    gg_dotplot <- if(!is.null(ref_exposures))
      exposure_dotplots(
        sample = sample,
        model = model,
        df_exposures_valid = df_exposures_valid,
        ref_exposures = ref_exposures,
        sigclass = sigclass
      )
    else
      NULL

    # Plot 5 most similar samples
    ls_similar_sample_plots <- if(!is.null(ref_tallies) & !is.null(df_similarity))
      similar_sample_plots(sample = sample, df_similarity = df_similarity, ref_tallies = ref_tallies, sigclass = sigclass)
    else
      NULL

    # Plot Umap
    gg_umap <- if(!is.null(df_umap)){
      df_metadata <- sigshared::bselect(df_umap, columns = "sample")
      df_metadata[["colour"]] <- dplyr::case_when(
          df_metadata[["sample"]] == sample ~ sample,
          .default = "other",
        )
      umap_colour_pal <- c("#D55E00", "#999999")
      names(umap_colour_pal) <- c(sample, "other")

    sigvis::sig_visualise_dimred(df_umap, metadata = df_metadata, col_colour = "colour", xlab = "UMAP1", ylab = "UMAP2") +
      ggplot2::scale_colour_manual(values = umap_colour_pal)
    }
    else
      NULL

    # Add Results to Tree
    result_tree[[sigclass]] <- list(
      collection_name = collection_name,
      df_tally = df_tally,
      df_exposures = df_exposures,
      df_exposures_valid = df_exposures_valid,
      df_umap = df_umap,
      n_comparison_samples = n_comparison_samples, # Number of reference samples checked for similarity to the current sample
      df_similarity = df_similarity,
      model=model,
      total_mutations = total_mutations,
      unexplained_mutations = unexplained_mutations,
      proportion_of_unexplained_mutations = unexplained_mutations/total_mutations,
      cosine_reconstructed_vs_observed = cosine_reconstructed_vs_tally,
      df_bootstraps = df_bootstraps,
      df_bootstrap_summary = df_bootstrap_summary,
      gg_tally = gg_tally,
      gg_reconstructed_vs_observed = gg_reconstructed_vs_observed,
      gg_signature_stability = gg_signature_stability,
      gg_dotplot = gg_dotplot,
      ls_similar_sample_plots = ls_similar_sample_plots,
      gg_umap = gg_umap,
      sigclass = sigclass,
      fitting_method = "Sigminer QP"
    )
  }

  saveRDS(result_tree, rds_outfile)
  return(result_tree)
}

delim_column_to_list <- function(char){
  if(!any(grepl(pattern = "|", fixed = TRUE, x = char)))
    return(as.numeric(char))
 strsplit(char, split = "|", fixed = TRUE) |>
    lapply(as.numeric)
}



sigminerUtils_expo_to_model <- function(df, signatures){
  df <- subset(df, signature %in% signatures)
  vals = df[["contribution"]]
  names(vals) <- df[["signature"]]
  return(vals)
}

read_sig_thresholds <- function(filepath){
  df = read.csv(filepath)
  l = as.list(df[[2]])
  names(l) <- df[[1]]
  return(l)
}

read_sig_collection_names <- function(filepath){
  df = read.csv(filepath)
  l = as.list(df[[2]])
  names(l) <- df[[1]]
  return(l)
}

read_bootstrap_summary <- function(filepath){
  assertions::assert_file_exists(filepath)
  df <- read.csv(filepath, header = TRUE)
  df <- tibble::as_tibble(df)
  return(df)
}

read_bootstraps <- function(filepath){
  assertions::assert_file_exists(filepath)
  df <- read.csv(filepath, header = TRUE)
  df <- tibble::as_tibble(df)
  return(df)
}

read_similarity <- function(filepath){
  if(is.na(filepath)) return(NULL);
  assertions::assert_file_exists(filepath)
  df <- read.csv(filepath, header = TRUE)
  df <- tibble::as_tibble(df)
  return(df)
}

read_umap <- function(filepath){
  if(is.na(filepath)) return(NULL);
  assertions::assert_file_exists(filepath)
  df <- read.csv(filepath, header = TRUE)
  df <- tibble::as_tibble(df)
  return(df)
}

read_tally <- function(filepath, convert_channels_to_cosmic = TRUE){
  assertions::assert_file_exists(filepath)
  df = read.csv(filepath, header = TRUE)
  sigshared::assert_catalogue(df)
  if(convert_channels_to_cosmic){
    df[["channel"]] <- sigstash::sig_convert_channel_name(channel = df[["channel"]], from = "sigminer", to = "cosmic")
  }
  return(df)
}

read_expo <- function(filepath, drop_type_and_method = TRUE){
  assertions::assert_file_exists(filepath)
  df = read.csv(filepath, header = TRUE)
  if(drop_type_and_method) { df <- df[!names(df) %in% c("Type", "Method")] }
  return(df)
}

parse_sigminer_utils_outputs <- function(filepath, exclude_na_filepaths = TRUE){
  filepaths <- dir(filepath, full.names = TRUE,)

  df_files <- tibble::tibble(
    sample = filepaths_to_sample_id(filepaths),
    type = filepaths_to_type(filepaths),
    subtype = filepaths_to_subtype(filepaths, type=type),
    sigclass = filepaths_to_sigclass(filepaths, type=type),
    filename = basename(filepaths),
    filepath = filepaths
  )
  sample <- unique(na.omit(df_files[["sample"]]))

  df_files[["sample"]] <- ifelse(df_files[["filename"]] == "signature_collections.csv", paste0(sample, collapse = ","), df_files[["sample"]])
  df_files[["sample"]] <- ifelse(df_files[["filename"]] == "thresholds.csv", paste0(sample, collapse = ","), df_files[["sample"]])

  df_files <- df_files[!is.na(df_files$sample),]
  return(df_files)
}

filepaths_to_subtype <- function(filepaths, type){
  filenames <- basename(filepaths)
  first <- extract_nth_component_from_filename(filenames,  n = 1, sep = ".")
  subtype <- sub(x=first, pattern = "^*.*_", replacement = "")

  if(!is.null(type)){
    subtype[type %in% c('log')] <- NA
  }

  return(subtype)
}

filepaths_to_sigclass <- function(filepaths, type = NULL){
  filenames <- basename(filepaths)
  first <- extract_nth_component_from_filename(filenames,  n = 1, sep = ".")
  class <- sub(x=first, pattern = "_.*", replacement = "")

  if(!is.null(type)){
    class[type %in% c('log')] <- NA
  }

  return(class)
}
filepaths_to_type <- function(filepaths){
  filenames <- basename(filepaths)

  last_extensions <- extract_nth_component_from_filename(filenames,  n = -1, sep = ".")
  second_last_extensions = extract_nth_component_from_filename(filenames,  n = -2, sep = ".")
  third_last_extensions = extract_nth_component_from_filename(filenames,  n = -3, sep = ".")

  type <- Map(f = function(last, secondlast, thirdlast){
    if(last == "log") {return("log")}
    else return(thirdlast)
    },last_extensions, second_last_extensions, third_last_extensions)

  type <- unname(unlist(type))

  return(type)
}

filepaths_to_sample_id <- function(filepaths){
  filenames <- basename(filepaths)
  is_log <- extract_nth_component_from_filename(filenames, n = -1) == "log"
  is_collection <- extract_nth_component_from_filename(filenames, n = 1) == "signature_collections"
  is_threshold <- extract_nth_component_from_filename(filenames, n = 1) == "thresholds"

  sample_ids <- ifelse(is_log,
                      yes = extract_nth_component_from_filename(filenames, n = 1),
                      no =  extract_nth_component_from_filename(filenames, n = 2)
                      )

  sample_ids <- ifelse(is_collection | is_threshold, yes = NA, no = sample_ids)

  return(sample_ids)
}


#' Extract the N-th component from filenames
#'
#' Splits each filename by a specified separator and returns the N-th component.
#' If the N-th component doesn't exist, returns `NA`.
#' Supports negative values for `n`, where `n = -1` returns the last component,
#' `n = -2` returns the second-to-last, and so on.
#'
#' @param filenames A character vector of filenames.
#' @param n An integer specifying which component to extract (default is 1).
#'   Negative values of `n` extract components from the end of the filename.
#' @param sep A character string used to split the filenames (default is ".").
#' @param fixed Logical, whether to treat `sep` as a literal string (default is TRUE).
#'
#' @return A character vector containing the N-th component from each filename,
#'   or `NA` if unavailable.
#'
#' @examples
#' # Example: Extracting the second component from filenames
#' extract_nth_component_from_filename(c("file.name.extension", "another.file.txt"), n = 2)
#' # Returns: c("name", "file")
#'
#' # Example: Extracting the last component from filenames
#' extract_nth_component_from_filename(c("file.name.extension", "another.file.txt"), n = -1)
#' # Returns: c("extension", "txt")
#'
#' @keywords internal
#' @noRd
extract_nth_component_from_filename <- function(filenames, n = 1, sep = ".", fixed = TRUE) {

  # Split each filename by the separator
  components_list <- strsplit(filenames, split = sep, fixed = fixed)

  # Extract the nth component from each split filename
  nth_components <- lapply(components_list, function(components) {

    # Positive n: Extract the nth component (if available)
    if (n > 0 && length(components) >= n) {
      return(components[n])
    }

    # Negative n: Extract from the end (-1 = last, -2 = second last, etc.)
    else if (n < 0 && length(components) >= abs(n)) {
      return(components[length(components) + n + 1])
    }

    # If the requested component doesn't exist, return NA
    return(NA)
  })

  # Convert the list of components to a character vector and return
  return(unlist(nth_components))
}

exposure_dotplots <- function(sample, model, df_exposures_valid, ref_exposures, sigclass){
  ref_exposures_parquet <- read_ref_exposures(ref_exposures)
  sigs = names(model)
  df_ref_exposures = ref_exposures_parquet |>
    dplyr::filter(
      Sig %in% sigs,
      class %in% sigclass,
      sample != !!sample
    ) |>
    dplyr::collect()

  df_exposures_valid[["sample"]] <- sample
  df_exposures_valid[["class"]] <- sigclass
  df_exposures <- dplyr::bind_rows(df_exposures_valid, df_ref_exposures)
  df_exposures[["colour"]] <- as.character(df_exposures[["sample"]] == sample)

  # Split and plot
  ls_signature_exposures <- split(df_exposures, df_exposures$Sig)
  gg_dotplots <- lapply(ls_signature_exposures, \(df_expo){
    sigvis::sig_visualise_dotplot(df_expo, col_sample = "sample", col_contribution = "ContributionRelative",
                                  sort_by = "frequency_fill",
                                  col_fill = "colour",
                                  # col_colour = "colour",
                                  palette_fill = c("TRUE"="red", "FALSE" = "black"),
                                  palette_colour = c("TRUE"="red", "FALSE" = "grey"),show_legend = FALSE
                                ) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
    })
  return(gg_dotplots)
}

# similar_samples a named vector where names are sample IDs and values are cosine similarity
similar_sample_plots <- function(sample, df_similarity, ref_tallies, sigclass){
  df_similarity_top5 = dplyr::slice_max(df_similarity,n = 5, order_by = .data[[colnames(df_similarity)[2]]], with_ties = FALSE)
  ref_tallies_parquet <- read_ref_exposures(ref_tallies)

  df_ref_tallies = ref_tallies_parquet |>
    dplyr::filter(
      sample %in% df_similarity_top5$sample,
      class %in% sigclass,
      sample != !!sample
    ) |>
    dplyr::collect()



  # return(df_ref_tallies)
#   df_exposures_valid[["sample"]] <- sample
#   df_exposures_valid[["class"]] <- sigclass
#   df_exposures <- dplyr::bind_rows(df_exposures_valid, df_ref_exposures)
#   df_exposures[["colour"]] <- as.character(df_exposures[["sample"]] == sample)
#
  # Split and plot
  fmt_round_2digits <- function() {sigvis::fmt_round(digits = 2)}

  ls_similar_samples <- split(df_ref_tallies, df_ref_tallies$sample)
  ls_similar_sample_plots <- lapply(df_similarity_top5$sample, \(curr_sample){
    df_sim = ls_similar_samples[[curr_sample]]
    cosine_sim = df_similarity_top5[[2]][match(curr_sample, df_similarity_top5[["sample"]])]
    sigvis::sig_visualise_minified(df_sim, proportion = cosine_sim, format = fmt_round_2digits, heights = c(0.85, 0.15))
    # minisig(df_sim, prop_contribution = cosine_sim)
  })
  names(ls_similar_sample_plots) <- names(ls_similar_samples)
  return(ls_similar_sample_plots)
}

fmt_pct <- function(proportion, digits = 2){
  paste(round(proportion*100, digits = digits), "%")
}

read_ref_exposures <- function(path){
  arrow::open_dataset(path)
  #dplyr::glimpse(ref_exposures_parquet)
}
# Should be moved to sigvis
minisig <- function(signature, prop_contribution){
  gg_sigplot <- sigvis::sig_visualise(signature, class = "signature") +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x.bottom= ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y.left = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(linewidth = 1)
    )
  gg_sigplot

  dd <- data.frame(
    contribution = c(prop_contribution, 1-prop_contribution),
    measure = c("contribution", "all")
  )
  bar = ggplot2::ggplot(dd, ggplot2::aes(x=contribution, y="")) +
    ggplot2::geom_col(fill = c("maroon", "grey90")) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::annotate(geom = "text", x = 1, y="", label=fmt_pct(prop_contribution, 0), hjust = 1.3) +
    ggplot2::theme_void()

  patchwork::wrap_plots(gg_sigplot,bar ,ncol = 1, heights = c(0.95, 0.1))
}
