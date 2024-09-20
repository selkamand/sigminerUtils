
#' Create Sigstory Inputs from Sigminer Outputs
#'
#'  This function takes the outputs of sigminerUtils and parses them into sigstory ready visualisations & analysies
#'
#' @return A list with all the information required to build a signature report
#' @export
#'
sigminer2sigstory <- function(signature_folder = "colo829_signature_results/COLO829v003T", rds_outfile = "result_tree.rds", sparsity_pvalue = 0.05){
  assertions::assert_directory_exists(signature_folder)
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
    collection_name=ls_sig_collection_names[[sigclass]]
    tally = as.list(df_files[df_files$sigclass == sigclass & df_files$type == "tally", .drop = FALSE][1,])
    expo = as.list(df_files[df_files$sigclass == sigclass & df_files$type == "expo", .drop = FALSE][1,])
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


    valid_sigs <- df_bootstrap_summary |>
      subset(experimental_pval < sparsity_pvalue, select=Sig, drop = TRUE)

    model <- sigminerUtils_expo_to_model(df_exposures, valid_sigs)

    total_mutations <- sum(df_tally$count)
    df_exposures_valid <- subset(df_exposures, Sig %in% valid_sigs)
    explained_mutations <- sum(df_exposures_valid[["Contribution"]])
    unexplained_mutations <- total_mutations - explained_mutations


    # Plot Reconstructed Vs Observed
    df_reconstructed <- sigstats::sig_combine(
      signatures = sigstash::sig_load(collection_name, format = "sigstash"),
      model = model,
      format = "signature"
    )
    cosine_reconstructed_vs_tally <- sigstats::sig_cosine_similarity(signature1 = df_reconstructed, df_tally)

    gg_reconstructed_vs_observed <-  sigvis::sig_visualise_compare_reconstructed_to_observed(
      signature = df_reconstructed,
      catalogue = df_tally,
      title = sigstats::sig_model_to_string(model, pair_sep = "  "),
      options =sigvis::vis_options(fontsize_title = 11),
    )

    # Bootstraps
    # Convert df_bootstraps to sigverse bootstrap format
    df_bootstraps <- sigshared::bselect(
      df_bootstraps,
      c(
        signature = "Sig",
        bootstrap = "Type",
        contribution_absolute = "Contribution",
        contribution = "ContributionRelative"
      )
    )
    # Plot bootstraps
    gg_signature_stability <- sigvis::sig_visualise_bootstraps(
      bootstraps = df_bootstraps,
      min_contribution_threshold = ls_thresholds$min_contribution_threshold,
      pvalue = sparsity_pvalue,
      horizontal = TRUE
    )

    # Plot Dotplot
    browser()
    # Add Results to Tree
    result_tree[[sigclass]] <- list(
      collection_name = collection_name,
      df_tally = df_tally,
      df_exposures = df_exposures,
      df_exposures_valid = df_exposures_valid,
      model=model,
      total_mutations = total_mutations,
      proportion_of_unexplained_mutations = unexplained_mutations/total_mutations,
      cosine_reconstructed_vs_observed = cosine_reconstructed_vs_tally,
      df_bootstraps = df_bootstraps,
      df_bootstrap_summary = df_bootstrap_summary,
      gg_reconstructed_vs_observed = gg_reconstructed_vs_observed,
      gg_signature_stability = gg_signature_stability,
      fitting_method = "Sigminer QP"
    )


  }

  #browser()
  #saveRDS(result_tree, rds_outfile)
  return(result_tree)
}

delim_column_to_list <- function(char){
  if(!any(grepl(pattern = "|", fixed = TRUE, x = char)))
    return(as.numeric(char))
 strsplit(char, split = "|", fixed = TRUE) |>
    lapply(as.numeric)
}



sigminerUtils_expo_to_model <- function(df, signatures){
  df <- subset(df, Sig %in% signatures)
  vals = df[["ContributionRelative"]]
  names(vals) <- df[["Sig"]]
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
