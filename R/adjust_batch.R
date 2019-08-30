#' Zero-inflated empirical Bayes adjustment of batch effect in compositional
#' feature abundance data
#'
#' @param feature_abd feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param batch name of the batch variable. This variable in data should be a
#' factor variable and will be converted to so with a warning if otherwise.
#' @param covariates name(s) of covariates to adjust for
#' in the batch correction model.
#' @param data data frame of metadata, columns must include batch and covariates
#' (if specified).
#' @param zero_inflation logical. Whether or not a zero-inflated model should be
#' run. Default to TRUE (zero-inflated model). If set to FALSE then ComBat
#' (with parametric adjustment) will be performed.
#' @param pseudo_count numeric. Pseudo count to add feature_abd before
#' log-transformation. Automatically set to half of minimal non-zero value if
#' not specified.
#' @param diagnostics logical. Whether or not diagnostic plots are generated.
#' @param verbose logical. Whether or not verbose modelling information is
#' printed.
#' @param control list of control parameters. See details.
#'
#' @return feature-by-sample matrix of adjusted feature abundance.
#' @export
adjust_batch <- function(feature_abd,
                         batch,
                         covariates = NULL,
                         data,
                         zero_inflation = TRUE,
                         pseudo_count = NULL,
                         diagnostics = FALSE,
                         verbose = TRUE,
                         control) {
  # Check data formats
  # Check feature abundance table
  feature_abd <- as.matrix(feature_abd)
  type_feature_abd <- check_feature_abd(feature_abd = feature_abd)
  if(verbose)
    message("feature_abd is ", type_feature_abd)
  # Check metadata data frame
  if(length(batch) > 1)
    stop("Only one batch variable is supported!")
  data <- as.data.frame(data)
  samples <- check_samples(feature_abd = feature_abd,
                           data = data)
  # Check batch and covariates are included in metadata data frame
  df_batch <- check_metadata(data = data,
                             variables = batch)
  df_covariates <- check_metadata(data = data,
                                  variables = covariates)
  # Check batch variable
  df_batch[[batch]] <- check_batch(df_batch[[batch]], min_n_batch = 2)
  n_batch <- nlevels(x = df_batch[[batch]])
  if(verbose)
    message("Found ", n_batch, " batches")

  # check and construct controls
  control <- match_control(default = control_adjust_batch,
                           control = control)

  # Construct batch and covariate model matrices. Check for confounding
  batchmod <- construct_design(data = df_batch)
  mod <- construct_design(data = df_covariates)
  if(!check_rank(design = mod))
    stop("Covariates are confounded!")
  design <- cbind(batchmod, mod)
  if(!check_rank(design = design))
    stop("Covariates and batch are confounded!")
  if(verbose)
    message("Adjusting for ", ncol(mod),
            " covariate(s) or covariate(s) level(s)")

  # Transform data for ComBat fit
  if(is.null(pseudo_count)) {
    pseudo_count <- set_pseudo(features = feature_abd)
    if(verbose)
      message("Pseudo count is not specified and set to half of minimal ",
              "non-zero value: ",
              format(pseudo_count, digits = 3, scientific = TRUE))
  }
  log_data <- transform_features(
    features = normalize_features(
      features = feature_abd,
      normalization = "TSS",
      pseudo_count = pseudo_count),
    transform = "LOG")

  # Identify data to adjust for
  l_ind <- construct_ind(feature_abd = feature_abd,
                         n_batch = n_batch,
                         design = design,
                         zero_inflation = zero_inflation)
  if(verbose)
    message("Adjusting for (after filtering) ", sum(l_ind$ind_feature),
            " features")

  # Standardize data across features
  if(verbose)
    message("Standardizing data across features")
  stand_fit <- fit_stand_feature(s_data = log_data,
                                 design = design,
                                 l_ind = l_ind)
  s_data <- stand_fit$s_data
  l_stand_feature <- stand_fit$l_stand_feature

  # Estimate per-batch location and scale parameters
  # and EB hyper-parameters
  if(verbose)
    message("Estimating batch difference parameters and EB priors")
  params_fit <- fit_EB(s_data = s_data, l_stand_feature = l_stand_feature,
                       batchmod = batchmod, n_batch = n_batch,
                       l_ind = l_ind)

  # Shrink per-batch location and scale parameters
  if(verbose)
    message("Performing shrinkage adjustments on batch difference parameters")
  params_shrinked <- fit_shrink(s_data = s_data, l_params = params_fit,
                                batchmod = batchmod, n_batch = n_batch,
                                l_ind = l_ind,
                                control = control)

  # Adjust the data
  if(verbose)
    message("Performing batch corrections")
  adj_data <- adjust_EB(s_data = s_data, l_params_shrink = params_shrinked,
                        l_stand_feature = l_stand_feature,
                        batchmod = batchmod, n_batch = n_batch,
                        l_ind = l_ind)

  # Transform adjusted data back to the original scale
  # For debugging only, this shouldn't happen
  ## FIXME
  if(any(is.na(adj_data)))
    stop("There are missing values in the adjusted data!")
  feature_abd_adj <- back_transform_abd(adj_data = adj_data,
                                        feature_abd = feature_abd,
                                        type_feature_abd = type_feature_abd)


  # If required, generate diagnostic plots
  if(diagnostics)
    diagnostics_adjust_batch(feature_abd = feature_abd,
                             feature_abd_adj = feature_abd_adj,
                             batch = df_batch[[batch]],
                             gamma_hat = params_fit$gamma_hat,
                             gamma_star = params_shrinked$gamma_star)


  return(feature_abd_adj)
}

