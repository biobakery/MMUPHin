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
#'
#' @return feature-by-sample matrix of adjusted feature abundance.
#' @export
#'
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
  type_feature_abd <- check_feature_abd(feature_abd)
  if(verbose)
    message("feature_abd is ", type_feature_abd)
  # Check metadata data frame
  if(length(batch) > 1)
    stop("Only one batch variable is supported!")
  data <- as.data.frame(data)
  samples <- check_samples(feature_abd, data)
  # Check batch and covariates are included in metadata data frame
  df_batch <- check_metadata(data, batch)
  df_covariates <- check_metadata(data, covariates)
  # Check batch variable
  df_batch[[batch]] <- check_batch(df_batch[[batch]])
  n_batch <- nlevels(df_batch[[batch]])
  if(verbose)
    message("Found ", n_batch, " batches")

  # check and construct controls
  control <- match_control(control_adjust_batch, control)

  # Construct batch and covariate model matrices. Check for confounding
  batchmod <- construct_design(df_batch)
  mod <- construct_design(df_covariates)
  if(qr(mod)$rank < ncol(mod))
    stop("Covariates are confounded!")
  design <- cbind(batchmod, mod)
  if(qr(design)$rank < ncol(design))
    stop("Covariates and batch are confounded!")
  if(verbose)
    message("Adjusting for ", ncol(mod),
            " covariate(s) or covariate(s) level(s)")

  # Transform data for ComBat fit
  if(is.null(pseudo_count)) {
    pseudo_count <- set_pseudo(feature_abd)
    if(verbose)
      message("Pseudo count is not specified and set to half of minimal ",
              "non-zero value: ",
              format(pseodu_count, digits = 3, scientific = TRUE))
  }
  log_data <- transform_features(
    normalize_features(
      feature_abd,
      normalization = "TSS",
      pseudo_count = pseudo_count),
    transform = "LOG")

  # Identify data to adjust for
  l_ind <- construct_ind(feature_abd, n_batch, design, zero_inflation)
  if(verbose)
    message("Adjusting for (after filtering) ", sum(l_ind$ind_feature),
            " features")

  # Standardize data across features
  if(verbose)
    message("Standardizing data across features.")
  stand_fit <- fit_stand_feature(s_data = log_data, design, l_ind)
  s_data <- stand_fit$s_data
  l_stand_feature <- stand_fit$l_stand_feature


  # Estimate per-batch location and scale parameters
  # and EB hyper-parameters
  if(verbose)
    message("Estimating batch difference parameters and EB priors")
  params_fit <- fit_EB(s_data, l_stand_feature = l_stand_feature,
                       batchmod, n_batch,
                       l_ind)

  # Shrink per-batch location and scale parameters
  if(verbose)
    message("Performing shrinkage adjustments on batch difference parameters")
  shrink_fit <- fit_shrink(s_data, l_params = params_fit,
                           batchmod, n_batch, l_ind,
                           control = control)

  # # Adjust the data
  # if(verbose) message("Performing batch corrections.")
  # adj.data <- s_data
  # for (i_batch in 1:n_batch){
  #   i_ind_feature <- !is.na(gamma_star[, i_batch]) & !is.na(delta_star[, i_batch])
  #   # For debugging, this shouldn't happen
  #   if(!all(i_ind_feature == ind_gamma[, i_batch])) stop("Something was wrong!") # FIXME
  #   for(i_feature in 1:nrow(adj.data)) {
  #     if(i_ind_feature[i_feature]) {
  #       i_ind.sample <- ind_data[i_feature, ] & as.logical(batchmod[, i_batch])
  #       adj.data[i_feature, i_ind.sample] <-
  #         (adj.data[i_feature, i_ind.sample] -
  #            gamma_star[i_feature, i_batch]) /
  #         sqrt(delta_star[i_feature, i_batch])
  #     }
  #   }
  # }
  #
  # for(i_feature in 1:nrow(adj.data)) {
  #   if(ind_feature[i_feature]) {
  #     standFit <- l_standFit[[i_feature]]
  #     adj.data[i_feature, ind_data[i_feature, ]] <-
  #       adj.data[i_feature, ind_data[i_feature, ]] *
  #       sqrt(standFit$var.pooled) +
  #       standFit$stand.mean
  #   }
  # }
  #
  # # For debugging only, this shouldn't happen
  # if(any(is.na(adj.data)))
  #   stop("There are missing values in the adjusted data!") # FIXME
  #
  # adj.data <- 2^adj.data
  # adj.data[feature_abd == 0] <- 0
  # adj.data <- normalizeFeatures(adj.data, normalization = "TSS")
  # feature_abd.adj <- t(t(adj.data) * read.depth)
  # dimnames(feature_abd.adj) <- dimnames(feature_abd)
  #
  # # If required, generate diagnostic plots
  # if(diagnostics)
  #   diagnostics.adjust.batch(feature_abd = feature_abd,
  #                            feature_abd.adj = feature_abd.adj,
  #                            batch = batch,
  #                            gamma_hat = gamma_hat,
  #                            gamma_star = gamma_star)
  #
  # return(feature_abd.adj)
}

