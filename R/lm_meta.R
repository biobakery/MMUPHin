#' Covariate adjusted meta-analytical differential abundance testing
#' 
#' \code{lm_meta} runs differential abundance models on microbial profiles 
#' within individual studies/batches, and aggregates per-batch effect sizes with
#' a meta-analysis fixed/random effects model. It takes as input a 
#' feature-by-sample microbial abundance table and the accompanying meta data 
#' data frame which should includes the batch indicator variable, the main 
#' exposure variable for differential abundance testing, and optional covariates
#' and random covariates. The function first runs 
#' \code{\link[Maaslin2]{Maaslin2}} models on the exposure with optional 
#' covariates/random covariates in each batch. The per-batch effect sizes are 
#' then aggregated with \code{\link[metafor]{rma.uni}} and reported as output.
#' Additional parameters, including those for both 
#' \code{\link[Maaslin2]{Maaslin2}} and \code{\link[metafor]{rma.uni}} can be
#' provided through \code{control} (see details).
#'
#' \code{control} should be provided as a named list of the following components
#' (can be a subset).
#' \describe{
#' \item{normalization}{
#' character. \code{normalization} parameter for Maaslin2. See 
#' \code{\link[Maaslin2]{Maaslin2}} for details and allowed values. Default to 
#' \code{"TSS"} (total sum scaling).
#' }
#' \item{transform}{
#' character. \code{transform} parameter for Maaslin2. See 
#' \code{\link[Maaslin2]{Maaslin2}} for details and allowed values. Default to 
#' \code{"LOG"} (log transformation).
#' }
#' \item{analysis_method}{
#' character. \code{analysis_method} parameter for Maaslin2. See 
#' \code{\link[Maaslin2]{Maaslin2}} for details and allowed values. Default to 
#' \code{"LM"} (linear modeling).
#' }
#' \item{rma_method}{
#' character. \code{method} parameter for rma.uni. See
#' \code{\link[metafor]{rma.uni}} for details and allowed values. Default to 
#' \code{"REML"} (estricted maximum-likelihood estimator).
#' }
#' \item{output}{
#' character. Output directory for intermediate Maaslin2 output and the optional
#' forest plots. Default to \code{"MMUPHin_lm_meta"}.
#' }
#' \item{forest_plot}{
#' character. Suffix in the name for the generated forest plots visualizing
#' significant meta-analyitical differential abundance effects. Default to 
#' \code{"forest.pdf"}. Can be set to \code{NULL} in which case no output will 
#' be generated.
#' }
#' \item{rma_conv}{
#' numeric. Convergence threshold for rma.uni (corresponds to 
#' \code{control$threshold}. See \code{\link[metafor]{rma.uni}} for details.
#' Default to 1e-4.
#' }
#' \item{rma_maxit}{
#' integer. Maximum number of iterations allowed for rma.uni (corresponds to 
#' \code{control$maxiter}. See \code{\link[metafor]{rma.uni}} for details.
#' Default to 1000.
#' }
#' \item{verbose}{
#' logical. Indicates whether or not verbose information will be printed.
#' }
#' }
#' 
#' @param feature_abd feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param batch name of the batch variable. This variable in data should be a
#' factor variable and will be converted to so with a warning if otherwise.
#' @param exposure name of the exposure variable for differential abundance
#' testing.
#' @param covariates names of covariates to adjust for in Maaslin2 
#' differential abundance testing models.
#' @param covariates_random names of random effects grouping covariates to
#' adjust for in Maaslin2 differential abundance testing models.
#' @param data data frame of metadata, columns must include exposure, batch,
#' and covariates and covariates_random (if specified).
#' @param control a named list of additional control parameters. See details.
#'
#' @return a list, with the following components:
#' \describe{
#' \item{meta_fits}{
#' data frame of per-feature meta-analyitical differential abundance results,
#' including columns for effect sizes, p-values and q-values, heterogeneity
#' statistics such as \eqn{\tau^2} and \eqn{I^2}, as well as weights for 
#' individual batches. Many of these statistics are explained in detail in
#' \code{\link[metafor]{rma.uni}}.
#' }
#' \item{maaslin_fits}{
#' list of data frames, each one corresponding to the fitted results of 
#' Maaslin2 in a individual batch. See \code{\link[Maaslin2]{Maaslin2}} on 
#' details of these output.
#' }
#' \item{control}{list of additional control parameters used in the function
#' call.
#' }
#' }
#' @export
#' @author Siyuan Ma, \email{siyuanma@@g.harvard.edu}
#' @examples
#' data("CRC_abd", "CRC_meta")
#' fit_meta <- lm_meta(feature_abd = CRC_abd,
#'                     exposure = "study_condition",
#'                     batch = "studyID",
#'                     covariates = c("gender", "age"),
#'                     data = CRC_meta)$meta_fits
lm_meta <- function(feature_abd,
                    batch,
                    exposure,
                    covariates = NULL,
                    covariates_random = NULL,
                    data,
                    control) {
  # Check and construct controls
  control <- match_control(default = control_lm_meta,
                           control = control)
  verbose <- control$verbose
  
  # Check data formats
  # Check feature abundance table
  feature_abd <- as.matrix(feature_abd)
  type_feature_abd <- check_feature_abd(feature_abd = feature_abd)
  # Check metadata data frame
  data <- as.data.frame(data)
  samples <- check_samples(feature_abd = feature_abd,
                           data = data)
  # Check variables are included in metadata data frame
  if(length(batch) > 1)
    stop("Only one batch variable is supported!")
  df_batch <- check_metadata(data = data,
                             variables = batch)
  df_meta <- check_metadata(data = data,
                            variables = c(exposure,
                                          covariates,
                                          covariates_random),
                            no_missing = FALSE)
  # Check batch variable
  var_batch <- check_batch(df_batch[[batch]], min_n_batch = 2)
  n_batch <- nlevels(var_batch)
  lvl_batch <- levels(var_batch)
  if(verbose)
    message("Found ", n_batch, " batches")
  
  
  # Determine if exposure can be fitted on each batch
  # First if exposure is character change to factor
  if(is.character(df_meta[[exposure]]))
    df_meta[[exposure]] <- as.factor(df_meta[[exposure]])
  ind_exposure <- check_exposure(df_meta[[exposure]], var_batch)
  if(any(!ind_exposure))
    warning("Exposure variable is missing or has only one non-missing value",
            " in the following batches; Maaslin2 won't be fitted on them\n",
            paste(lvl_batch[!ind_exposure], collapse = ", "))
  
  # Determine if/which covariates can be fitted on each batch
  ind_covariates <- check_covariates(df_meta[covariates], var_batch)
  for(covariate in covariates) {
    if(any(ind_exposure & !ind_covariates[, covariate]))
      warning("Covariate ", covariate,
              " is missing or has only one non-missing value",
              " in the following batches; will be excluded from model for",
              " these batches:\n",
              paste(lvl_batch[ind_exposure & !ind_covariates[, covariate]],
                    collapse = ", "))
  }
  
  # Determine if/which random covariates can be fitted on each batch
  ind_covariates_random <- check_covariates_random(df_meta[covariates_random], 
                                                   var_batch)
  for(covariate in covariates_random) {
    if(!any(ind_exposure & ind_covariates_random[, covariate]))
      warning("Random covariate ", covariate, 
              " has no clustered observations!")
    else if(verbose)
      message("Random covariate ", covariate,
              "will be fitted for the following batches:\n",
              paste(lvl_batch[ind_exposure & 
                                ind_covariates_random[, covariate]],
                    collapse = ", "))
  }
  
  # Create temporary output for Maaslin output files
  dir.create(control$output, recursive = TRUE, showWarnings = FALSE)
  
  # Fit individual models
  maaslin_fits <- list()
  for(i in seq_len(n_batch)) {
    i_batch <- lvl_batch[i]
    if(!ind_exposure[i_batch]) next
    if(verbose) message("Fitting Maaslin2 on batch ", i_batch, "...")
    i_feature_abd <- feature_abd[, var_batch == i_batch, drop = FALSE]
    i_data <- df_meta[var_batch == i_batch, , drop = FALSE]
    i_covariates <- covariates[ind_covariates[i_batch, , drop = TRUE]]
    i_covariates_random <- covariates_random[
      ind_covariates_random[i_batch, , drop = TRUE]]
    i_output <- paste0(control$output, "/", i_batch)
    dir.create(i_output, showWarnings = FALSE)
    i_maaslin <- Maaslin2_wrapper(
      feature_abd = i_feature_abd,
      data = i_data,
      exposure = exposure,
      covariates = i_covariates,
      covariates_random = i_covariates_random,
      output = i_output,
      normalization = control$normalization,
      transform = control$transform,
      analysis_method = control$analysis_method
    )
    maaslin_fits[[i_batch]] <- i_maaslin
    maaslin_fits[[i_batch]]$batch <- i_batch
  }
  
  # Fit fixed/random effects models
  if(verbose) message("Fitting meta-analysis model.")
  meta_fits <- rma_wrapper(maaslin_fits, 
                           method = control$rma_method,
                           output = control$output,
                           forest_plot = control$forest_plot, 
                           rma_conv = control$rma_conv,
                           rma_maxit = control$rma_maxit,
                           verbose = verbose)
  
  return(list(meta_fits = meta_fits, 
              maaslin_fits = maaslin_fits,
              control = control))
}
