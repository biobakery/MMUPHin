#' maaslin3-based meta-analytical differential abundance testing
#' 
#' \code{maaslin_meta} runs differential abundance models on microbial profiles 
#' within individual studies/batches using maaslin3, and aggregates per-batch 
#' effect sizes with a meta-analysis fixed/random effects model. It takes as 
#' input a feature-by-sample microbial abundance table and the accompanying 
#' metadata data frame which should include the batch indicator variable, the 
#' main exposure variable for differential abundance testing, and optional 
#' covariates and random covariates. The function first runs 
#' \code{\link[maaslin3]{maaslin3}} models on the exposure in each batch. The 
#' per-batch effect sizes for both abundance (linear) and prevalence (logistic) 
#' components are then aggregated with \code{\link[metafor]{rma.uni}} and reported 
#' as output. Additional parameters for both packages can be provided through 
#' \code{control} (see details).
#'
#' \code{control} should be provided as a named list of the following components 
#' (can be a subset).
#' \describe{
#' \item{normalization}{
#' character. \code{normalization} parameter for maaslin3. See 
#' \code{\link[maaslin3]{maaslin3}} for details and allowed values. Default to 
#' \code{"TSS"} (total sum scaling).
#' }
#' \item{transform}{
#' character. \code{transform} parameter for maaslin3. See 
#' \code{\link[maaslin3]{maaslin3}} for details and allowed values. Default to 
#' \code{"LOG"} (log transformation).
#' }
#' \item{rma_method}{
#' character. \code{method} parameter for rma.uni. See
#' \code{\link[metafor]{rma.uni}} for details and allowed values. Default to 
#' \code{"REML"} (restricted maximum-likelihood estimator).
#' }
#' \item{output}{
#' character. Output directory for intermediate maaslin3 output and optional 
#' visualizations. Default to \code{"maaslin_meta_output"}.
#' }
#' \item{rma_conv}{
#' numeric. Convergence threshold for rma.uni (corresponds to 
#' \code{control$threshold}). See \code{\link[metafor]{rma.uni}} for details.
#' Default to 1e-4.
#' }
#' \item{rma_maxit}{
#' integer. Maximum number of iterations allowed for rma.uni (corresponds to 
#' \code{control$maxiter}). See \code{\link[metafor]{rma.uni}} for details.
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
#' @param covariates names of covariates to adjust for in maaslin3 
#' differential abundance testing models.
#' @param covariates_random names of random effects grouping covariates to 
#' adjust for in maaslin3 differential abundance testing models.
#' @param data data frame of metadata, columns must include exposure, batch, 
#' and covariates and covariates_random (if specified).
#' @param control a named list of additional control parameters. See details.
#'
#' @return a list, with the following components:
#' \describe{
#' \item{meta_fits_abundance}{
#' data frame of per-feature meta-analytical results for the abundance component, 
#' including columns for effect sizes, p-values and q-values, heterogeneity 
#' statistics such as \eqn{\tau^2} and \eqn{I^2}.
#' }
#' \item{meta_fits_prevalence}{
#' data frame of per-feature meta-analytical results for the prevalence component 
#' (logistic modeling).
#' }
#' \item{maaslin_fits}{
#' list of lists, each one corresponding to the fitted results of 
#' maaslin3 in an individual batch.
#' }
#' \item{control}{list of additional control parameters used in the function 
#' call.
#' }
#' }
#' @export
#' @author Siyuan Ma, \email{syma.research@@gmail.com}
#' @examples
#' # Assuming CRC_abd and CRC_meta are available
#' fit_meta <- maaslin_meta(feature_abd = CRC_abd,
#'                          exposure = "study_condition",
#'                          batch = "studyID",
#'                          covariates = c("gender", "age"),
#'                          data = CRC_meta)
maaslin_meta <- function(feature_abd,
                         batch,
                         exposure,
                         covariates = NULL,
                         covariates_random = NULL,
                         data,
                         control = list()) {
  # Check and construct controls
  control <- match_control(default = control_maaslin_meta, control)
  verbose <- control$verbose
  
  feature_abd <- as.matrix(feature_abd)
  data <- as.data.frame(data)
  samples <- check_samples(feature_abd = feature_abd, data = data)
  
  if(length(batch) > 1) stop("Only one batch variable is supported!")
  df_batch <- check_metadata(data = data, variables = batch)
  df_meta <- check_metadata(data = data,
                            variables = c(exposure, covariates, covariates_random),
                            no_missing = FALSE)
  
  var_batch <- check_batch(df_batch[[batch]], min_n_batch = 2)
  n_batch <- nlevels(var_batch)
  lvl_batch <- levels(var_batch)
  if(verbose) message("Found ", n_batch, " batches")
  
  if(is.character(df_meta[[exposure]]))
    df_meta[[exposure]] <- factor(df_meta[[exposure]], 
                                  levels = stringr::str_sort(unique(df_meta[[exposure]])))
  
  ind_exposure <- check_exposure(df_meta[[exposure]], var_batch)
  if(any(!ind_exposure))
    warning("Exposure variable is missing or constant in batches: ",
            paste(lvl_batch[!ind_exposure], collapse = ", "))
  
  # check for covariates
  ind_covariates <- check_covariates(df_meta[covariates], var_batch)
  for(covariate in covariates) {
    if(any(ind_exposure & !ind_covariates[, covariate]))
      warning("Covariate ", covariate, " constant in batches: ",
              paste(lvl_batch[ind_exposure & !ind_covariates[, covariate]], collapse = ", "))
  }
  
  # check for random covariates
  ind_covariates_random <- check_covariates_random(df_meta[covariates_random], var_batch)
  for(covariate in covariates_random) {
    if(!any(ind_exposure & ind_covariates_random[, covariate]))
      warning("Random covariate ", covariate, " has no clustered observations!")
  }
  
  dir.create(control$output, recursive = TRUE, showWarnings = FALSE)
  
  maaslin_fits <- list()
  for(i in seq_len(n_batch)) {
    i_batch <- lvl_batch[i]
    if(!ind_exposure[i_batch]) next
    if(verbose) message("Fitting maaslin3 on batch ", i_batch, "...")
    
    # subsetting
    i_covariates <- covariates[ind_covariates[i_batch, , drop = TRUE]]
    i_covariates_random <- covariates_random[ind_covariates_random[i_batch, , drop = TRUE]]
    i_output <- paste0(control$output, "/", i_batch)
    dir.create(i_output, showWarnings = FALSE)
    
    maaslin_fits[[i_batch]] <- Maaslin3_wrapper(
      feature_abd = feature_abd[, var_batch == i_batch, drop = FALSE],
      data = df_meta[var_batch == i_batch, , drop = FALSE],
      exposure = exposure,
      covariates = i_covariates,
      covariates_random = i_covariates_random,
      output = i_output,
      control = control
    )
    maaslin_fits[[i_batch]]$batch <- i_batch
  }
  
  if(verbose) message("Fitting meta-analysis models.")
  meta_abundance <- rma_wrapper(maaslin_fits, stream = "abundance", control = control)
  meta_prevalence <- rma_wrapper(maaslin_fits, stream = "prevalence", control = control)
  
  return(list(meta_fits_abundance = meta_abundance,
              meta_fits_prevalence = meta_prevalence,
              maaslin_fits = maaslin_fits,
              control = control))
}