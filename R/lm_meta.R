#' Covariate adjusted meta-analytical differential abundance testing
#'
#' @param feature_abd feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param exposure name of the exposure variable for differential abundance
#' testing.
#' @param batch name of the batch variable. This variable in data should be a
#' factor variable and will be converted to so with a warning if otherwise.
#' @param covariates name(s) of covariates to adjust for in differential
#' abundance testing.
#' @param covariates_random name(s) of random effects grouping covariates to
#' adjust for in differential abundance testing.
#' @param data data frame of metadata, columns must include exposure, batch,
#' and covariates and covariates_random (if specified).
#' @param control list of control parameters. See details.
#'
#' @return a list with component meta.results for per-feature meta-analysis
#' results, and component l.Maaslin.fit which itself is a list of results from
#' fitting Maaslin2 in individual studies.
#' @export
lm_meta <- function(feature_abd,
                    exposure,
                    batch,
                    covariates = NULL,
                    covariates_random = NULL,
                    data,
                    control) {
  # Check and construct controls
  control <- match_control(default = control_adjust_batch,
                           control = control)
  verbose <- control$verbose

  # Check data formats
  # Check feature abundance table
  feature_abd <- as.matrix(feature_abd)
  type_feature_abd <- check_feature_abd(feature_abd = feature_abd)
  # Check metadata data frame
  if(length(batch) > 1)
    stop("Only one batch variable is supported!")
  data <- as.data.frame(data)
  samples <- check_samples(feature_abd = feature_abd,
                           data = data)
  # Check variables are included in metadata data frame
  df_batch <- check_metadata(data = data,
                             variables = batch)
  df_meta <- check_metadata(data = data,
                            variables = c(exposure,
                                          covariates,
                                          covariates_random))
  # Check batch variable
  df_batch[[batch]] <- check_batch(df_batch[[batch]], min_n_batch = 2)
  n_batch <- nlevels(df_batch[[batch]])
  lvl_batch <- levels(df_batch[[batch]])
  if(verbose)
    message("Found ", n_batch, " batches")



  # Determine which covariates can be fitted on which datasets.
  ind.covariate <- NULL
  if(!is.null(covariates)) {
    ind.covariate <- sapply(covariates, function(covariate) {
      tapply(data[, covariate, drop = TRUE], batch, function(x) {
        length(unique(x[!is.na(x)])) > 1
      })
    })
    for(covariate in covariates) {
      if(any(ind_exposure & !ind.covariate[, covariate]))
        message("Covariate ", covariate,
                " is missing or has only one non-missing value",
                " in the following batches; will be excluded from model for these",
                " batches:\n",
                paste(lvl_batch[ind_exposure & !ind.covariate[, covariate]],
                      collapse = ", "))
    }
  }

  # Check for random covariates
  ind.random <- NULL
  if(!is.null(covariates_random)) {
    if(length(covariates_random) > 1)
      stop("Multiple random covariates currently not supported.\n")
    ind.random <- sapply(covariates_random,
                         function(covariate.random) {
                           tapply(data[, covariate.random, drop = TRUE],
                                  batch, function(x) {
                                    length(unique(x[!is.na(x)])) > 1 & any(table(x) > 1)
                                  }) })
    if(!any(ind.random)) stop("Random covariates are provided ",
                              "but no batch has clustered observations!")
    message("Random covariates are provided, and will be fitted for ",
            "the following batches:\n",
            paste(lvl_batch[ind_exposure & apply(ind.random, 1, any)],
                  collapse = ", "))
  }

  # Create temporary output for Maaslin output files
  output.Maaslin <- paste0(output, "/Maaslin/")
  dir.create(output.Maaslin, recursive = TRUE)

  # Fit individual models
  l.Maaslin.fit <- list()
  for(i in 1:n.batch) {
    i.batch <- lvl_batch[i]
    if(!ind_exposure[i.batch]) next
    if(verbose) message("Fitting Maaslin2 on batch ", i.batch, "...")
    i.feature_abd <- feature_abd[, batch == i.batch]
    i.data <- data[batch == i.batch, ]
    i.covariates <- covariates[ind.covariate[i.batch, , drop = TRUE]]
    i.covariates_random <- covariates_random[ind.random[i.batch, , drop = TRUE]]
    i.output.Maaslin <- paste0(output.Maaslin, i.batch)
    dir.create(i.output.Maaslin)
    i.Maaslin.fit <- Maaslin2.wrapper(
      feature_abd = i.feature_abd,
      data = i.data,
      exposure = exposure,
      covariates = i.covariates,
      covariates_random = i.covariates_random,
      output = i.output.Maaslin,
      normalization = normalization,
      transform = transform,
      analysis_method = analysis_method
    )
    l.Maaslin.fit[[i.batch]] <- i.Maaslin.fit
    l.Maaslin.fit[[i.batch]]$batch <- i.batch
  }

  # Fit fixed/random effects models
  if(verbose) message("Fitting meta-analysis model.")
  meta.results <- rma.wrapper(l.Maaslin.fit, method = rma_method,
                              forest_plots = forest_plots, output = output,
                              rma_threshold = rma_threshold,
                              rma_maxiter = rma_maxiter)

  return(list(meta.results = meta.results, ind.results = l.Maaslin.fit))
}
