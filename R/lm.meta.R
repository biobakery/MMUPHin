#' Covariate-adjusted meta-analysis of per-feature associations in compositional data
#'
#' @param feature.abd feature*sample matrix of feature abundance (counts preferred).
#' @param exposure name of the exposure of interest variable.
#' @param batch  name of the batch variable.
#' @param covariates name(s) of additional covariates to adjust for in the Maaslin2 model.
#' @param covariates.random name(s) of random covariates in the Maaslin2 model.
#' @param data data frame of metadata, must contain exposure, batch, covariates (if present),
#' and covariates.random (if present).
#' @param normalization normalization parameter for Maaslin2.
#' @param transform transformation parameter for Maaslin2.
#' @param analysis_method analysis method parameter for Maaslin2.
#' @param rma.method method parameter for rma.
#' @param forest.plots should the function generate forest plots for significant results?
#' Deafault to TRUE.
#' @param verbose should verbose modelling information should be printed? Default to TRUE.
#' @param output output directory (for Maaslin2 output and forest plots).
#'
#' @return a list with component meta.results for per-feature meta-analysis results, and
#' component l.Maaslin.fit which itself is a list of results from fitting Maaslin2 in
#' individual studies.
#' @export
lm.meta <- function(feature.abd,
                    exposure,
                    batch,
                    covariates = NULL,
                    covariates.random = NULL,
                    data,
                    normalization = "TSS",
                    transform = "AST",
                    analysis_method = "LM",
                    rma.method = "REML",
                    forest.plots = TRUE,
                    output = "./MMUPHin_lm.meta/",
                    verbose = TRUE) {
  # Ensure data formats are as expected
  feature.abd <- as.matrix(feature.abd)
  if(any(is.na(feature.abd)))
    stop("Found missing values in the feature table!")
  if(any(feature.abd < 0))
    stop("Found negative values in the feature table!")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if(!all(c(batch, exposure, covariates, covariates.random) %in% names(data)))
    stop("Batch/covariate variable not found in data.")
  if(!all(apply(data[, c(batch, exposure, covariates, covariates.random),
                     drop = FALSE],
                2, class) %in% c("character", "numeric")))
    stop("Covariates must be of either character or numeric class!")

  # Data dimensions need to agree with each other
  if(ncol(feature.abd) != nrow(data))
    stop("Dimensions of feature table and metadata table do not agree!")

  # Check that sample names agree between the feature and metadata table
  # And assign row and column names if emppty
  if(is.null(colnames(feature.abd))) colnames(feature.abd) <-
    paste0("Sample",
           1:ncol(feature.abd))
  if(is.null(rownames(feature.abd))) rownames(feature.abd) <-
    paste0("Feature",
           1:nrow(feature.abd))
  if(is.null(rownames(data))) rownames(data) <-
    paste0("Sample",
           1:ncol(feature.abd))
  if(any(colnames(feature.abd) != rownames(data)))
    stop("Sample names in feature.abd and data don't agree!")

  # Check batch variable and identify groups
  batch <- data[, batch]
  if(any(is.na(batch))) stop("Found missing values in the batch variable!")
  batch <- as.factor(batch)
  n.batch <- length(unique(batch))
  if(n.batch < 2) stop("Batch variable has only one level!")
  if(verbose) message("Found ", n.batch, " batches")
  lvl.batch <- levels(batch)

  # Check for exposure variables
  ind.exposure <- tapply(data[, exposure, drop = TRUE], batch,
                         function(x) length(setdiff(unique(x), NA)) > 1)
  if(any(!ind.exposure) & verbose)
    message("Exposure variable is missing or has only one non-missing value",
            " in the following batches; Maaslin2 won't be fitted on them:\n",
            paste(lvl.batch[!ind.exposure], collapse = ", "))

  # Factor exposures must have common levels across batches
  lvl.exposure <- NULL
  if(is.character(data[, exposure, drop = TRUE]))
    data[, exposure] <- factor(data[, exposure, drop = TRUE])
  if(is.factor(data[, exposure, drop = TRUE])) {
    lvl.exposure <- levels(data[, exposure, drop = TRUE])
    ind.exposure.cat <- tapply(data[, exposure, drop = TRUE], batch,
                               function(x) all(lvl.exposure %in% x))
    if(any(ind.exposure & !ind.exposure.cat))
      stop("Exposure is categorical and does not have common levels ",
           "in the following batches.\n",
           paste(lvl.batch[ind.exposure & !ind.exposure.cat], collapse = ", "))
  }

  # Determine which covariates can be fitted on which datasets.
  ind.covariate <- NULL
  if(!is.null(covariates)) {
    ind.covariate <- sapply(covariates, function(covariate) {
      tapply(data[, covariate, drop = TRUE], batch, function(x) {
        length(unique(x[!is.na(x)])) > 1
      })
    })
    for(covariate in covariates) {
      if(any(ind.exposure & !ind.covariate[, covariate]))
        message("Covariate ", covariate,
                " is missing or has only one non-missing value",
                " in the following batches; will be excluded from model for these",
                " batches:\n",
                paste(lvl.batch[ind.exposure & !ind.covariate[, covariate]],
                      collapse = ", "))
    }
  }

  # Check for random covariates
  ind.random <- NULL
  if(!is.null(covariates.random)) {
    if(length(covariates.random) > 1)
      stop("Multiple random covariates currently not supported.\n")
    ind.random <- sapply(covariates.random,
                         function(covariate.random) {
                           tapply(data[, covariate.random, drop = TRUE],
                                  batch, function(x) {
                                    length(unique(x[!is.na(x)])) > 1 & any(table(x) > 1)
                                  }) })
    if(!any(ind.random)) stop("Random covariates are provided ",
                              "but no batch has clustered observations!")
    message("Random covariates are provided, and will be fitted for ",
            "the following batches:\n",
            paste(lvl.batch[ind.exposure & apply(ind.random, 1, any)],
                  collapse = ", "))
  }

  # Create temporary output for Maaslin output files
  output.Maaslin <- paste0(output, "/Maaslin/")
  dir.create(output.Maaslin, recursive = TRUE)

  # Fit individual models
  l.Maaslin.fit <- list()
  for(i in 1:n.batch) {
    i.batch <- lvl.batch[i]
    if(!ind.exposure[i.batch]) next
    if(verbose) message("Fitting Maaslin2 on batch ", i.batch, "...")
    i.feature.abd <- feature.abd[, batch == i.batch]
    i.data <- data[batch == i.batch, ]
    i.covariates <- covariates[ind.covariate[i.batch, , drop = TRUE]]
    i.covariates.random <- covariates.random[ind.random[i.batch, , drop = TRUE]]
    i.output.Maaslin <- paste0(output.Maaslin, i.batch)
    dir.create(i.output.Maaslin)
    i.Maaslin.fit <- Maaslin2.wrapper(
      feature.abd = i.feature.abd,
      data = i.data,
      exposure = exposure,
      covariates = i.covariates,
      covariates.random = i.covariates.random,
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
  meta.results <- rma.wrapper(l.Maaslin.fit, method = rma.method,
                              forest.plots = forest.plots, output = output)

  return(list(meta.results = meta.results, ind.results = l.Maaslin.fit))
}
