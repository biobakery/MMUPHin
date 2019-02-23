#' Main function for batch effect adjustment
#'
#' @param feature.count Feature x sample matrix of feature abundance
#' @param batch Name of the batch variable
#' @param exposure Name of the primary exposure of interest variable
#' @param covariates Additional covariates for adjustment in individual regressions
#' @param covariates.random Group indicator for random effects
#' @param data Data frame for metadata.
#' @param normalization Normalization parameter for Maaslin2
#' @param transform Transformation parameter for Maaslin2
#' @param analysis_method Model parameter for Maaslin2
#' @param forest.plots Flag for whether or not forest plots figures should be generated.
#' Deafault to TRUE.
#' @param directory Directory for Maaslin output/generated forest plots
#' @param verbose Flag for whether or not verbose modelling information should be printed.
#' Default to yes.
#' @return a list
#' @export
#'
lm.meta <- function(feature.count,
                    batch,
                    exposure,
                    covariates = NULL,
                    covariates.random = NULL,
                    data,
                    normalization = "TSS",
                    transform = "AST",
                    analysis_method = "LM",
                    rma_method = "REML",
                    forest.plots = TRUE,
                    directory = "./MMUPHin_lm.meta/",
                    verbose = TRUE) {
  ## Ensure data formatts are as expected
  feature.count <- as.matrix(feature.count)
  if(any(is.na(feature.count)))
    stop("Found missing values in the feature table!")
  if(any(feature.count < 0))
    stop("Found negative values in the feature table!")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if(!all(c(batch, exposure, covariates, covariates.random) %in% names(data)))
    stop("Batch/covariate variable not found in data.")
  if(!all(apply(data[, c(batch, exposure, covariates, covariates.random),
                     drop = FALSE],
                2, class) %in% c("character", "numeric")))
    stop("Covariates must be of either character or numeric class!")

  ## Data dimensions need to agree with each other
  if(ncol(feature.count) != nrow(data))
    stop("Dimensions of feature table and metadata table do not agree!")

  ## Check that sample names agree between the feature and metadata table
  ## And assign row and column names if emppty
  if(is.null(colnames(feature.count))) colnames(feature.count) <-
    paste0("Sample",
           1:ncol(feature.count))
  if(is.null(rownames(feature.count))) rownames(feature.count) <-
    paste0("Feature",
           1:nrow(feature.count))
  if(is.null(rownames(data))) rownames(data) <-
    paste0("Sample",
           1:ncol(feature.count))
  if(any(colnames(feature.count) != rownames(data)))
    stop("Sample names in feature.count and data don't agree!")

  # Check batch variable and identify groups
  batch <- data[, batch]
  if(any(is.na(batch))) stop("Found missing values in the batch variable!")
  batch <- as.factor(batch)
  n.batch <- length(unique(batch))
  if(n.batch < 2) stop("Batch variable has only one level!")
  if(verbose) message("Found ", n.batch, " batches")
  lvl.batch <- levels(batch)

  ## Check for exposure variables
  ind.exposure <- tapply(data[, exposure, drop = TRUE], batch,
                         function(x) length(setdiff(unique(x), NA)) > 1)
  if(any(!ind.exposure) & verbose)
    message("Exposure variable is missing or has only one non-missing value",
            " in the following batches; Maaslin2 won't be fitted on them:\n",
            paste(lvl.batch[!ind.exposure], collapse = ", "))

  ## Factor exposures must have common levels across batches
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

  ## Determine which covariates can be fit on which datasets.
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

  ## Check for random covariates
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

  ## Create temporary directory for Maaslin output files
  directory.Maaslin <- paste0(directory, "/MMUPHin_Maaslin_tmp/")
  dir.create(directory.Maaslin, recursive = TRUE)

  ## Fit individual models
  l.Maaslin.fit <- list()
  for(i in 1:n.batch) {
    i.batch <- lvl.batch[i]
    if(!ind.exposure[i.batch]) next
    if(verbose) message("Fitting Maaslin2 on batch ", i.batch, "...")
    i.feature.count <- feature.count[, batch == i.batch]
    i.data <- data[batch == i.batch, ]
    i.covariates <- covariates[ind.covariate[i.batch, , drop = TRUE]]
    i.covariates.random <- covariates.random[ind.random[i.batch, , drop = TRUE]]
    i.directory.Maaslin <- paste0(directory.Maaslin, i.batch)
    dir.create(i.directory.Maaslin)
    i.Maaslin.fit <- Maaslin2.wrapper(
      feature.count = i.feature.count,
      data = i.data,
      exposure = exposure,
      covariates = i.covariates,
      covariates.random = i.covariates.random,
      directory = i.directory.Maaslin,
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
                              forest.plots = forest.plots, directory = directory)
  if(!is.null(cbind(ind.covariate, ind.random))) {
    meta.results.mod <- rma.mod.wrapper(l.Maaslin.fit, method = rma_method,
                                        data.moderator = cbind(ind.covariate, ind.random))
    meta.results <- dplyr::left_join(meta.results,
                                     meta.results.mod,
                                     by = c("feature", "exposure"),
                                     suffix = c("", "_moderator"))
  }
  return(list(meta.results = meta.results, ind.results = l.Maaslin.fit))
}
