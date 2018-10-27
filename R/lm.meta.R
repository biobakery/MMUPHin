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
                    forest.plots = TRUE,
                    directory = "./MMUPHin_lm.meta/",
                    verbose = TRUE,
                    ...) {
  ### Ensure data formatts are as expected
  feature.count <- as.matrix(feature.count)
  if(any(feature.count < 0, na.rm = TRUE))
    stop("Found negative values in the feature table!")
  if(any(is.na(feature.count)))
    stop("Found missing values in the feature table!")
  data <- as.data.frame(data)
  if(!all(c(batch, exposure, covariates, covariates.random) %in% names(data)))
    stop("Batch/covariate variable not found in data.")
  batch <- data[, batch]
  batch <- as.factor(batch)
  if(any(is.na(batch))) stop("Found missing values in the batch variable!")

  ## Data dimensions need to agree with each other
  if(ncol(feature.count) != nrow(data))
    stop("Dimensions of feature table and metadata table do not agree!")

  ## Check that sample names agree between the feature and metadata table
  ## And assign row and column names if emppty
  ## Shouldn't happen if the data is imported through interface
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

  ## Identify batch variables
  n.batch <- length(unique(batch))
  if(n.batch < 2) stop("Batch variable has only one level!")
  if(verbose) message("Found ", n.batch, " batches")
  lvl.batch <- levels(batch)

  ## Check that factor exposures have common levels across studies
  if(is.character(data[, exposure]) || is.factor(data[, exposure])) {
    if(any(table(batch, data[, exposure]) == 0))
      stop("Exposure is factor and does not have common levels ",
           "across batches!")
  }

  ## Check for random covariates
  ind.random <- NULL
  if(!is.null(covariates.random)) {
    ind.random <- sapply(covariates.random,
                         function(covariate.random) {
                           sapply(lvl.batch, function(i.batch) {
                             anyDuplicated(data[batch == i.batch, covariate.random]) > 0
                           })
                         })
    if(!any(ind.random)) stop("Random covariates are provided ",
                              "but no batch has clustered observations!")
    message("Random covariates are provided, and will be fitted for ",
            "the following batches: ",
            paste(lvl.batch[apply(ind.random, 1, any)], collapse = ", "))
  }

  ## Create temporary directory for Maaslin output files
  directory.Maaslin <- paste0(directory, "/MMUPHin_Maaslin_tmp/")
  dir.create(directory.Maaslin, recursive = TRUE)

  ## Fit individual models
  if(verbose) message("Fitting individual regressions.")
  if(is.character(data[, exposure]) || is.factor(data[, exposure])) {
    if(any(table(batch, data[, exposure]) == 0))
      stop("Exposure is factor and does not have common levels ",
           "across batches!")
  }
  l.Maaslin.fit <- list()
  for(i in 1:n.batch) {
    i.batch <- lvl.batch[i]
    if(verbose) message("Fitting Maaslin2 on batch ", i.batch, "...")
    i.feature.count <- feature.count[, batch == i.batch]
    i.data <- data[batch == i.batch, ]
    i.covariates.random <- NULL
    if(!is.null(ind.random))
      if(any(ind.random[lvl.batch == i.batch, ]))
        i.covariates.random <- covariates.random[ind.random[lvl.batch == i.batch, ]]
    i.directory.Maaslin <- paste0(directory.Maaslin, i.batch)
    dir.create(i.directory.Maaslin)
    Maaslin.fit <- Maaslin2.wrapper(
      taxa = i.feature.count,
      metadata = i.data,
      variables = c(exposure, covariates),
      covariates.random = i.covariates.random,
      directory = i.directory.Maaslin,
      normalization = normalization,
      analysis_method = analysis_method
    )
    l.Maaslin.fit[[i]] <- Maaslin.fit[[exposure]]
    l.Maaslin.fit[[i]]$batch <- i.batch
  }
  ind.results <- Reduce("rbind", l.Maaslin.fit)
  ind.results <- dplyr::arrange(ind.results, Feature, Value, batch)

  # Fit fixed/random effects models
  if(verbose) message("Fitting meta-analysis model.")
  exposure.values <- unique(unlist(
    lapply(l.Maaslin.fit, function(i.fit) i.fit$Value)
  ))
  meta.results <- list()
  for(exposure.value in exposure.values) {
    i.result <- data.frame(Feature = rownames(feature.count),
                           Exposure = exposure.value,
                           Coefficient = NA,
                           Standard.error = NA,
                           P.val = NA,
                           tau2 = NA,
                           I2 = NA,
                           H2 = NA,
                           QEp = NA,
                           QMp = NA)
    rownames(i.result) <- i.result$Feature
    if(forest.plots) pdf(paste0(directory, exposure.value, ".pdf"),
                         width = 4,
                         height = 4 + ifelse(n.batch > 4,
                                             (n.batch - 4) * 0.5,
                                             0))
    betas <- sapply(l.Maaslin.fit, function(Maaslin.fit) {
      Maaslin.fit[Maaslin.fit$Value == exposure.value, "Coefficient"]
    })
    sds <- sapply(l.Maaslin.fit, function(Maaslin.fit) {
      Maaslin.fit[Maaslin.fit$Value == exposure.value, "stderr"]
    })
    rownames(betas) <- rownames(sds) <- rownames(feature.count)
    ind.feature <- apply(!is.na(betas) & !is.na(sds), 1, sum) >= 2
    for(feature in rownames(feature.count)) {
      if(ind.feature[feature]) {
        tmp.rma.fit <- metafor::rma.uni(yi = betas[feature, ],
                                        sei = sds[feature, ],
                                        slab = lvl.batch,
                                        ...)
        i.result[feature, c("Coefficient",
                            "Standard.error",
                            "P.val",
                            "tau2",
                            "I2",
                            "H2",
                            "QEp",
                            "QMp")] <- unlist(tmp.rma.fit[c("beta",
                                                            "se",
                                                            "pval",
                                                            "tau2",
                                                            "I2",
                                                            "H2",
                                                            "QEp",
                                                            "QMp")])
        if(tmp.rma.fit$pval < 0.05 & forest.plots)
          metafor::forest(tmp.rma.fit, xlab = feature)
      }
    }
    dev.off()
    i.result$P.bonf <- p.adjust(i.result$P.val, method = "bonf")
    i.result$Q.fdr <- p.adjust(i.result$P.val, method = "fdr")

    meta.results[[exposure.value]] <- i.result
  }

  return(list(meta.results = meta.results, ind.results = ind.results))
}
