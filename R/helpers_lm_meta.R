check_exposure <- function(exposure, batch){
  ind_exposure <- tapply(df_meta[[exposure]], df_batch[[batch]],
                         function(x) length(setdiff(unique(x), NA)) > 1)

}

# Check for exposure variables
if(any(!ind_exposure) & verbose)
  message("Exposure variable is missing or has only one non-missing value",
          " in the following batches; Maaslin2 won't be fitted on them:\n",
          paste(lvl_batch[!ind_exposure], collapse = ", "))

# Factor exposures must have common levels across batches
lvl.exposure <- NULL
if(is.character(data[, exposure, drop = TRUE]))
  data[, exposure] <- factor(data[, exposure, drop = TRUE])
if(is.factor(data[, exposure, drop = TRUE])) {
  lvl.exposure <- levels(data[, exposure, drop = TRUE])
  ind_exposure.cat <- tapply(data[, exposure, drop = TRUE], batch,
                             function(x) all(lvl.exposure %in% x))
  if(any(ind_exposure & !ind_exposure.cat))
    stop("Exposure is categorical and does not have common levels ",
         "in the following batches.\n",
         paste(lvl_batch[ind_exposure & !ind_exposure.cat], collapse = ", "))
}

#' Wrapper function for Maaslin2
#'
#' @param feature.abd feature*sample matrix of feature abundance.
#' @param data data frame of metadata.
#' @param exposure name of exposure variable.
#' @param covariates name of covariates.
#' @param covariates.random name of random covariates.
#' @param output directory for Maaslin2.
#' @param normalization normalization parameter for Maaslin2.
#' @param transform transformation parameter for Maaslin2.
#' @param analysis_method analysis method parameter for Maaslin2.
#'
#' @return a data frame recording per-feature coefficients, p-values, etc. from running Maaslin2.
Maaslin2.wrapper <- function(feature.abd,
                             data,
                             exposure,
                             covariates = NULL,
                             covariates.random = NULL,
                             output = "./",
                             normalization = "TSS",
                             transform = "AST",
                             analysis_method = "LM") {
  # Create temporary feature/sample/covariate names to avoid
  # Weird scenarios
  feature.abd.rename <- feature.abd
  data.rename <- data[, c(exposure, covariates, covariates.random), drop = FALSE]
  features.rename <- rename.Maaslin(rownames(feature.abd.rename), prefix = "T")
  samples.rename <- rename.Maaslin(colnames(feature.abd.rename), prefix = "S")
  exposure.rename <- rename.Maaslin(exposure, prefix = "E")
  covariates.rename <- rename.Maaslin(covariates, prefix = "X")
  covariates.random.rename <- rename.Maaslin(covariates.random, prefix = "RX")
  dimnames(feature.abd.rename) <- list(features.rename, samples.rename)
  dimnames(data.rename) <- list(samples.rename,
                                c(exposure.rename,
                                  covariates.rename,
                                  covariates.random.rename))
  # subset so that don't run into issues with all-zero features
  ind.feature <- apply(feature.abd.rename > 0, 1, any)

  # Run Maaslin2
  log.Maaslin <- suppressWarnings(
    capture.output(
      res.rename <- Maaslin2::Maaslin2(input_data = feature.abd.rename[ind.feature, , drop = TRUE],
                                       input_metadata = data.rename,
                                       output = output,
                                       min_abundance = 0,
                                       min_prevalence = 0,
                                       normalization = normalization,
                                       transform = transform,
                                       analysis_method = analysis_method,
                                       max_significance = 1,
                                       random_effects = covariates.random.rename,
                                       fixed_effects = c(exposure.rename, covariates.rename),
                                       standardize = FALSE,
                                       plot_heatmap = FALSE,
                                       plot_scatter = FALSE)$results
    ))
  # cat(paste(log.Maaslin, collapse = "\n"),
  #     file = file.path(output, "Maaslin2_warnings.log"))

  # Read Maaslin results
  lvl.exposure <- NULL
  if(is.factor(data[, exposure, drop = TRUE]))
    lvl.exposure <- levels(data[, exposure, drop = TRUE])
  suppressWarnings(table.Maaslin <-
                     dplyr::left_join(data.frame(feature = names(features.rename),
                                                 feature.rename = features.rename,
                                                 stringsAsFactors = FALSE),
                                      create.table.Maaslin(features.rename,
                                                           exposure.rename,
                                                           lvl.exposure),
                                      by = c("feature.rename" = "feature")))

  res <- dplyr::left_join(table.Maaslin, res.rename,
                          by = c("feature.rename" = "feature",
                                 "metadata",
                                 "value"))
  res <- dplyr::select(res, -feature.rename, -name)
  res$metadata <- exposure
  if(all(res$value == exposure.rename)) res$value <- exposure
  # Maaslin adjust p-values for all coefficients, modify to be for only the exposure
  res$qval <- p.adjust(res$pval, method = "fdr")

  return(res)
}

#' Utility for temporarily renaming samples/features for Maaslin2 run to bypass the rare
#' cases where unconventional names can cause exceptions
#'
#' @param old.names vector of names.
#' @param prefix prefix for the replacement (new numbered names).
#'
#' @return vector of new names - numbered vector with same length as old names and with the specified
#' prefix.
rename.Maaslin <- function(old.names, prefix) {
  if(is.null(old.names) | length(old.names) == 0) return(NULL)
  new.names <- paste0(prefix, seq_along(old.names))
  names(new.names) <- old.names
  return(new.names)
}

#' Utility for generating empty Maaslin2 results table
#'
#' @param features name of the features fitted to Maaslin2.
#' @param exposure the exposure variable.
#' @param lvl.exposure levels of the exposure variable, if a factor.
#'
#' @return a table for each feature-exposure value pai; reference level of exposure, if a factor,
#' is taken out because is absorbed into the intercept term in Maaslin2 regression.
create.table.Maaslin <- function(features, exposure, lvl.exposure) {
  if(is.null(lvl.exposure))
    values.exposure <- exposure
  else
    values.exposure <- lvl.exposure[-1]
  table.Maaslin <- expand.grid(features, exposure, values.exposure, stringsAsFactors = FALSE)
  names(table.Maaslin) <- c("feature", "metadata", "value")
  return(table.Maaslin)
}

#' Wrapper for fitting fixed/random effects meta-analysis model using metafor
#'
#' @param l.Maaslin.fit list of Maaslin2 result data frames, outputted from Maaslin2.wrapper.
#' @param method meta-analysis model to run, options provided in metafor::rma.
#' @param forest.plots should forest plots be generated (for the significant associations).
#' @param output directory for the output forest plots.
#' @param rma.threshold rma fit threshold control
#' @param rma.maxiter rma fit maximum iteration control
#'
#' @return a data frame recording per-feature meta-analysis association results.
#' (coefficients, p-values, etc.)
rma.wrapper <- function(l.Maaslin.fit, method = "REML",
                        forest.plots = TRUE, output,
                        rma.threshold = 1e-6,
                        rma.maxiter = 1000) {
  lvl.batch <- names(l.Maaslin.fit)
  n.batch <- length(lvl.batch)
  exposure <- unique(l.Maaslin.fit[[1]]$metadata)
  values.exposure <- unique(l.Maaslin.fit[[1]]$value)
  features <- unique(l.Maaslin.fit[[1]]$feature)
  l.results <- list()
  for(value.exposure in values.exposure) {
    i.result <- data.frame(matrix(NA,
                                  nrow = length(features),
                                  ncol = 11 + length(lvl.batch)))
    colnames(i.result) <- c("feature",
                            "exposure",
                            "coef",
                            "stderr",
                            "pval",
                            "k",
                            "tau2",
                            "stderr.tau2",
                            "pval.tau2",
                            "I2",
                            "H2",
                            paste0("weight_", lvl.batch))
    i.result$feature <- features
    i.result$exposure <- value.exposure
    rownames(i.result) <- i.result$feature
    if(forest.plots) pdf(paste0(output, exposure, "_", value.exposure, ".pdf"),
                         width = 6,
                         height = 4 + ifelse(n.batch > 4,
                                             (n.batch - 4) * 0.5,
                                             0))
    # sanity check
    if(any(features != l.Maaslin.fit[[2]][l.Maaslin.fit[[2]]$value == value.exposure, "feature"]))
      stop("Feature names don't match between l.Maaslin.fit components!")
    betas <- sapply(l.Maaslin.fit, function(i.Maaslin.fit) {
      i.Maaslin.fit[i.Maaslin.fit$value == value.exposure, "coef"]
    })
    sds <- sapply(l.Maaslin.fit, function(i.Maaslin.fit) {
      i.Maaslin.fit[i.Maaslin.fit$value == value.exposure, "stderr"]
    })
    pvals <- sapply(l.Maaslin.fit, function(i.Maaslin.fit) {
      i.Maaslin.fit[i.Maaslin.fit$value == value.exposure, "pval"]
    })
    rownames(betas) <- rownames(sds) <- rownames(pvals) <- features
    ind.feature <- !is.na(betas) & !is.na(sds) & (sds != 0)
    count.feature <- apply(ind.feature, 1, sum)
    for(feature in features) {
      if(count.feature[feature] >= 2) {
        tmp.rma.fit <- try(metafor::rma.uni(yi = betas[feature, ind.feature[feature, ]],
                                            sei = sds[feature, ind.feature[feature, ]],
                                            slab = lvl.batch[ind.feature[feature, ]],
                                            method = method,
                                            control = list(threshold = rma.threshold,
                                                           maxiter = rma.maxiter)),

                           silent = TRUE) # FIXME
        if("try-error" %in% class(tmp.rma.fit))
          next
        wts <- metafor::weights.rma.uni(tmp.rma.fit)
        i.result[feature, c("coef",
                            "stderr",
                            "pval",
                            "k",
                            "tau2",
                            "stderr.tau2",
                            "pval.tau2",
                            "I2",
                            "H2",
                            paste0("weight_",
                                   names(wts))
        )] <- c(unlist(tmp.rma.fit[c("beta",
                                     "se",
                                     "pval",
                                     "k",
                                     "tau2",
                                     "se.tau2",
                                     "QEp",
                                     "I2",
                                     "H2")]),
                wts)
        if(tmp.rma.fit$pval < 0.05 & forest.plots)
          metafor::forest(tmp.rma.fit,
                          xlab = shorten.name(feature, cutoff = 10),
                          slab = shorten.name(lvl.batch[ind.feature[feature, ]], cutoff = 6))
      }
      if(count.feature[feature] == 1) {
        tmp.ind.feature <- ind.feature[feature, ]
        tmp.batch <- lvl.batch[tmp.ind.feature]
        i.result[feature, c("coef",
                            "stderr",
                            "pval",
                            "k",
                            paste0("weight_",
                                   tmp.batch)
        )] <- c(betas[feature, tmp.ind.feature],
                sds[feature, tmp.ind.feature],
                pvals[feature, tmp.ind.feature],
                1,
                100)
      }
    }
    if(forest.plots) dev.off()
    i.result$pval.bonf <- p.adjust(i.result$pval, method = "bonf")
    i.result$qval.fdr <- p.adjust(i.result$pval, method = "fdr")

    l.results[[value.exposure]] <- i.result
  }
  results <- Reduce("rbind", l.results)
  return(results)
}

## This interface could be opened up in future versions?
#' Wrapper for fitting rma models with a moderator parameter. This allows to
#' analyze an interaction model with meta-analysis effects
#'
#' @param l.Maaslin.fit list of Maaslin2 result data frames, outputted from Maaslin2.wrapper.
#' @param data.moderator data frame recording the moderator variables. Each row corresponds to
#' a single study, and should have the same number of rows as l.Maaslin.fit.
#' @param method meta-analysis model to run, options provided in metafor::rma.
#' @param rma.threshold rma fit threshold control
#' @param rma.maxiter rma fit maximum iteration control
#'
#' @return a data frame recording per-feature/moderator value meta-analysis association results.
#' (coefficients, p-values, etc.)
rma.mod.wrapper <- function(l.Maaslin.fit,
                            data.moderator,
                            method = "REML",
                            rma.threshold = 1e-6,
                            rma.maxiter = 1000){
  lvl.batch <- names(l.Maaslin.fit)
  if(!all(lvl.batch %in% rownames(data.moderator)))
    stop("data.moderator must have all the batches fitted in Maaslin!")
  data.moderator <- data.moderator[lvl.batch, , drop = FALSE]
  exposure <- unique(l.Maaslin.fit[[1]]$metadata)
  values.exposure <- unique(l.Maaslin.fit[[1]]$value)
  features <- unique(l.Maaslin.fit[[1]]$feature)
  l.results <- list()
  for(value.exposure in values.exposure) {
    i.result <- data.frame(feature = features,
                           exposure = value.exposure,
                           tau2 = NA,
                           se.tau2 = NA,
                           p.tau2 = NA,
                           p.moderator = NA,
                           I2 = NA,
                           H2 = NA,
                           R2 = NA, stringsAsFactors = FALSE)
    rownames(i.result) <- i.result$feature
    # sanity check
    if(any(features != l.Maaslin.fit[[2]][l.Maaslin.fit[[2]]$value == value.exposure, "feature"]))
      stop("Feature names don't match between l.Maaslin.fit components!")
    betas <- sapply(l.Maaslin.fit, function(i.Maaslin.fit) {
      i.Maaslin.fit[i.Maaslin.fit$value == value.exposure, "coef"]
    })
    sds <- sapply(l.Maaslin.fit, function(i.Maaslin.fit) {
      i.Maaslin.fit[i.Maaslin.fit$value == value.exposure, "stderr"]
    })
    rownames(betas) <- rownames(sds) <- features
    ind.feature <- !is.na(betas) & !is.na(sds) & (sds != 0)
    count.feature <- apply(ind.feature, 1, sum)
    for(feature in features) {
      if(count.feature[feature] <= 1) next
      suppressWarnings(tmp.rma.fit <-
                         try(metafor::rma.uni(yi = betas[feature, ind.feature[feature, ]],
                                              sei = sds[feature, ind.feature[feature, ]],
                                              mod = ~.,
                                              data = data.moderator[ind.feature[feature, ], ,
                                                                    drop = FALSE],
                                              method = method,
                                              control = list(threshold = rma.threshold,
                                                             maxiter = rma.maxiter)),
                             silent = TRUE)) # FIXME
      if("try-error" %in% class(tmp.rma.fit))
        next
      if(is.null(tmp.rma.fit$R2))
        next
      i.result[feature, c("tau2",
                          "se.tau2",
                          "p.tau2",
                          "p.moderator",
                          "I2",
                          "H2",
                          "R2")] <-
        unlist(tmp.rma.fit[c("tau2",
                             "se.tau2",
                             "QEp",
                             "QMp",
                             "I2",
                             "H2",
                             "R2")])
    }
    l.results[[value.exposure]] <- i.result
  }
  results <- Reduce("rbind", l.results)
  results$R2[is.na(results$R2) & !is.na(results$tau2)] <- 0
  return(results)
}
