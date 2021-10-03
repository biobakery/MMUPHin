#' Check exposure variable
#'
#' @param exposure exposure variable.
#' @param batch batch variable.
#'
#' @return vector of indicators per batch for whether or not the exposure can
#' be fitted within the batches
#' @keywords internal
check_exposure <- function(exposure, batch){
  ind_exposure <- as.vector(
    tapply(exposure, batch,
           function(x) length(setdiff(unique(x), NA)) > 1)
  )
  names(ind_exposure) <- levels(batch)
  
  # Factor exposures must have common levels across batches
  if(is.factor(exposure)) {
    lvl_exposure <- levels(exposure)
    ind_exposure_cat <- as.vector(
      tapply(exposure, batch,
             function(x) all(lvl_exposure %in% x))
    )
    if(any(ind_exposure & !ind_exposure_cat))
      stop("Exposure is character/factor and does not have common levels ",
           "in the following batches.\n",
           paste(lvl_batch[ind_exposure & !ind_exposure_cat], collapse = ", "))
  }
  
  return(ind_exposure)
}

#' Check covariates
#'
#' @param data_covariates data frame of covariates.
#' @param batch batch variable.
#'
#' @return vector of indicators per batch for if/which covariates can be fitted 
#' within the batches
#' @keywords internal
check_covariates <- function(data_covariates, batch){
  ind_covariates <- matrix(NA, 
                           nrow = nlevels(batch), 
                           ncol = ncol(data_covariates))
  dimnames(ind_covariates) <- list(levels(batch), names(data_covariates))
  ind_covariates[] <- vapply(
    data_covariates, 
    function(covariate)
      as.vector(tapply(
        covariate, batch, 
        function(x) {
          length(unique(x[!is.na(x)])) > 1})),
    rep_len(TRUE, nlevels(batch))
  )
  return(ind_covariates)
}

#' Check random covariates
#'
#' @param data_covariates data frame of random covariates.
#' @param batch batch variable.
#'
#' @return vector of indicators per batch for if/which random covariates can be 
#' fitted within the batches
#' @keywords internal
check_covariates_random <- function(data_covariates, batch){
  ind_covariates_random <- matrix(NA, 
                                  nrow = nlevels(batch), 
                                  ncol = ncol(data_covariates))
  dimnames(ind_covariates_random) <- list(levels(batch), names(data_covariates))
  ind_covariates_random[] <- vapply(
    data_covariates, 
    function(covariate)
      as.vector(tapply(
        covariate, batch, 
        function(x) {
          length(unique(x[!is.na(x)])) > 1 & any(table(x) > 1)
        })),
    rep_len(TRUE, nlevels(batch))
    )
  if(all(!ind_covariates_random) & ncol(ind_covariates_random) > 0) 
    stop("Random covariates are provided,",
         " but no batch has clustered observations!")
  return(ind_covariates_random)
}

#' Wrapper function for Maaslin2
#'
#' @param feature_abd feature*sample matrix of feature abundance.
#' @param data data frame of metadata.
#' @param exposure name of exposure variable.
#' @param covariates name of covariates.
#' @param covariates_random name of random covariates.
#' @param output directory for Maaslin2.
#' @param normalization normalization parameter for Maaslin2.
#' @param transform transformation parameter for Maaslin2.
#' @param analysis_method analysis method parameter for Maaslin2.
#'
#' @return a data frame recording per-feature coefficients, p-values, etc. from 
#' running Maaslin2.
#' @keywords internal
#' @importFrom utils capture.output
Maaslin2_wrapper <- function(feature_abd,
                             data,
                             exposure,
                             covariates = NULL,
                             covariates_random = NULL,
                             output = tempdir(),
                             normalization = "TSS",
                             transform = "AST",
                             analysis_method = "LM") {
  # Create temporary feature/sample/covariate names to avoid weird scenarios
  feature_abd_rename <- feature_abd
  data_rename <- data[, c(exposure, covariates, covariates_random), 
                      drop = FALSE]
  features_rename <- rename_maaslin(rownames(feature_abd_rename), prefix = "T")
  samples_rename <- rename_maaslin(colnames(feature_abd_rename), prefix = "S")
  exposure_rename <- rename_maaslin(exposure, prefix = "E")
  covariates_rename <- rename_maaslin(covariates, prefix = "X")
  covariates_random_rename <- rename_maaslin(covariates_random, prefix = "RX")
  dimnames(feature_abd_rename) <- list(features_rename, samples_rename)
  dimnames(data_rename) <- list(samples_rename,
                                c(exposure_rename,
                                  covariates_rename,
                                  covariates_random_rename))
  
  # subset so that don't run into issues with all-zero features
  ind_features <- apply(feature_abd_rename> 0, 1, any)
  
  # Run Maaslin2
  message_maaslin <- capture.output(
    log_maaslin <- catchToList(Maaslin2::Maaslin2(
      input_data = feature_abd_rename[ind_features, , drop = TRUE],
      input_metadata = data_rename,
      output = output,
      min_abundance = 0,
      min_prevalence = 0,
      normalization = normalization,
      transform = transform,
      analysis_method = analysis_method,
      max_significance = 1,
      random_effects = covariates_random_rename,
      fixed_effects = c(exposure_rename, covariates_rename),
      standardize = FALSE,
      plot_heatmap = FALSE,
      plot_scatter = FALSE)$results))
  
  if(!is.null(log_maaslin$error)) {
    ch_error <- log_maaslin$error
    # switch back to original variable names
    variables_rename <- c(exposure_rename, covariates_rename, covariates_random_rename)
    for(i_variable in names(variables_rename)) {
      i_pattern <- paste0("\'", variables_rename[i_variable], "\'")
      i_pattern_replace <- paste0("\'", i_variable, "\'")
      ch_error <- gsub(i_pattern, i_pattern_replace, x = ch_error,
                       fixed = TRUE)
    }
    
    stop("Internal Maaslin run error!\n", ch_error)
  }
    
  
  res_rename <- log_maaslin$value
  
  # Read Maaslin results
  lvl_exposure <- NULL
  if(is.factor(data[[exposure]]))
    lvl_exposure <- levels(data[[exposure]])
  table_maaslin <- dplyr::left_join(
    data.frame(feature = names(features_rename),
               feature_rename = features_rename,
               stringsAsFactors = FALSE),
    create_table_maaslin(features_rename,
                         exposure_rename,
                         lvl_exposure),
    by = c("feature_rename" = "feature"))
  
  res <- dplyr::left_join(table_maaslin, res_rename,
                          by = c("feature_rename" = "feature",
                                 "metadata",
                                 "value"))
  res <- dplyr::select(res, -feature_rename, -name)
  res$metadata <- exposure
  if(all(res$value == exposure_rename)) res$value <- exposure
  # Maaslin adjust p-values for all coefficients, modify to be for only the 
  # exposure
  res$qval <- p.adjust(res$pval, method = "fdr")
  
  return(res)
}

#' Utility for temporarily renaming samples/features for Maaslin2 run to bypass 
#' the rare cases where unconventional names can cause exceptions
#'
#' @param old_names vector of names.
#' @param prefix prefix for the replacement (new numbered names).
#'
#' @return vector of new names - numbered vector with same length as old names 
#' and with the specified prefix
#' @keywords internal
rename_maaslin <- function(old_names, prefix) {
  if(is.null(old_names) | length(old_names) == 0) return(NULL)
  new_names <- paste0(prefix, seq_along(old_names))
  names(new_names) <- old_names
  return(new_names)
}

#' Utility for generating empty Maaslin2 results table
#'
#' @param features name of the features fitted to Maaslin2.
#' @param exposure the exposure variable.
#' @param lvl_exposure levels of the exposure variable, if a factor.
#'
#' @return a table for each feature-exposure value pai; reference level of 
#' exposure, if a factor, is taken out because is absorbed into the intercept 
#' term in Maaslin2 regression
#' @keywords internal
create_table_maaslin <- function(features, exposure, lvl_exposure) {
  if(is.null(lvl_exposure))
    values_exposure <- exposure
  else
    values_exposure <- lvl_exposure[-1]
  names(features) <- NULL
  table_maaslin <- expand.grid(features, exposure, values_exposure, 
                               stringsAsFactors = FALSE)
  names(table_maaslin) <- c("feature", "metadata", "value")
  return(table_maaslin)
}

#' Wrapper for fitting fixed/random effects meta-analysis model using metafor
#'
#' @param maaslin_fits list of Maaslin2 result data frames, outputted from 
#' Maaslin2_wrapper.
#' @param method meta-analysis model to run, options provided in metafor::rma.
#' @param forest_plot logical. should forest plots be generated for 
#' the significant associations.
#' @param output directory for the output forest plots.
#' @param rma_conv rma threshold control.
#' @param rma_maxit rma maximum iteration control.
#' @param verbose should verbose information be printed.
#'
#' @return a data frame recording per-feature meta-analysis association results.
#' (coefficients, p-values, etc.)
#' @keywords internal
#' @importFrom grDevices dev.off pdf
#' @importFrom stats p.adjust
rma_wrapper <- function(maaslin_fits, 
                        method = "REML",
                        output = tempdir(),
                        forest_plot = NULL, 
                        rma_conv = 1e-6,
                        rma_maxit = 1000,
                        verbose = TRUE) {
  lvl_batch <- names(maaslin_fits)
  n_batch <- length(lvl_batch)
  exposure <- unique(maaslin_fits[[1]]$metadata)
  values_exposure <- unique(maaslin_fits[[1]]$value)
  features <- unique(maaslin_fits[[1]]$feature)
  l_results <- list()
  for(value_exposure in values_exposure) {
    i_result <- data.frame(matrix(NA,
                                  nrow = length(features),
                                  ncol = 11 + length(lvl_batch)))
    colnames(i_result) <- c("feature",
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
                            paste0("weight_", lvl_batch))
    i_result$feature <- features
    i_result$exposure <- value_exposure
    rownames(i_result) <- i_result$feature
    if(!is.null(forest_plot)) 
      pdf(paste0(output, "/",
                 exposure, "_", value_exposure, "_",
                 forest_plot),
          width = 6,
          height = 4 + ifelse(n_batch > 4,
                              (n_batch - 4) * 0.5,
                              0))
    # sanity check
    if(any(features != maaslin_fits[[2]][
      maaslin_fits[[2]]$value == value_exposure, 
      "feature"]))
      stop("Feature names don't match between maaslin_fits components!")
    betas <- vapply(
      maaslin_fits, 
      function(i_maaslin_fit)
        i_maaslin_fit[i_maaslin_fit$value == value_exposure, 
                      "coef", drop = TRUE],
      rep_len(0.0, length(features))
    )
    sds <- vapply(
      maaslin_fits, 
      function(i_maaslin_fit)
        i_maaslin_fit[i_maaslin_fit$value == value_exposure, 
                      "stderr", drop = TRUE],
      rep_len(0.0, length(features))
    )
    pvals <- vapply(
      maaslin_fits, 
      function(i_maaslin_fit)
        i_maaslin_fit[i_maaslin_fit$value == value_exposure, 
                      "pval", drop = TRUE],
      rep_len(0.0, length(features))
    )
    rownames(betas) <- rownames(sds) <- rownames(pvals) <- features
    ind_features <- !is.na(betas) & !is.na(sds) & (sds != 0)
    count_feature <- apply(ind_features, 1, sum)
    for(feature in features) {
      if(count_feature[feature] >= 2) {
        i_log <- catchToList(
          metafor::rma.uni(yi = betas[feature, ind_features[feature, ]],
                           sei = sds[feature, ind_features[feature, ]],
                           slab = lvl_batch[ind_features[feature, ]],
                           method = method,
                           control = list(threshold = rma_conv,
                                          maxiter = rma_maxit))
        )
        if(!is.null(i_log$error)) {
          warning("Fitting rma on feature ", feature, ";\n",
                  i_log$error)
          next
        }
        if(!is.null(i_log$warnings))
          warning("Fitting rma on feature ", feature, ";\n",
                  i_log$warnings)
        i_rma_fit <- i_log$value
        wts <- metafor::weights.rma.uni(i_rma_fit)
        i_result[feature, c("coef",
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
        )] <- c(unlist(i_rma_fit[c("beta",
                                   "se",
                                   "pval",
                                   "k",
                                   "tau2",
                                   "se.tau2",
                                   "QEp",
                                   "I2",
                                   "H2")]),
                wts)
        if(i_rma_fit$pval < 0.05 & !is.null(forest_plot))
          metafor::forest(
            i_rma_fit,
            xlab = shorten_name(feature, cutoff = 5),
            slab = shorten_name(lvl_batch[ind_features[feature, ]], cutoff = 3))
      }
      if(count_feature[feature] == 1) {
        i_ind_features <- ind_features[feature, ]
        tmp_batch <- lvl_batch[i_ind_features]
        i_result[feature, c("coef",
                            "stderr",
                            "pval",
                            "k",
                            paste0("weight_",
                                   tmp_batch)
        )] <- c(betas[feature, i_ind_features],
                sds[feature, i_ind_features],
                pvals[feature, i_ind_features],
                1,
                100)
      }
    }
    if(!is.null(forest_plot)) dev.off()
    i_result$pval.bonf <- p.adjust(i_result$pval, method = "bonf")
    i_result$qval.fdr <- p.adjust(i_result$pval, method = "fdr")
    
    l_results[[value_exposure]] <- i_result
  }
  results <- Reduce("rbind", l_results)
  return(results)
}

## This interface could be opened up in future versions?
#' Wrapper for fitting rma models with a moderator parameter. This allows to
#' analyze an interaction model with meta-analysis effects
#'
#' @param maaslin_fits list of Maaslin2 result data frames, outputted from Maaslin2_wrapper.
#' @param data.moderator data frame recording the moderator variables. Each row corresponds to
#' a single study, and should have the same number of rows as maaslin_fits.
#' @param method meta-analysis model to run, options provided in metafor::rma.
#' @param rma_conv rma fit threshold control
#' @param rma_maxit rma fit maximum iteration control
#'
#' @return a data frame recording per-feature/moderator value meta-analysis association results.
#' (coefficients, p-values, etc.)
#' @keywords internal
# rma.mod.wrapper <- function(maaslin_fits,
#                             data.moderator,
#                             method = "REML",
#                             rma_conv = 1e-6,
#                             rma_maxit = 1000){
#   lvl_batch <- names(maaslin_fits)
#   if(!all(lvl_batch %in% rownames(data.moderator)))
#     stop("data.moderator must have all the batches fitted in Maaslin!")
#   data.moderator <- data.moderator[lvl_batch, , drop = FALSE]
#   exposure <- unique(maaslin_fits[[1]]$metadata)
#   values_exposure <- unique(maaslin_fits[[1]]$value)
#   features <- unique(maaslin_fits[[1]]$feature)
#   l_results <- list()
#   for(value_exposure in values_exposure) {
#     i_result <- data.frame(feature = features,
#                            exposure = value_exposure,
#                            tau2 = NA,
#                            se.tau2 = NA,
#                            p.tau2 = NA,
#                            p.moderator = NA,
#                            I2 = NA,
#                            H2 = NA,
#                            R2 = NA, stringsAsFactors = FALSE)
#     rownames(i_result) <- i_result$feature
#     # sanity check
#     if(any(features != maaslin_fits[[2]][maaslin_fits[[2]]$value == value_exposure, "feature"]))
#       stop("Feature names don't match between maaslin_fits components!")
#     betas <- sapply(maaslin_fits, function(i_maaslin_fit) {
#       i_maaslin_fit[i_maaslin_fit$value == value_exposure, "coef"]
#     })
#     sds <- sapply(maaslin_fits, function(i_maaslin_fit) {
#       i_maaslin_fit[i_maaslin_fit$value == value_exposure, "stderr"]
#     })
#     rownames(betas) <- rownames(sds) <- features
#     ind_features <- !is.na(betas) & !is.na(sds) & (sds != 0)
#     count_feature <- apply(ind_features, 1, sum)
#     for(feature in features) {
#       if(count_feature[feature] <= 1) next
#       suppressWarnings(i_rma_fit <-
#                          try(metafor::rma.uni(yi = betas[feature, ind_features[feature, ]],
#                                               sei = sds[feature, ind_features[feature, ]],
#                                               mod = ~.,
#                                               data = data.moderator[ind_features[feature, ], ,
#                                                                     drop = FALSE],
#                                               method = method,
#                                               control = list(threshold = rma_conv,
#                                                              maxiter = rma_maxit)),
#                              silent = TRUE)) # FIXME
#       if("try-error" %in% class(i_rma_fit))
#         next
#       if(is.null(i_rma_fit$R2))
#         next
#       i_result[feature, c("tau2",
#                           "se.tau2",
#                           "p.tau2",
#                           "p.moderator",
#                           "I2",
#                           "H2",
#                           "R2")] <-
#         unlist(i_rma_fit[c("tau2",
#                            "se.tau2",
#                            "QEp",
#                            "QMp",
#                            "I2",
#                            "H2",
#                            "R2")])
#     }
#     l_results[[value_exposure]] <- i_result
#   }
#   results <- Reduce("rbind", l_results)
#   results$R2[is.na(results$R2) & !is.na(results$tau2)] <- 0
#   return(results)
# }
