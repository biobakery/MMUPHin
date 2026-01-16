#' Check exposure variable consistency across batches
#'
#' @param exposure Vector of the exposure variable.
#' @param batch Factor vector indicating batch/study membership.
#'
#' @return A logical vector indicating if the exposure has variance (>1 level) 
#'         within each batch.
#' @keywords internal
check_exposure <- function(exposure, batch) {
  ind_exposure <- as.vector(
    tapply(exposure, batch, function(x) length(setdiff(unique(x), NA)) > 1)
  )
  names(ind_exposure) <- levels(batch)
  
  # Factor exposures must have common levels across batches to allow meta-analysis
  if (is.factor(exposure)) {
    lvl_exposure <- levels(exposure)
    ind_exposure_cat <- as.vector(
      tapply(exposure, batch, function(x) all(lvl_exposure %in% x))
    )
    if (any(ind_exposure & !ind_exposure_cat)) {
      stop("Exposure is character/factor and does not have common levels ",
           "in the following batches:\n",
           paste(levels(batch)[ind_exposure & !ind_exposure_cat], collapse = ", "))
    }
  }
  
  return(ind_exposure)
}

#' Check fixed-effect covariates for variance within batches
#'
#' @param data_covariates Data frame containing only covariate columns.
#' @param batch Factor vector of batch IDs.
#'
#' @return A logical matrix (batch x covariate) where TRUE indicates the 
#'         covariate can be fitted in that batch.
#' @keywords internal
check_covariates <- function(data_covariates, batch) {
  ind_covariates <- matrix(NA, nrow = nlevels(batch), ncol = ncol(data_covariates))
  dimnames(ind_covariates) <- list(levels(batch), names(data_covariates))
  ind_covariates[] <- vapply(
    data_covariates, 
    function(covariate) {
      as.vector(tapply(covariate, batch, function(x) length(unique(x[!is.na(x)])) > 1))
    },
    rep_len(TRUE, nlevels(batch))
  )
  return(ind_covariates)
}

#' Check random-effect covariates for clustering within batches
#'
#' @param data_covariates Data frame of random effect variables.
#' @param batch Factor vector of batch IDs.
#'
#' @return A logical matrix indicating if random effects have clustered 
#'         observations within each batch.
#' @keywords internal
check_covariates_random <- function(data_covariates, batch) {
  ind_covariates_random <- matrix(NA, nrow = nlevels(batch), ncol = ncol(data_covariates))
  dimnames(ind_covariates_random) <- list(levels(batch), names(data_covariates))
  ind_covariates_random[] <- vapply(
    data_covariates, 
    function(covariate) {
      as.vector(tapply(covariate, batch, function(x) {
        length(unique(x[!is.na(x)])) > 1 & any(table(x) > 1)
      }))
    },
    rep_len(TRUE, nlevels(batch))
  )
  
  if (all(!ind_covariates_random) && ncol(ind_covariates_random) > 0) {
    stop("Random covariates are provided, but no batch has clustered observations!")
  }
  return(ind_covariates_random)
}

#' Internal wrapper for Maaslin 3 to ensure name parity and handle errors
#'
#' @param feature_abd Feature-by-sample matrix.
#' @param data Metadata data frame.
#' @param exposure Character; primary exposure variable.
#' @param covariates Character vector; fixed effects.
#' @param covariates_random Character vector; random effects.
#' @param control List of maaslin3 parameters.
#'
#' @return The standard Maaslin 3 fit object with restored original names.
#' @keywords internal
Maaslin3_wrapper <- function(feature_abd, data, exposure, covariates = NULL, 
                             covariates_random = NULL, 
                             output = tempdir(),
                             control) {
  
  # Sanitize names to prevent MaAsLin3 formula parsing errors
  feature_abd_rename <- feature_abd
  data_rename <- data[, c(exposure, covariates, covariates_random), drop = FALSE]
  features_rename <- rename_maaslin(rownames(feature_abd_rename), prefix = "T")
  samples_rename <- rename_maaslin(colnames(feature_abd_rename), prefix = "S")
  exposure_rename <- rename_maaslin(exposure, prefix = "E")
  covariates_rename <- rename_maaslin(covariates, prefix = "X")
  covariates_random_rename <- rename_maaslin(covariates_random, prefix = "RX")
  
  dimnames(feature_abd_rename) <- list(features_rename, samples_rename)
  dimnames(data_rename) <- list(samples_rename, 
                                c(exposure_rename, covariates_rename, covariates_random_rename))
  
  ind_features <- apply(feature_abd_rename > 0, 1, any)
  
  # Run MaAsLin3 while catching potential internal errors
  message_maaslin <- capture.output(
    log_maaslin <- catchToList(
      maaslin3::maaslin3(
        input_data = feature_abd_rename[ind_features, , drop = FALSE],
        input_metadata = data_rename, 
        output = output,
        min_abundance = as.numeric(control$min_abundance), 
        min_prevalence = as.numeric(control$min_prevalence),
        median_comparison_abundance = TRUE,
        median_comparison_prevalence = FALSE,
        normalization = control$normalization, 
        transform = control$transform,
        random_effects = covariates_random_rename, 
        fixed_effects = c(exposure_rename, covariates_rename),
        standardize = FALSE, 
        subtract_median = TRUE, 
        augment = control$augment, 
        verbosity = "ERROR",
        plot_associations = FALSE, 
        plot_summary_plot = FALSE, 
        cores = control$cores
      )
    )
  )
  
  # Your existing error handling logic
  if(!is.null(log_maaslin$error)) {
    ch_error <- log_maaslin$error
    # Switch back to original variable names for the error message
    variables_rename <- c(exposure_rename, covariates_rename, covariates_random_rename)
    for(i_variable in names(variables_rename)) {
      i_pattern <- paste0("\'", variables_rename[i_variable], "\'")
      i_pattern_replace <- paste0("\'", i_variable, "\'")
      ch_error <- gsub(i_pattern, i_pattern_replace, x = ch_error, fixed = TRUE)
    }
    
    stop("Internal Maaslin run error!\n", ch_error)
  }
  
  fit <- log_maaslin$value
  lvl_exposure <- if (is.factor(data[[exposure]])) levels(data[[exposure]]) else NULL
  
  # Internal helper to map back to original feature and metadata names
  restore_res <- function(res_df) {
    if (is.null(res_df)) return(NULL)
    table_maaslin <- dplyr::left_join(
      data.frame(feature = names(features_rename), feature_rename = features_rename, stringsAsFactors = FALSE),
      create_table_maaslin(features_rename, exposure_rename, lvl_exposure),
      by = c("feature_rename" = "feature"))
    
    res <- dplyr::left_join(table_maaslin, res_df, 
                            by = c("feature_rename" = "feature", "metadata", "value"))
    res <- dplyr::select(res, -feature_rename)
    res$metadata <- exposure
    
    return(res)
  }
  
  fit$fit_data_abundance$results <- restore_res(fit$fit_data_abundance$results)
  fit$fit_data_prevalence$results <- restore_res(fit$fit_data_prevalence$results)
  return(fit)
}

#' Utility to shorten names for plotting
#' @keywords internal
shorten_name <- function(x, cutoff = 15) {
  ifelse(nchar(x) > cutoff, paste0(substr(x, 1, cutoff), ".."), x)
}

#' Wrapper for fitting fixed/random effects meta-analysis model using metafor
#'
#' @param maaslin_fits list of results from individual batches.
#' @param stream stream to aggregate ("abundance" or "prevalence").
#' @param control control parameters including output path and plot settings.
#'
#' @return a data frame recording per-feature meta-analysis association results.
#' @keywords internal
rma_wrapper <- function(maaslin_fits, stream = "abundance", control) {
  # Extract correct results list based on stream
  fits_list <- lapply(maaslin_fits, function(x) {
    if(stream == "abundance") x$fit_data_abundance$results else x$fit_data_prevalence$results
  })
  
  fits_list <- fits_list[!vapply(fits_list, is.null, logical(1))]
  if (length(fits_list) == 0) return(NULL)
  
  lvl_batch <- names(fits_list)
  n_batch   <- length(lvl_batch)
  features  <- unique(unlist(lapply(fits_list, function(x) x$feature)))
  exp_vals  <- unique(unlist(lapply(fits_list, function(x) x$value)))
  exposure  <- unique(fits_list[[1]]$metadata) # Restore from your old logic
  
  l_results <- list()
  for(val in exp_vals) {
    i_res <- data.frame(matrix(NA, nrow = length(features), ncol = 11 + length(lvl_batch)))
    colnames(i_res) <- c("feature", "exposure", "coef", "stderr", "pval", "k", 
                         "tau2", "stderr.tau2", "pval.tau2", "I2", "H2", paste0("weight_", lvl_batch))
    i_res$feature <- features; i_res$exposure <- val; rownames(i_res) <- features
    
    # RESTORE: Forest plot PDF initialization
    if(!is.null(control$forest_plot)) {
      pdf(file.path(control$output, paste0(exposure, "_", val, "_", stream, "_", control$forest_plot)),
          width = 6,
          height = 4 + ifelse(n_batch > 4, (n_batch - 4) * 0.5, 0))
    }
    
    betas <- vapply(fits_list, function(x) x$coef[x$value == val][match(features, x$feature[x$value == val])], numeric(length(features)))
    pvals <- vapply(fits_list, function(x) x$pval_individual[x$value == val][match(features, x$feature[x$value == val])], numeric(length(features)))
    if(stream == "abundance") {
      # invert p-values from maaslin3 to obtain "true" sds if assuming betas follow normal distribution
      # pvalues of 1 causes infinite sds
      # replace those; this does not matter for meta-analyses as such studies will have
      # near-zero weights regardless
      pval_max_default <- 0.999999
      max_pval <- max(setdiff(pvals, 1), na.rm = TRUE)
      if(max_pval < pval_max_default) max_pval <- pval_max_default
      pvals[pvals %in% 1] <- max_pval
      sds   <- vapply(seq(1, ncol(betas)), 
                      function(i) abs(betas[, i]) / qnorm(1 - pvals[, i] / 2), numeric(length(features)))
    } else {
      sds   <- vapply(fits_list, function(x) x$stderr[x$value == val][match(features, x$feature[x$value == val])], numeric(length(features)))
    }
    
    rownames(betas) <- rownames(sds) <- rownames(pvals) <- features
    
    ind_f <- !is.na(betas) & !is.na(sds) & (sds != 0)
    count_f <- apply(ind_f, 1, sum)
    
    for(f in features) {
      if(count_f[f] >= 2) {
        m_fit <- try(metafor::rma.uni(yi = betas[f, ind_f[f, ]], 
                                      sei = sds[f, ind_f[f, ]],
                                      slab = lvl_batch[ind_f[f, ]],
                                      method = control$rma_method,
                                      control = list(threshold = control$rma_conv, 
                                                     maxiter = control$rma_maxit)), silent = TRUE)
        
        if(!inherits(m_fit, "try-error")) {
          wts <- metafor::weights.rma.uni(m_fit)
          i_res[f, c("coef", "stderr", "pval", "k", "tau2", "stderr.tau2", 
                     "pval.tau2", "I2", "H2", paste0("weight_", names(wts)))] <- 
            c(m_fit$beta, m_fit$se, m_fit$pval, m_fit$k, m_fit$tau2, 
              m_fit$se.tau2, m_fit$QEp, m_fit$I2, m_fit$H2, wts)
          
          # Forest plot generation for significant hits
          if(m_fit$pval < 0.05 & !is.null(control$forest_plot)) {
            metafor::forest(
              m_fit,
              xlab = shorten_name(f, cutoff = 15),
              slab = shorten_name(lvl_batch[ind_f[f, ]], cutoff = 5))
          }
        }
      } else if(count_f[f] == 1) {
        idx_b <- ind_f[f, ]
        i_res[f, c("coef", "stderr", "pval", "k", paste0("weight_", lvl_batch[idx_b]))] <- 
          c(betas[f, idx_b], sds[f, idx_b], pvals[f, idx_b], 1, 100)
      }
    }
    
    if(!is.null(control$forest_plot)) dev.off()
    
    i_res$qval.fdr <- p.adjust(i_res$pval, method = "fdr")
    l_results[[val]] <- i_res
  }
  return(do.call(rbind, l_results))
}

#' Temporary renaming of samples/features for formula safety
#' @keywords internal
rename_maaslin <- function(old_names, prefix) {
  if (is.null(old_names) || length(old_names) == 0) return(NULL)
  new_names <- paste0(prefix, seq_along(old_names))
  names(new_names) <- old_names
  return(new_names)
}

#' Generate dummy grid for Maaslin results table
#' @keywords internal
create_table_maaslin <- function(features, exposure, lvl_exposure) {
  vals <- if (is.null(lvl_exposure)) exposure else lvl_exposure[-1]
  names(features) <- NULL
  expand.grid(feature = features, metadata = exposure, value = vals, stringsAsFactors = FALSE)
}


#' Capture values, warnings, and errors in a list
#' @keywords internal
catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e
    NULL
  }
  val <- withCallingHandlers(tryCatch(expr, error = eHandler), warning = wHandler)
  list(value = val, warnings = myWarnings, error = myError)
}