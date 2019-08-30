#' Construct a non-intercept design model matrix given a metadata data frame
#'
#' @param data metadata data frame.
#'
#' @return design matrix.
construct_design <- function(data) {
  # Returns NULL if data is NULL. This happens if the covariate data frame is
  # NULL (when no covariates are provided)
  if(is.null(data)) return(NULL)

  # Construct the matrix using all variables in data but excluding the intercept
  # term (because batch dummy variables encompass the intercept term)
  model.matrix(~ . - 1, data = data)
}

#' Check if a design matrix is full rank
#'
#' @param design design matrix.
#'
#' @return TRUE/FALSE for whether or not the design matrix is full rank.
check_rank <- function(design) {
  # a zero-column matrix is full rank
  if(is.null(design)) return(TRUE)

  qr(design)$rank == ncol(design)
}


#' Create indicator matrices for which feature/batch/samples to adjust. This is
#' relevant for zero_inflation is TRUE and only non-zero values are adjusted.
#'
#' @param feature_abd feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param n_batch number of batches in the data.
#' @param design design matrix.
#' @param zero_inflation zero inflation flag.
#'
#' @return list of indicator matrices needed by fitting in adjust_batch.
construct_ind <- function(feature_abd, n_batch, design, zero_inflation) {
  # which feature table values are zero
  ind_data <- matrix(TRUE, nrow(feature_abd), ncol(feature_abd))
  # which feature x batch pairs are adjustable
  ind_gamma <- matrix(TRUE, nrow(feature_abd), n_batch)
  # covariates are always adjusted for
  ind_mod <- rep(TRUE, ncol(design) - n_batch)
  if(zero_inflation) {
    ind_data[feature_abd == 0] <- FALSE
    for(i_feature in 1:nrow(feature_abd)) {
      # subset design to non-zero samples for i_feature
      i_design <- design[ind_data[i_feature, ], , drop = FALSE]
      # indicate whether each batch has non-zero values for i_feature
      i_check_batch <- apply(i_design[, 1:n_batch, drop = FALSE] == 1, 2, any)
      i_design <- i_design[, c(i_check_batch, ind_mod), drop = FALSE]
      if(
        # should have at least two batches to adjust for
        sum(i_check_batch) > 1 &&
        # design matrix should be full rank
        qr(i_design)$rank == ncol(i_design) &&
        # design matrix cannot give exact fit in linear regression
        nrow(i_design) > ncol(i_design)
      ) {
        ind_gamma[i_feature, ] <- i_check_batch
      } else ind_gamma[i_feature, ] <- FALSE
    }
  }
  # Batch has to have more than one adjustable feature to make EB estimates
  ind_gamma[, apply(ind_gamma, 2, sum) < 2] <- FALSE
  ind_feature <- apply(ind_gamma, 1, any)
  return(list(ind_data = ind_data,
              ind_gamma = ind_gamma,
              ind_mod = ind_mod,
              ind_feature = ind_feature))
}

#' Fit lm and standardize all features
#'
#' @param s_data feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param design design matrix.
#' @param l_ind list of indicator matrices, as returned by construct_ind.
#'
#' @return list of two componet: the standardized feature abundance matrix, and
#' a list of per-feature standardization fits.
fit_stand_feature <- function(s_data, design, l_ind) {
  l_stand_feature <- list()
  for(i_feature in 1:nrow(s_data)) {
    if(l_ind$ind_feature[i_feature]) {
      i_design <- design[l_ind$ind_data[i_feature, ],
                         c(l_ind$ind_gamma[i_feature, ],
                           l_ind$ind_mod),
                         drop = FALSE]
      # For debugging, this shouldn't happen
      if(nrow(i_design) <= 1 | ncol(i_design) <= 1)
        stop("Something wrong happened!" ) # FIXME
      stand_fit <- standardize_feature(
        y = s_data[i_feature, l_ind$ind_data[i_feature, ]],
        i_design = i_design,
        n_batch = sum(l_ind$ind_gamma[i_feature, ])
      )
      s_data[i_feature, l_ind$ind_data[i_feature, ]] <- stand_fit$y_stand
      l_stand_feature[[i_feature]] <- stand_fit
    } else l_stand_feature[[i_feature]] <- NULL
  }
  return(list(s_data = s_data,
              l_stand_feature = l_stand_feature))
}

#' Centralize (by design matrix) and standardize (by pooled variance across all
#' batches) feature abundances for empirical Bayes fit
#'
#' @param y vector of non-zero abundance of a single feature (if zero-inflated
#' is true).
#' @param i_design design matrix for the feature; samples with zeros are taken
#' out (if zero-inflated is true).
#' @param n_batch number of batches in the data.
#'
#' @return a list with component: y_stand for vector of centralized and
#' standardized feature abundance, and stand_mean/varpooled for the location and
#' scale factor (these are used later to back transform the batch-shrinked
#' feature abundance).
standardize_feature <- function(y,
                                i_design,
                                n_batch) {
  beta_hat <- solve(crossprod(i_design),
                    crossprod(i_design, y))
  grand_mean <- mean(i_design[, 1:n_batch] %*%
                       beta_hat[1:n_batch, ])

  var_pooled <- var(y - (i_design %*% beta_hat)[, 1])
  stand_mean <- rep(grand_mean, length(y))
  if(ncol(i_design) > n_batch){
    stand_mean <- stand_mean +
      (i_design[, -(1:n_batch), drop = FALSE] %*%
         beta_hat[-(1:n_batch), ])[, 1]
  }
  y_stand <- (y - stand_mean) / sqrt(var_pooled)
  return(list(y_stand = y_stand,
              stand_mean = stand_mean,
              var_pooled = var_pooled))
}

#' Parametric estimation of per-batch location and scale parameters, and
#' Empirical Bayes estimation of their priors
#'
#' @param s_data feature-by-sample matrix of standardized abundances.
#' @param l_stand_feature list of per-feature standardization fits, as returned
#' by fit_stand_feature.
#' @param batchmod design matrix for batch variables.
#' @param n_batch number of batches in the data.
#' @param l_ind list of indicator matrices, as returned by construct_ind.
#'
#' @return list of parameter estimations.
fit_EB <- function(s_data, l_stand_feature, batchmod, n_batch, l_ind) {

  if(n_batch != ncol(batchmod))
    stop("n_batch does not agree with batchmod!")

  gamma_hat <-
    delta_hat <-
    matrix(NA, nrow = nrow(s_data), ncol = n_batch)

  # estimate per-feature per-batch location and scale parameters
  for(i_feature in 1:nrow(s_data)) {
    if(l_ind$ind_feature[i_feature]) {
      i_s_data_batch <- s_data[i_feature, l_ind$ind_data[i_feature, ]] *
        batchmod[l_ind$ind_data[i_feature, ], l_ind$ind_gamma[i_feature, ],
                 drop = FALSE]
      # For debugging, this shouldn't happen
      if(
        # less than two samples are non-zero to correct for the feature
        nrow(batchmod[l_ind$ind_data[i_feature, ],
                      l_ind$ind_gamma[i_feature, ],
                      drop = FALSE]) <= 1 |
        # less than two batches are eligible to correct for the feature
        ncol(batchmod[l_ind$ind_data[i_feature, ],
                      l_ind$ind_gamma[i_feature, ],
                      drop = FALSE]) <= 1)
        stop("Something wrong happened!" ) ## FIXME
      i_gamma <- apply(i_s_data_batch, 2, mean)
      i_delta <- apply(i_s_data_batch, 2, sd)
      i_delta[is.na(i_delta)] <- 1
      i_delta[i_delta == 0] <- 1
      gamma_hat[i_feature, l_ind$ind_gamma[i_feature, ]] <- i_gamma
      delta_hat[i_feature, l_ind$ind_gamma[i_feature, ]] <- i_delta
    }
  }

  # EM hyper-parameter estimations
  gamma_bar <- apply(gamma_hat, 2, mean, na.rm = TRUE)
  t2 <- apply(gamma_hat, 2, var, na.rm = TRUE)
  a_prior <- apply(delta_hat, 2, aprior, na.rm = TRUE)
  b_prior <- apply(delta_hat, 2, bprior, na.rm = TRUE)

  # For debugging, this shouldn't happen
  # If a batch has only one feature with valid location/scale parameters
  # Will cause problem for hyper-parameter estimation
  ## FIXME
  if(any(apply(!is.na(gamma_hat), 2, sum) < 2) |
     any(apply(!is.na(delta_hat), 2, sum) < 2))
    stop("One batch has only one feature with valid parameter estimate!")

  return(list(gamma_hat = gamma_hat,
              delta_hat = delta_hat,
              gamma_bar = gamma_bar,
              t2 = t2,
              a_prior = a_prior,
              b_prior = b_prior))
}

#' EB prior estimation for scale parameters
#'
#' @param delta_hat frequentist per-batch scale estimations.
#' @param na.rm whether or not missing values should be removed.
#'
#' @return shape hyper parameter
aprior <- function(delta_hat, na.rm = FALSE) {
  m <- mean(delta_hat, na.rm = na.rm)
  s2 <- var(delta_hat, na.rm = na.rm)
  (2*s2 + m^2) / s2
}

#' EB prior estimation for scale parameters
#'
#' @param delta_hat frequentist per-batch location estimations.
#' @param na.rm whether or not missing values should be removed.
#'
#' @return scale hyper parameter
bprior <- function(delta_hat, na.rm = FALSE){
  m <- mean(delta_hat, na.rm = na.rm)
  s2 <- var(delta_hat, na.rm = na.rm)
  (m*s2 + m^3) / s2
}

#' A posteriori shrink per-batch location and scale parameters towards their EB
#' priors
#'
#' @param s_data feature-by-sample matrix of standardized abundances.
#' @param l_params list of parameter fits, as returned by fit_EB.
#' @param batchmod design matrix for batch variables.
#' @param n_batch number of batches in the data.
#' @param l_ind list of indicator matrices, as returned by construct_ind.
#' @param control list of control parameters (passed on to it_sol)
#'
#' @return list of shrinked per-batch location and scale parameters.
fit_shrink <- function(s_data, l_params, batchmod, n_batch, l_ind, control) {

  if(n_batch != ncol(batchmod))
    stop("n_batch does not agree with batchmod!")

  gamma_star <-
    delta_star <-
    matrix(NA, nrow = nrow(s_data), ncol = n_batch)

  results <- lapply(1:n_batch, function(i_batch) {
    i_s_data <- s_data
    # set all zeros to NA
    i_s_data[!l_ind$ind_data] <- NA
    # set all other batches to NA
    i_s_data[, !as.logical(batchmod[, i_batch])] <- NA
    # set features not adjustable to NA
    i_s_data[!l_ind$ind_gamma[, i_batch], ] <- NA
    temp <- it_sol(s_data = i_s_data,
                   g_hat = l_params$gamma_hat[, i_batch],
                   d_hat = l_params$delta_hat[, i_batch],
                   g_bar = l_params$gamma_bar[i_batch],
                   t2 = l_params$t2[i_batch],
                   a = l_params$a_prior[i_batch],
                   b = l_params$b_prior[i_batch],
                   control = control)
    gamma_star <- temp[1, ]
    delta_star <- temp[2, ]
    list(gamma_star=gamma_star, delta_star=delta_star)
  })
  for (i_batch in 1:n_batch) {
    gamma_star[, i_batch] <- results[[i_batch]]$gamma_star
    delta_star[, i_batch] <- results[[i_batch]]$delta_star
  }

  return(list(gamma_star = gamma_star,
              delta_star = delta_star))
}

#' Iteratively solve for one feature's shrinked location and scale parameters
#'
#' @param s_data the feature's standardized abundances.
#' @param g_hat the feature's location parameter frequentist estimations.
#' @param d_hat the feature's scale parameter frequentist estimations.
#' @param g_bar EB estimation of location hyper parameters.
#' @param t2 EB estimation of location hyper parameters.
#' @param a EB estimation of scale hyper parameters.
#' @param b EB estimation of scale hyper parameters.
#' @param control list of control parameters
#'
#' @return matrix of shrinked location and scale parameters.
it_sol  <- function(s_data,
                    g_hat,
                    d_hat,
                    g_bar,
                    t2,
                    a,
                    b,
                    control){
  n <- rowSums(!is.na(s_data))
  g.old <- g_hat
  d.old <- d_hat
  change <- 1
  count <- 0
  while(change>control$conv){
    g.new <- postmean(g_hat, g_bar, n, d.old, t2)
    sum2 <- rowSums((s_data - g.new %*% t(rep(1,ncol(s_data))))^2, na.rm=TRUE)
    sum2[sum2 == 0] <- NA
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new-g.old) / g.old, abs(d.new-d.old) / d.old,
                  na.rm=TRUE)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
    if(count > control$maxit)
      stop("Maximum iteration reached!")
  }
  ## cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g_star","d_star")
  adjust
}

postmean <- function(g_hat,g_bar,n,d_star,t2){
  (t2*n*g_hat + d_star*g_bar) / (t2*n + d_star)
}

postvar <- function(sum2,n,a,b){
  (.5*sum2 + b) / (n/2 + a - 1)
}

#' Perform batch adjustment on standardized feature abundances, based on EB
#' shrinked per-batch location and scale parameters
#'
#' @param s_data feature-by-sample matrix of standardized abundances.
#' @param l_params_shrink list of shrinked parameters, as returned by
#' fit_shrink.
#' @param l_stand_feature list of per-feature standardization fits, as returned
#' by fit_stand_feature.
#' @param batchmod design matrix for batch variables.
#' @param n_batch number of batches in the data.
#' @param l_ind list of indicator matrices, as returned by construct_ind.
#'
#' @return feature-by-sample matrix of batch-adjusted feature abundances.
adjust_EB <- function(s_data, l_params_shrink, l_stand_feature,
                      batchmod, n_batch,
                      l_ind) {
  if(n_batch != ncol(batchmod))
    stop("n_batch does not agree with batchmod!")
  if(n_batch != ncol(l_params_shrink[[1]]))
    stop("n_batch does not agree with l_params_shrink!")

  adj_data <- relocate_scale(s_data, l_params_shrink,
                             batchmod, n_batch,
                             l_ind)
  adj_data <- add_back_covariates(adj_data, l_stand_feature,
                                  l_ind)

  return(adj_data)
}

#' Relocate and scale feature abundances to correct for batch effects, given
#' shrinked per-batch location and scale parameters
#'
#' @param s_data feature-by-sample matrix of standardized abundances.
#' @param l_params_shrink list of shrinked parameters, as returned by
#' fit_shrink.
#' @param batchmod design matrix for batch variables.
#' @param n_batch number of batches in the data.
#' @param l_ind list of indicator matrices, as returned by construct_ind.
#'
#' @return feature-by-sample matrix of batch-adjusted feature abundances
#' (but without covariate effects).
relocate_scale <- function(s_data, l_params_shrink,
                           batchmod, n_batch,
                           l_ind) {
  adj_data <- s_data
  for (i_batch in 1:n_batch) {
    i_ind_feature <-
      !is.na(l_params_shrink$gamma_star[, i_batch]) &
      !is.na(l_params_shrink$delta_star[, i_batch])
    # For debugging, this shouldn't happen
    # Features with valid shrinked parameters in the batch should agree
    # with the ones determined to be eligible for batch estimation in the
    # first place
    ## FIXME
    if(!all(i_ind_feature == l_ind$ind_gamma[, i_batch]))
      stop("Features determined to be eligible for batch estimation do not ",
           "agree with the ones with valid per-batch shrinked parameters!")

    for(i_feature in 1:nrow(adj_data)) {
      if(i_ind_feature[i_feature]) {
        i_ind_sample <-
          l_ind$ind_data[i_feature, ] &
          as.logical(batchmod[, i_batch])
        # relocate and scale
        adj_data[i_feature, i_ind_sample] <-
          (adj_data[i_feature, i_ind_sample] -
             l_params_shrink$gamma_star[i_feature, i_batch]) /
          sqrt(l_params_shrink$delta_star[i_feature, i_batch])
      }
    }
  }

  return(adj_data)
}

#' Add back covariate effects to batch-corrected feature abundance data
#'
#' @param adj_data feature-by-sample matrix of batch-adjusted feature abundances
#' (but without covariate effects), as returned by relocate_scale.
#' @param l_stand_feature list of per-feature standardization fits, as returned
#' by fit_stand_feature.
#' @param l_ind list of indicator matrices, as returned by construct_ind.
#'
#' @return feature-by-sample matrix of batch-adjusted feature abundances
#' with covariate effects retained.
add_back_covariates <- function(adj_data, l_stand_feature,
                                l_ind) {
  for(i_feature in 1:nrow(adj_data)) {
    if(l_ind$ind_feature[i_feature]) {
      i_stand_feature <- l_stand_feature[[i_feature]]
      adj_data[i_feature, l_ind$ind_data[i_feature, ]] <-
        adj_data[i_feature, l_ind$ind_data[i_feature, ]] *
        sqrt(i_stand_feature$var_pooled) +
        i_stand_feature$stand_mean
    }
  }
  return(adj_data)
}

#' Transform batch adjusted feature abundances back to the original scale in
#' feature_abd
#'
#' @param adj_data feature-by-sample matrix of batch-adjusted feature abundances
#' with covariate effects retained.
#' @param feature_abd original feature-by-sample matrix of abundances
#' (proportions or counts).
#' @param type_feature_abd type of feature abundance table (counts or
#' proportions). If counts, the final output will be rounded into counts as
#' well.
#'
#' @return feature-by-sample matrix of batch-adjusted feature abundances,
#' with covariate effects retained and scales consistent with original abundance
#' matrix.
back_transform_abd <- function(adj_data, feature_abd, type_feature_abd) {
  adj_data <- 2^adj_data
  adj_data[feature_abd == 0] <- 0
  adj_data <- normalize_features(adj_data, normalization = "TSS")
  adj_data <- t(t(adj_data) * apply(feature_abd, 2, sum))
  dimnames(adj_data) <- dimnames(feature_abd)

  if(type_feature_abd == "counts")
    adj_data <- round(adj_data)

  return(adj_data)
}
