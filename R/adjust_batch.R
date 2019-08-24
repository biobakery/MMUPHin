#' Zero-inflated empirical Bayes adjustment of batch effect in compositional
#' feature abundance data
#'
#' @param feature_abd feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param batch name of the batch variable. This variable in data should be a
#' factor variable and will be converted to so with a warning if otherwise.
#' @param covariates name(s) of covariates to adjust for
#' in the batch correction model.
#' @param data data frame of metadata, columns must include batch and covariates
#' (if specified).
#' @param zero_inflation logical. Whether or not a zero-inflated model should be
#' run. Default to TRUE (zero-inflated model). If set to FALSE then ComBat
#' (with parametric adjustment) will be performed.
#' @param pseudo_count numeric. Pseudo count to add feature_abd before
#' log-transformation. Automatically set to half of minimal non-zero value if
#' not specified.
#' @param diagnostics logical. Whether or not diagnostic plots are generated.
#' @param verbose logical. Whether or not verbose modelling information is
#' printed.
#'
#' @return feature-by-sample matrix of adjusted feature abundance.
#' @export
#'
adjust_batch <- function(feature_abd,
                         batch,
                         covariates = NULL,
                         data,
                         zero_inflation = TRUE,
                         pseudo_count = NULL,
                         diagnostics = FALSE,
                         verbose = TRUE) {
  # Check data formats
  # Check feature abundance table
  feature_abd <- as.matrix(feature_abd)
  type_feature_abd <- check_feature_abd(feature_abd)
  if(verbose)
    message("feature_abd is ", type_feature_abd)
  # Check metadata data frame
  if(length(batch) > 1)
    stop("Only one batch variable is supported!")
  data <- as.data.frame(data)
  samples <- check_samples(feature_abd, data)
  # Check batch and covariates are included in metadata data frame
  df_batch <- check_metadata(data, batch)
  df_covariates <- check_metadata(data, covariates)
  # Check batch variable
  df_batch[[batch]] <- check_batch(df_batch[[batch]])


  # Construct batch and covariate model matrices. Check for confounding
  batchmod <- construct_design(df_batch)
  mod <- construct_design(df_covariates)
  if(qr(mod)$rank < ncol(mod))
    stop("Covariates are confounded!")
  design <- cbind(batchmod, mod)
  if(qr(design)$rank < ncol(design))
    stop("Covariates and batch are confounded!")
  if(verbose) {
    message("Found ", nlevels(df_batch[[batch]]), " batches")
    message("Adjusting for ", ncol(mod),
            " covariate(s) or covariate(s) level(s)")
  }

  # Transform data for ComBat fit
  read.depth <- apply(feature_abd, 2, sum)
  if(all(read.depth <= 1)) {
    warning("Feature table appears to be on the relative abundance scale!")
    if(pseudo_count == 0.5) {
      warning("pseudo_count was set to 0.5 which is only appropriate to count data,",
              " setting pseudo_count to be min(feature_abd)/2.")
      pseudo_count <- min(setdiff(feature_abd, 0)) / 2
    }
    # In this case don't renormalize the feature
    log.data <- transformFeatures(feature_abd,
                                  transform = "LOG",
                                  pseudo_count = pseudo_count)
  }
  log.data <- transformFeatures(normalizeFeatures(feature_abd,
                                                  normalization = "TSS",
                                                  pseudo_count = pseudo_count),
                                transform = "LOG")

  # Identify data to adjust for
  ind.data <- matrix(TRUE, nrow(feature_abd), ncol(feature_abd)) # which non-missing feature table values are zero
  ind.data[, !ind.sample] <- FALSE
  ind.gamma <- matrix(TRUE, nrow(feature_abd), n.batch) # which feature x batch pairs are adjustable
  ind.mod <- rep(TRUE, ncol(design) - n.batch) # covariates are always adjusted for
  if(zero_inflation) {
    ind.data[feature_abd == 0] <- FALSE
    for(i.feature in 1:nrow(feature_abd)) {
      i.design <- design[ind.data[i.feature, ], , drop = FALSE]
      i.check.batch <- apply(i.design[, 1:n.batch, drop = FALSE] == 1, 2, any)
      i.design <- i.design[, c(i.check.batch, ind.mod), drop = FALSE]
      if(sum(i.check.batch) > 1 && # should have at least two batches to adjust for
         qr(i.design)$rank == ncol(i.design) && # design matrix should be full rank &&
         nrow(i.design) > ncol(i.design) # design matrix cannot give exact fit in linear regression
      ) {
        ind.gamma[i.feature, ] <- i.check.batch
      } else ind.gamma[i.feature, ] <- FALSE
    }
  }
  # Batch has to have more than one adjustable feature to make EB estimates
  ind.gamma[, apply(ind.gamma, 2, sum) < 2] <- FALSE # FIXME?
  ind.feature <- apply(ind.gamma, 1, any)
  if(verbose) message("(After filtering) adjusting for ", sum(ind.feature), " feature(s).")

  # Standardize data across features
  if(verbose) message("Standardizing Data across features.")
  s.data <- log.data
  l.standFit <- list()
  for(i.feature in 1:nrow(s.data)) {
    if(ind.feature[i.feature]) {
      i.design <- design[ind.data[i.feature, ],
                         c(ind.gamma[i.feature, ], ind.mod),
                         drop = FALSE] # FIXME
      # For debugging, this shouldn't happen
      if(nrow(i.design) <= 1 | ncol(i.design) <= 1) stop("Something wrong happened!" ) # FIXME
      standFit <- standardize_feature(y = s.data[i.feature, ind.data[i.feature, ]],
                                      i.design = i.design,
                                      n.batch = sum(ind.gamma[i.feature, ]))
      s.data[i.feature, ind.data[i.feature, ]] <- standFit$y.stand
      l.standFit[[i.feature]] <- standFit
    } else l.standFit[[i.feature]] <- NULL
  }

  # Estimate per-batch mean/variance
  if(verbose) message("Estimating batch difference parameters and EB priors.")
  gamma.hat <-
    delta.hat <-
    gamma.star <-
    delta.star <-
    matrix(NA, nrow = nrow(s.data), ncol = n.batch)
  for(i.feature in 1:nrow(s.data)) {
    if(ind.feature[i.feature]) {
      i.s.data.batch <- s.data[i.feature, ind.data[i.feature, ]] *
        batchmod[ind.data[i.feature, ], ind.gamma[i.feature, ],
                 drop = FALSE] # FIXME
      # For debugging, this shouldn't happen
      if(nrow(batchmod[ind.data[i.feature, ], ind.gamma[i.feature, ],
                       drop = FALSE]) <= 1 |
         ncol(batchmod[ind.data[i.feature, ], ind.gamma[i.feature, ],
                       drop = FALSE]) <= 1) stop("Something wrong happened!" ) # FIXME
      i.gamma <- apply(i.s.data.batch, 2, mean)
      i.delta <- apply(i.s.data.batch, 2, sd)
      i.delta[is.na(i.delta)] <- 1
      i.delta[i.delta == 0] <- 1
      gamma.hat[i.feature, ind.gamma[i.feature, ]] <- i.gamma
      delta.hat[i.feature, ind.gamma[i.feature, ]] <- i.delta
    }
  }
  t2 <- apply(gamma.hat, 2, var, na.rm = TRUE)
  gamma.bar <- apply(gamma.hat, 2, mean, na.rm = TRUE)
  a.prior <- apply(delta.hat, 2, aprior, na.rm = TRUE)
  b.prior <- apply(delta.hat, 2, bprior, na.rm = TRUE)

  # For debugging, this shouldn't happen
  if(any(apply(!is.na(gamma.hat), 2, sum) < 2) |
     any(apply(!is.na(delta.hat), 2, sum) < 2))
    stop("One batch has only one feature with valid parameter estimate!") # FIXME

  # Shrink per-batch mean/variance
  if(verbose) message("Performing shrinkage adjustments on batch difference parameters.")
  results <- lapply(1:n.batch, function(i.batch) {
    i.s.data <- s.data
    i.s.data[!ind.data] <- NA # set all zeros to NA
    i.s.data[, !as.logical(batchmod[, i.batch])] <- NA # set all other batches to NA
    i.s.data[!ind.gamma[, i.batch], ] <- NA # set features not adjustable to NA
    temp <- it.sol(sdat = i.s.data,
                   g.hat = gamma.hat[, i.batch],
                   d.hat = delta.hat[, i.batch],
                   g.bar = gamma.bar[i.batch],
                   t2 = t2[i.batch],
                   a = a.prior[i.batch],
                   b = b.prior[i.batch])
    gamma.star <- temp[1, ]
    delta.star <- temp[2, ]
    list(gamma.star=gamma.star, delta.star=delta.star)
  })
  for (i.batch in 1:n.batch) {
    gamma.star[, i.batch] <- results[[i.batch]]$gamma.star
    delta.star[, i.batch] <- results[[i.batch]]$delta.star
  }

  # Adjust the data
  if(verbose) message("Performing batch corrections.")
  adj.data <- s.data
  for (i.batch in 1:n.batch){
    i.ind.feature <- !is.na(gamma.star[, i.batch]) & !is.na(delta.star[, i.batch])
    # For debugging, this shouldn't happen
    if(!all(i.ind.feature == ind.gamma[, i.batch])) stop("Something was wrong!") # FIXME
    for(i.feature in 1:nrow(adj.data)) {
      if(i.ind.feature[i.feature]) {
        i.ind.sample <- ind.data[i.feature, ] & as.logical(batchmod[, i.batch])
        adj.data[i.feature, i.ind.sample] <-
          (adj.data[i.feature, i.ind.sample] -
             gamma.star[i.feature, i.batch]) /
          sqrt(delta.star[i.feature, i.batch])
      }
    }
  }

  for(i.feature in 1:nrow(adj.data)) {
    if(ind.feature[i.feature]) {
      standFit <- l.standFit[[i.feature]]
      adj.data[i.feature, ind.data[i.feature, ]] <-
        adj.data[i.feature, ind.data[i.feature, ]] *
        sqrt(standFit$var.pooled) +
        standFit$stand.mean
    }
  }

  # For debugging only, this shouldn't happen
  if(any(is.na(adj.data)))
    stop("There are missing values in the adjusted data!") # FIXME

  adj.data <- 2^adj.data
  adj.data[feature_abd == 0] <- 0
  adj.data <- normalizeFeatures(adj.data, normalization = "TSS")
  feature_abd.adj <- t(t(adj.data) * read.depth)
  dimnames(feature_abd.adj) <- dimnames(feature_abd)

  # If required, generate diagnostic plots
  if(diagnostics)
    diagnostics.adjust.batch(feature_abd = feature_abd,
                             feature_abd.adj = feature_abd.adj,
                             batch = batch,
                             gamma.hat = gamma.hat,
                             gamma.star = gamma.star)

  return(feature_abd.adj)
}

