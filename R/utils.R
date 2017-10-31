#' Wrapper for different ways of data normalization
#'
#' @param physeq phyloseq object to import
#' @param method method for normalization
#' @param rarefy_size the library size used for rarefying. Only meaningful if method = 'Rarefy
#' @return
#' @export
#'
#' @import phyloseq magrittr edgeR metagenomeSeq
#' @importFrom edgeR calcNormFactors
#' @examples
normalize <- function(physeq, method = 'none', rarefy_size = 5000) {
  mat_otu <- otu_table(physeq)@.Data
  physeq_norm <- physeq
  sample_data(physeq_norm)$library_size <- apply(mat_otu, 2, sum)

  if (method == 'TSS') {
    normfactor_TSS <- sample_data(physeq_norm)$library_size
    otu_table(physeq_norm)@.Data <- t(t(mat_otu) / normfactor_TSS)
  }
  else if (method == 'Rarefy') {
    physeq_norm <- rarefy_even_depth(physeq_norm, sample.size = rarefy_size,
                                     rngseed = 711, replace = T,
                                     trimOTUs = F)
  }
  else if (method == 'TMM') {
    normfactor_TMM <- edgeR::calcNormFactors(mat_otu,
                                             method = 'TMM')
    otu_table(physeq_norm)@.Data <- t(t(mat_otu) / normfactor_TMM)
  }
  else if (method == 'CSS') {
    p <- cumNormStatFast(mat_otu %>% newMRexperiment)
    otu_table(physeq_norm)@.Data <- cumNorm(mat_otu %>% newMRexperiment, p = p) %>%
      MRcounts(norm = TRUE, log = FALSE)
  }
  else cat('Set to default (no normalization).\n')
  return(physeq_norm)
}

#' Wrapper for performing PERMANOVA on microbial compositional data
#'
#' @param physeq_before phyloseq object before adjustment
#' @param physeq_after phyloseq object after adjustment
#' @param batch_var batch variable used for correction
#' @param adj_var additional variables to adjust for
#' @param distance distance metric used for PERMANOVA
#' @param nperm number of permutations to use in PERMANOVA
#'
#' @return
#' @export
#'
#' @importFrom vegan adonis
#' @import phyloseq magrittr
#' @examples
PERMANOVA <- function(physeq, dist = 'bray', model, nperm = 99) {
  if(is(dist, 'character')) {
    cat('Computing distance... \n')
    dist <- distance(physeq, method = dist)
  }
  else if(!is(dist, 'dist')) stop('dist must be either a distance object or dissimilarity measure!')

  df_meta <- sample_data(physeq)
  class(df_meta) <- 'data.frame'

  cat('Running PERMANOVA... \n')
  fit_PERMANOVA <- model %>% as.character %>% c('dist', .) %>% paste(collapse = '') %>%
    as.formula %>% adonis(data = df_meta, permutations = nperm)

  return(fit_PERMANOVA)
}

#' Batch effect adjustment, only for simulation
#'
#' @param physeq phyloseq object with otu table and metadata containing variables
#' to adjust for in normalizing
#' @param batch batch variable
#' @param adj_model formula for covariates to adjust for
#'
#' @return
#' @export
#' @import phyloseq magrittr
#' @examples
batch_adj_forSim <- function(physeq, batch, adj_model=NULL) {
  mat_otu <- otu_table(physeq)@.Data
  batch <- as.factor(batch)

  # model for adjustment
  log_library_size <- apply(mat_otu, 2, sum) %>% log
  df_meta <- sample_data(physeq)
  class(df_meta) <- 'data.frame'
  df_model <- data.frame(df_meta, log_library_size = log_library_size)
  mod <- matrix(log_library_size, ncol = 1)
  if (missing(adj_model)) {
    # cat("No adjustment covariates provided, adjusting for library size only...")
  }
  else mod <- cbind(mod, model.matrix(adj_model, df_model))

  # matrix for ComBat to adjust
  no_adj_index <- apply(mat_otu == 0, 1, function(feature_abd) {
    tapply(feature_abd, batch, all) %>% any
  })
  log_dat <- (mat_otu[!no_adj_index, ] + 0.5) %>% log

  # Run ComBat
  cat("Running ComBat...\n")
  sim_ComBat <- ComBat_mod_forSim(log_dat, batch, mod)
  log_dat_adj <- sim_ComBat$bayesdata
  mat_otu_adj <- mat_otu
  mat_otu_adj[!no_adj_index, ] <- exp(log_dat_adj)
  mat_otu_adj[mat_otu == 0] <- 0

  physeq_return <- physeq
  otu_table(physeq_return) <- otu_table(mat_otu_adj, taxa_are_rows = T)
  return(list(physeq = physeq_return, sim_ComBat = sim_ComBat))
}


#' Batch effect adjustment using original ComBat, only for simulation for the paper
#'
#' @param physeq
#' @param batch
#' @param adj_model
#'
#' @return
#' @export
#'
#' @examples
batch_adj_ComBat_forSim <- function(physeq, batch, adj_model=NULL) {
  mat_otu <- otu_table(physeq)@.Data
  batch <- as.factor(batch)

  # model for adjustment
  log_library_size <- apply(mat_otu, 2, sum) %>% log
  df_meta <- sample_data(physeq)
  class(df_meta) <- 'data.frame'
  df_model <- data.frame(df_meta, log_library_size = log_library_size)
  mod <- matrix(log_library_size, ncol = 1)
  if (missing(adj_model)) {
    # cat("No adjustment covariates provided, adjusting for library size only...")
  }
  else mod <- cbind(mod, model.matrix(adj_model, df_model))

  # matrix for ComBat to adjust
  no_adj_index <- apply(mat_otu == 0, 1, function(feature_abd) {
    tapply(feature_abd, batch, all) %>% any
  })
  log_dat <- (mat_otu[!no_adj_index, ] + 0.5) %>% log

  # Run ComBat
  cat("Running ComBat...\n")
  sim_ComBat <- ComBat_forSim(log_dat, batch, mod)
  log_dat_adj <- sim_ComBat$bayesdata
  mat_otu_adj <- mat_otu
  mat_otu_adj[!no_adj_index, ] <- exp(log_dat_adj)
  mat_otu_adj[mat_otu == 0] <- 0

  physeq_return <- physeq
  otu_table(physeq_return) <- otu_table(mat_otu_adj, taxa_are_rows = T)
  return(list(physeq = physeq_return, r))
}


#' Running modified ComBat model, only for simulation
#'
#' @param dat
#' @param batch
#' @param mod
#' @param mean.only
#' @param BPPARAM
#'
#' @return
#'
#' @examples
ComBat_mod_forSim <- function (dat, batch, mod = NULL, mean.only = FALSE,
                               BPPARAM = bpparam("SerialParam"))
{
  if (mean.only == TRUE) {
    message("Using the 'mean only' version of ComBat")
  }
  if (length(dim(batch)) > 1) {
    stop("This version of ComBat only allows one batch variable")
  }
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1 + batch)
  ref <- NULL

  message("Found", nlevels(batch), "batches")
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  }
  n.batches <- sapply(batches, length)
  if (any(n.batches == 1)) {
    mean.only = TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  message("Adjusting for", ncol(design) - ncol(batchmod),
          "covariate(s) or covariate level(s)")
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n.batch + 1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if (ncol(design) > (n.batch + 1)) {
      if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[,
                                                          -c(1:n.batch)]))) {
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
      }
      else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  NAs <- any(is.na(dat))
  if (NAs) {
    message(c("Found", sum(is.na(dat)), "Missing Data Values"),
            sep = " ")
  }
  cat("Standardizing Data across genes\n")
  if (!NAs) {
    B.hat <- solve(crossprod(design), tcrossprod(t(design),
                                                 as.matrix(dat)))
  }
  else {
    B.hat <- apply(dat, 1, Beta.NA, design)
  }

  grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,
                                                   ])

  if (!NAs) {
    var.pooled <- ((dat - t(design %*% B.hat))^2) %*%
      rep(1/n.array, n.array)
  }
  else {
    var.pooled <- rowVars(dat - t(design %*% B.hat),
                          na.rm = TRUE)
  }
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1,
                                                           n.array)))
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs) {
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  }
  else {
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design)
  }
  delta.hat <- NULL
  for (i in batches) {
    if (mean.only == TRUE) {
      delta.hat <- rbind(delta.hat, rep(1, nrow(s.data)))
    }
    else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[, i],
                                            na.rm = TRUE))
    }
  }
  t2 <- rowVars(t(gamma.hat))

  gamma.star <- delta.star <- matrix(NA, nrow = n.batch, ncol = nrow(s.data))
  message("Finding parametric adjustments")
  results <- bplapply(1:n.batch, function(i) {
    # if (mean.only) {
    #   # need to figure out why delta is not estimated in this case
    #   # and why effective sample size is treated as 1
    #   gamma.star <- postmean(gamma.hat[i, ], 0,
    #            1, 1, t2[i])
    #   delta.star <- rep(1, nrow(s.data))
    # }
    {
      temp <- it.sol_mod(s.data[, batches[[i]]],
                         gamma.hat[i, ],
                         delta.hat[i, ],
                         0, t2)
      gamma.star <- temp[1, ]
      delta.star <- temp[2, ]
    }
    list(gamma.star = gamma.star, delta.star = delta.star)
  }, BPPARAM = BPPARAM)
  for (i in 1:n.batch) {
    gamma.star[i, ] <- results[[i]]$gamma.star
    delta.star[i, ] <- results[[i]]$delta.star
  }
  message("Adjusting the Data\n")
  bayesdata <- s.data
  j <- 1
  for (i in batches) {
    bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i,
                                                       ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1,
                                                                                                           n.batches[j])))
    j <- j + 1
  }
  bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1,
                                                        n.array)))) + stand.mean
  return(list(bayesdata = bayesdata, gamma.hat = gamma.hat, gamma.star = gamma.star,
              delta.hat = delta.hat, delta.star = delta.star))
}


#' Modified posterior mean function from sva, used for modified ComBat
#'
#' @param g.hat
#' @param g.bar
#' @param n
#' @param d.star
#' @param t2
#'
#' @return
#'
#' @examples
postmean_mod <- function (g.hat, g.bar, n, d.star, t2)
{
  (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
}

#' Modified posterior var function from sva, used for modified ComBat
#'
#' @param sum2
#' @param n
#' @param a
#' @param b
#'
#' @return
#'
#' @examples
postvar_mod <- function (sum2, n, a, b)
{
  sum2 / n
}

#' Modified iterated solving function from sva, used for modified ComBat
#'
#' @param sdat
#' @param g.hat
#' @param d.hat
#' @param g.bar
#' @param t2
#' @param conv
#'
#' @return
#'
#' @examples
it.sol_mod <- function(sdat, g.hat, d.hat, g.bar, t2, conv = 1e-04) {
  n <- rowSums(!is.na(sdat))
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while (change > conv) {
    g.new <- postmean_mod(g.hat, g.bar, n, d.old, t2)
    sum2 <- rowSums((sdat - g.new %*% t(rep(1, ncol(sdat))))^2,
                    na.rm = TRUE)
    d.new <- postvar_mod(sum2, n)
    change <- max(abs(g.new - g.old)/g.old, abs(d.new - d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count + 1
  }
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star", "d.star")
  adjust
}


#' Running original ComBat model, only for simulation
#'
#' @param dat
#' @param batch
#' @param mod
#' @param par.prior
#' @param prior.plots
#' @param mean.only
#' @param ref.batch
#' @param BPPARAM
#'
#' @return
#'
#' @examples
ComBat_forSim <- function (dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                           mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam"))
{
  if (mean.only == TRUE) {
    message("Using the 'mean only' version of ComBat")
  }
  if (length(dim(batch)) > 1) {
    stop("This version of ComBat only allows one batch variable")
  }
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1 + batch)
  if (!is.null(ref.batch)) {
    if (!(ref.batch %in% levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    cat("Using batch =", ref.batch, "as a reference batch (this batch won't change)\n")
    ref <- which(levels(as.factor(batch)) == ref.batch)
    batchmod[, ref] <- 1
  }
  else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  }
  n.batches <- sapply(batches, length)
  if (any(n.batches == 1)) {
    mean.only = TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  design <- cbind(batchmod, mod)
  check <- apply(design, 2, function(x) all(x == 1))
  if (!is.null(ref)) {
    check[ref] <- FALSE
  }
  design <- as.matrix(design[, !check])
  message("Adjusting for", ncol(design) - ncol(batchmod),
          "covariate(s) or covariate level(s)")
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n.batch + 1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if (ncol(design) > (n.batch + 1)) {
      if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[,
                                                          -c(1:n.batch)]))) {
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
      }
      else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  NAs <- any(is.na(dat))
  if (NAs) {
    message(c("Found", sum(is.na(dat)), "Missing Data Values"),
            sep = " ")
  }
  cat("Standardizing Data across genes\n")
  if (!NAs) {
    B.hat <- solve(crossprod(design), tcrossprod(t(design),
                                                 as.matrix(dat)))
  }
  else {
    B.hat <- apply(dat, 1, Beta.NA, design)
  }
  if (!is.null(ref.batch)) {
    grand.mean <- t(B.hat[ref, ])
  }
  else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,
                                                     ])
  }
  if (!NAs) {
    if (!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat - t(design[batches[[ref]],
                                         ] %*% B.hat))^2) %*% rep(1/n.batches[ref], n.batches[ref])
    }
    else {
      var.pooled <- ((dat - t(design %*% B.hat))^2) %*%
        rep(1/n.array, n.array)
    }
  }
  else {
    if (!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat - t(design[batches[[ref]],
                                               ] %*% B.hat), na.rm = TRUE)
    }
    else {
      var.pooled <- rowVars(dat - t(design %*% B.hat),
                            na.rm = TRUE)
    }
  }
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1,
                                                           n.array)))
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs) {
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  }
  else {
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design)
  }
  delta.hat <- NULL
  for (i in batches) {
    if (mean.only == TRUE) {
      delta.hat <- rbind(delta.hat, rep(1, nrow(s.data)))
    }
    else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[, i],
                                            na.rm = TRUE))
    }
  }
  gamma.bar <- rowMeans(gamma.hat)
  t2 <- rowVars(gamma.hat)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  if (prior.plots && par.prior) {
    par(mfrow = c(2, 2))
    tmp <- density(gamma.hat[1, ])
    plot(tmp, type = "l", main = expression(paste("Density Plot of First Batch ",
                                                  hat(gamma))))
    xx <- seq(min(tmp$x), max(tmp$x), length = 100)
    lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
    qqnorm(gamma.hat[1, ], main = expression(paste("Normal Q-Q Plot of First Batch ",
                                                   hat(gamma))))
    qqline(gamma.hat[1, ], col = 2)
    tmp <- density(delta.hat[1, ])
    xx <- seq(min(tmp$x), max(tmp$x), length = 100)
    tmp1 <- list(x = xx, y = dinvgamma(xx, a.prior[1], b.prior[1]))
    plot(tmp, typ = "l", ylim = c(0, max(tmp$y, tmp1$y)),
         main = expression(paste("Density Plot of First Batch ",
                                 hat(delta))))
    lines(tmp1, col = 2)
    invgam <- 1/qgamma(1 - ppoints(ncol(delta.hat)), a.prior[1],
                       b.prior[1])
    qqplot(invgam, delta.hat[1, ], main = expression(paste("Inverse Gamma Q-Q Plot of First Batch ",
                                                           hat(delta))), ylab = "Sample Quantiles", xlab = "Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
  }
  gamma.star <- delta.star <- matrix(NA, nrow = n.batch, ncol = nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        gamma.star <- postmean(gamma.hat[i, ], gamma.bar[i],
                               1, 1, t2[i])
        delta.star <- rep(1, nrow(s.data))
      }
      else {
        temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i,
                                                         ], delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                       b.prior[i])
        gamma.star <- temp[1, ]
        delta.star <- temp[2, ]
      }
      list(gamma.star = gamma.star, delta.star = delta.star)
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i, ] <- results[[i]]$gamma.star
      delta.star[i, ] <- results[[i]]$delta.star
    }
  }
  else {
    message("Finding nonparametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i, ], delta.hat[i, ])
      list(gamma.star = temp[1, ], delta.star = temp[2,
                                                     ])
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i, ] <- results[[i]]$gamma.star
      delta.star[i, ] <- results[[i]]$delta.star
    }
  }
  if (!is.null(ref.batch)) {
    gamma.star[ref, ] <- 0
    delta.star[ref, ] <- 1
  }
  message("Adjusting the Data\n")
  bayesdata <- s.data
  j <- 1
  for (i in batches) {
    bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i,
                                                       ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1,
                                                                                                           n.batches[j])))
    j <- j + 1
  }
  bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1,
                                                        n.array)))) + stand.mean
  if (!is.null(ref.batch)) {
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  return(list(bayesdata = bayesdata, gamma.hat = gamma.hat, gamma.star = gamma.star,
              delta.hat = delta.hat, delta.star = delta.star))
}


#' Prior for alpha from sva, used for original ComBat
#'
#' @param gamma.hat
#'
#' @return
#'
#' @examples
aprior <- function (gamma.hat)
{
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (2 * s2 + m^2)/s2
}

#' Prior for beta from sva, used for original ComBat
#'
#' @param gamma.hat
#'
#' @return
#'
#' @examples
bprior <- function (gamma.hat)
{
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (m * s2 + m^3)/s2
}

#' Posterior mean function from sva, used for original ComBat
#'
#' @param g.hat
#' @param g.bar
#' @param n
#' @param d.star
#' @param t2
#'
#' @return
#'
#' @examples
postmean <- function (g.hat, g.bar, n, d.star, t2)
{
  (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
}

#' Posterior var function from sva, used for original ComBat
#'
#' @param sum2
#' @param n
#' @param a
#' @param b
#'
#' @return
#'
#' @examples
postvar <- function (sum2, n, a, b)
{
  (0.5 * sum2 + b)/(n/2 + a - 1)
}

#' Iterated solving function from sva, used for original ComBat
#'
#' @param sdat
#' @param g.hat
#' @param d.hat
#' @param g.bar
#' @param t2
#' @param a
#' @param b
#' @param conv
#'
#' @return
#'
#' @examples
it.sol <- function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 1e-04)
{
  n <- rowSums(!is.na(sdat))
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while (change > conv) {
    g.new <- postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- rowSums((sdat - g.new %*% t(rep(1, ncol(sdat))))^2,
                    na.rm = TRUE)
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new - g.old)/g.old, abs(d.new -
                                                  d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count + 1
  }
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star", "d.star")
  adjust
}

#' Row variance function from sva
#'
#' @param x
#' @param na.rm
#'
#' @return
#'
#' @examples
rowVars <- function (x,na.rm = TRUE)
{
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}



#' Modified linear model from metagenomeSeq
#'
#' @param obj
#' @param mod
#' @param coef
#' @param szero
#' @param spos
#'
#' @return
#' @export
#'
#' @examples
fitFeatureModel2 <- function (obj, mod, coef = 2, szero = FALSE, spos = TRUE)
{
  stopifnot(is(obj, "MRexperiment"))
  if (any(is.na(normFactors(obj))))
    stop("At least one NA normalization factors")
  if (any(is.na(libSize(obj))))
    stop("Calculate the library size first!")
  if (any(is.na(normFactors(obj)))) {
    stop("Calculate the normalization factors first!")
  }
  mmCount = cbind(mod, log(normFactors(obj)/median(normFactors(obj))))
  colnames(mmCount)[ncol(mmCount)] = "scalingFactor"
  if (ncol(mmCount) > 3) {
    stop("Can't analyze currently.")
  }
  i = permuttedFits = NULL
  fitzeroln = fitZeroLogNormal2(obj, mmCount, coef = coef,
                                szero = szero, spos = spos)
  zscore = fitzeroln$logFC/fitzeroln$se

  pvals = 2 * (1 - pnorm(abs(zscore)))
  res = list(call = match.call(), fitZeroLogNormal = fitzeroln,
             design = mmCount, taxa = rownames(obj), counts = MRcounts(obj),
             pvalues = pvals, permuttedFits = permuttedFits)
}




#' Modified log normal model from metagenomeSeq
#'
#' @param obj
#' @param mod
#' @param coef
#' @param szero
#' @param spos
#'
#' @return
#'
#' @examples
fitZeroLogNormal2 <- function (obj, mod, coef = 2, szero = TRUE, spos = TRUE)
{
  positiveMod = mod[, -ncol(mod)]
  zeroMod = mod
  mat <- MRcounts(obj, norm = FALSE, log = FALSE, sl = 1)
  posIndices = mat > 0
  nr = nrow(mat)
  nc = ncol(mat)
  exclude = zeroExclude = tauZero = tauPos = posRidge = zeroRidge = NULL
  results = array(NA, dim = c(nr, 3))
  rownames(results) = rownames(mat)
  colnames(results) = c("logFC", "adjFactor", "se")
  fitln = calcPosComponent(mat, positiveMod, posIndices)
  zeros2 = which(fitln[, "s2"] == 0)
  rs = rowsum(t(1 - (1 - posIndices)), positiveMod[, coef])
  exclude = union(which(rs[1, ] <= 1), which(rs[2, ] <= 1))
  zeroExclude = which(colSums(rs) >= (nc - 3))
  exclude = union(zeros2, exclude)
  if (length(exclude) == 0)
    exclude = NULL
  if (length(zeroExclude) == 0)
    zeroExclude = NULL
  sdensity = density(fitln[, "s2"], na.rm = TRUE)
  smode = sdensity$x[which.max(sdensity$y)]
  if (length(zeros2) > 0)
    fitln[zeros2, "s2"] = smode
  if (spos == TRUE) {
    shrinkPos <- calcShrinkParameters(fitln, coef, smode,
                                      exclude)
    tauPos = shrinkPos$tau
    vpost = shrinkPos$v.post
    fitln[, "s2"] = vpost
    posRidge = sapply(seq(nr), function(i) {
      k = which(posIndices[i, ])
      y = log(mat[i, k])
      x = positiveMod[k, ]
      l = vpost[i]/(nrow(x) * tauPos)
      if (i %in% exclude)
        return(matrix(rep(NA, ncol(positiveMod))))
      ridge = glmnet(y = y, x = x, lambda = l, alpha = 0)
      as.matrix(coefficients(ridge)[colnames(positiveMod),
                                    ])
    })
    posFittedCoefficients = t(posRidge)
    rownames(posFittedCoefficients) = rownames(mat)
    fitln[rownames(posFittedCoefficients), 1:ncol(positiveMod)] = posFittedCoefficients
  }
  fitzero = calcZeroComponent(mat, zeroMod, posIndices)
  sdensity = density(fitzero[, "s2"], na.rm = TRUE)
  smode = sdensity$x[which.max(sdensity$y)]
  if (length(exclude) > 0)
    fitzero[exclude, "s2"] = smode
  if (szero == TRUE) {
    shrinkZero <- calcShrinkParameters(fitzero, coef, smode,
                                       exclude)
    tauZero = shrinkZero$tau
    vpostZero = shrinkZero$v.post
    fitzero[, "s2"] = vpostZero
    zeroRidge = sapply(1:nr, function(i) {
      y = posIndices[i, ]
      l = 1/(nc * tauZero)
      if (i %in% c(zeroExclude, exclude))
        return(matrix(rep(NA, ncol(zeroMod))))
      ridge = glmnet(y = y, x = zeroMod, lambda = l, family = "binomial",
                     alpha = 0, penalty.factor = c(rep(1, (ncol(zeroMod) -
                                                             1)), 0))
      as.matrix(coefficients(ridge))[colnames(zeroMod),
                                     ]
    })
    zeroFittedCoefficients = t(zeroRidge)
    rownames(zeroFittedCoefficients) = rownames(mat)
    fitzero[rownames(zeroFittedCoefficients), 1:ncol(zeroMod)] = zeroFittedCoefficients
  }
  se = calcStandardError(zeroMod, fitln, fitzero, coef = coef,
                         exclude = union(exclude, zeroExclude))
  se[zeroExclude] = sqrt(fitln[zeroExclude, "s2"])
  adjFactor = calcZeroAdjustment(fitln, fitzero, zeroMod,
                                 coef, exclude = exclude)
  adjFactor[zeroExclude] = 0
  logFC <- fitln[, coef] + adjFactor
  list(logFC = logFC, adjFactor = adjFactor, se = se, fitln = fitln,
       fitzero = fitzero, zeroRidge = zeroRidge, posRidge = posRidge,
       tauPos = tauPos, tauZero = tauZero, exclude = exclude,
       zeroExclude = zeroExclude)
}

#' Utility function from metagenomeSeq
#'
#' @param mat
#' @param mod
#' @param weights
#'
#' @return
#'
#' @examples
calcPosComponent <- function (mat, mod, weights)
{
  fitln <- lmFit(log(mat), mod, weights = weights)
  b = coefficients(fitln)
  df = fitln$df
  res = residuals(fitln, log(mat))
  s2 = sapply(seq(nrow(res)), function(i) {
    sum(res[i, which(weights[i, ])]^2, na.rm = TRUE)/df[i]
  })
  fitln <- data.frame(b = b, s2 = s2, df = df)
  rownames(fitln) = rownames(mat)
  fitln
}

#' Utility function from metagenomeSeq
#'
#' @param fit
#' @param coef
#' @param mins2
#' @param exclude
#'
#' @return
#'
#' @examples
calcShrinkParameters <- function (fit, coef, mins2, exclude = NULL)
{
  if (is.null(exclude)) {
    shrunkVar <- limma::squeezeVar(fit[, "s2"], fit[, "df"])
    v.post = shrunkVar$var.post
    tau <- var(fit[, coef], na.rm = TRUE)
  }
  else {
    v.post = rep(mins2, nrow(fit))
    shrunkVar <- limma::squeezeVar(fit[-exclude, "s2"], fit[-exclude,
                                                            "df"])
    v.post[-exclude] <- shrunkVar$var.post
    tau <- var(fit[-exclude, coef], na.rm = TRUE)
  }
  list(tau = tau, v.post = v.post)
}

#' Utility function from metagenomeSeq
#'
#' @param mat
#' @param mod
#' @param weights
#'
#' @return
#'
#' @examples
calcZeroComponent <- function (mat, mod, weights)
{
  fitzero <- sapply(seq(nrow(mat)), function(i) {
    fit <- glm.fit(mod, weights[i, ], family = binomial())
    cf = coefficients(fit)
    df = fit$df.residual
    mc = exp(mod %*% cf)
    s2 = sum((weights[i, ] - t(mc/(1 + mc)))^2)/df
    c(beta = cf, s2 = s2, df = df)
  })
  fitzero <- data.frame(t(fitzero))
  rownames(fitzero) = rownames(mat)
  fitzero
}

#' Utility function from metagenomeSeq
#'
#' @param mod
#' @param fitln
#' @param fitzero
#' @param coef
#' @param exclude
#'
#' @return
#'
#' @examples
calcStandardError <- function (mod, fitln, fitzero, coef = 2, exclude = NULL)
{
  mod0 = mod1 = mod
  mod1[, coef] <- 1
  mod0[, coef] <- 0
  ve = rep(NA, nrow(fitln))
  features = seq(nrow(fitln))
  if (length(exclude) > 0)
    features = features[-exclude]
  fullvar = sapply(features, function(i) {
    beta = fitzero[i, 1:ncol(mod)]
    b = fitln[i, 1:(ncol(mod) - 1)]
    s = as.numeric(fitln[i, "s2"])
    mu0 = as.vector(exp(mod0[, -ncol(mod)] %*% t(b) + 0.5 *
                          s))
    mu1 = as.vector(exp(mod1[, -ncol(mod)] %*% t(b) + 0.5 *
                          s))
    theta <- mod %*% t(beta)
    theta1 <- mod1 %*% t(beta)
    theta0 <- mod0 %*% t(beta)
    p <- t(exp(theta)/(1 + exp(theta)))
    p1 <- t(exp(theta1)/(1 + exp(theta1)))
    p0 <- t(exp(theta0)/(1 + exp(theta0)))
    checkInverse <- function(m) {
      class(try(qr.solve(m), silent = T)) == "matrix"
    }
    Dp2 <- diag(length(p)) * as.vector(p * (1 - p))
    infz = t(mod) %*% Dp2 %*% mod
    Dp <- diag(length(p)) * as.vector(p)
    infln = t(mod[, -ncol(mod)]) %*% Dp %*% mod[, -ncol(mod)]
    if (checkInverse(infz)) {
      invinf_z <- qr.solve(infz)
    }
    else {
      return(NA)
    }
    if (checkInverse(infln)) {
      invinf_ln <- as.numeric(s) * qr.solve(infln)
    }
    else {
      return(NA)
    }
    invInfFull = as.matrix(bdiag(invinf_z, invinf_ln, (2 *
                                                         s^2/sum(p))))
    logRatioBeta0 <- (mean(p1 * (1 - p1) * mu0)/mean(p1 *
                                                       mu0)) - (mean(p0 * (1 - p0) * mu0)/mean(p0 * mu0))
    logRatioBeta1 <- mean(p1 * (1 - p1) * mu0)/mean(p1 *
                                                      mu0)
    logRatioBeta2 <- (mean(mod[, 3] * p1 * (1 - p1) * mu0)/mean(p1 *
                                                                  mu0)) - (mean(mod[, 3] * p0 * (1 - p0) * mu0)/mean(p0 *
                                                                                                                       mu0))
    logRatioFull = t(c(logRatioBeta0, logRatioBeta1, logRatioBeta2,
                       0, 1, 0))
    logRatioVar = logRatioFull %*% invInfFull %*% t(logRatioFull)
    logRatioVar
  })
  if (!is.null(exclude)) {
    if (length(features) > 0) {
      ve[features] = fullvar
    }
  }
  else {
    ve = fullvar
  }
  sqrt(ve)
}

#' Utility function from metagenomeSeq
#'
#' @param fitln
#' @param fitzero
#' @param mod
#' @param coef
#' @param exclude
#'
#' @return
#'
#' @examples
calcZeroAdjustment <- function (fitln, fitzero, mod, coef, exclude = NULL)
{
  b = fitln[, 1:(ncol(mod) - 1)]
  beta = fitzero[, 1:ncol(mod)]
  mod1 <- mod
  mod1[, coef] <- 1
  theta1 <- mod1 %*% t(beta)
  p1 <- exp(theta1)/(1 + exp(theta1))
  p1 <- t(p1)
  if (ncol(b) > 2)
    p1 = p1 * exp(t(mod[, 3:(ncol(mod) - 1)] %*% t(b[, 3:ncol(b)])))
  mean_p1 <- rowMeans(p1)
  mod0 <- mod
  mod0[, coef] <- 0
  theta0 <- mod0 %*% t(beta)
  p0 <- exp(theta0)/(1 + exp(theta0))
  p0 <- t(p0)
  if (ncol(b) > 2)
    p0 = p0 * exp(t(mod[, 3:(ncol(mod) - 1)] %*% t(b[, 3:ncol(b)])))
  mean_p0 <- rowMeans(p0)
  adjFactor <- log(mean_p1/mean_p0)
  if (!is.null(exclude))
    adjFactor[exclude] = NA
  adjFactor
}
