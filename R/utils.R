#' Wrapper for different ways of data normalization
#'
#' @param physeq phyloseq object to import
#' @param method method for normalization
#' @param rarefy_size the library size used for rarefying. Only meaningful if method = 'Rarefy'
#' @return phyloseq object with normalized abundance table
#' @export
#' @import phyloseq magrittr edgeR
#' @importFrom  metagenomeSeq newMRexperiment
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
    p <- metagenomeSeq::cumNormStatFast(mat_otu %>% metagenomeSeq::newMRexperiment)
    otu_table(physeq_norm)@.Data <- metagenomeSeq::cumNorm(mat_otu %>% metagenomeSeq::newMRexperiment, p = p) %>%
      metagenomeSeq::MRcounts(norm = TRUE, log = FALSE)
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
#' @return fitted PERMANOVA result
#' @export
#'
#' @import phyloseq magrittr vegan
PERMANOVA <- function(physeq, dist = 'bray', model, nperm = 99) {
  if(is(dist, 'character')) {
    cat('Computing distance... \n')
    dist <- phyloseq::distance(physeq, method = dist)
  }
  else if(!is(dist, 'dist')) stop('dist must be either a distance object or dissimilarity measure!')

  df_meta <- sample_data(physeq)
  class(df_meta) <- 'data.frame'

  cat('Running PERMANOVA... \n')
  fit_PERMANOVA <- model %>% as.character %>% c('dist', .) %>% paste(collapse = '') %>%
    as.formula %>% vegan::adonis(data = df_meta, permutations = nperm)

  return(fit_PERMANOVA)
}

#' Batch effect adjustment, only for simulation
#'
#' @param physeq phyloseq object with otu table and metadata containing variables
#' to adjust for in normalizing
#' @param batch batch variable
#' @param adj_model formula for covariates to adjust for
#'
#' @return a list with normalized phyloseq object and simulation details
#' @export
#' @import phyloseq magrittr
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

#' Batch effect adjustment using original ComBat, only for simulation
#'
#' @param physeq phyloseq object with otu table and metadata containing variables
#' to adjust for in normalizing
#' @param batch batch variable
#' @param adj_model formula for covariates to adjust for
#'
#' @return a list with normalized phyloseq object and simulation details
#' @export
#' @import phyloseq magrittr
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
    if (mean.only) {
      # need to figure out why delta is not estimated in this case
      # and why effective sample size is treated as 1
      gamma.star <- postmean(gamma.hat[i, ], 0,
               1, 1, t2[i])
      delta.star <- rep(1, nrow(s.data))
    }
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
postmean_mod <- function (g.hat, g.bar, n, d.star, t2)
{
  (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
}

#' Modified posterior var function from sva, used for modified ComBat
postvar_mod <- function (sum2, n, a, b)
{
  sum2 / n
}

#' Modified iterated solving function from sva, used for modified ComBat
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
#' @import sva
ComBat_forSim <- function (dat, batch, mod = NULL,
                           mean.only = FALSE, BPPARAM = bpparam("SerialParam"))
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

  grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,
                                                   ])

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

  gamma.star <- delta.star <- matrix(NA, nrow = n.batch, ncol = nrow(s.data))

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

  bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]

  return(list(bayesdata = bayesdata, gamma.hat = gamma.hat, gamma.star = gamma.star,
              delta.hat = delta.hat, delta.star = delta.star))
}


#' #' Prior for alpha from sva, used for original ComBat
#' #'
#' #' @param gamma.hat
#' #'
#' #' @return
#' #'
#' #' @examples
#' aprior <- function (gamma.hat)
#' {
#'   m <- mean(gamma.hat)
#'   s2 <- var(gamma.hat)
#'   (2 * s2 + m^2)/s2
#' }
#'
#' #' Prior for beta from sva, used for original ComBat
#' #'
#' #' @param gamma.hat
#' #'
#' #' @return
#' #'
#' #' @examples
#' bprior <- function (gamma.hat)
#' {
#'   m <- mean(gamma.hat)
#'   s2 <- var(gamma.hat)
#'   (m * s2 + m^3)/s2
#' }
#'
#' #' Posterior mean function from sva, used for original ComBat
#' #'
#' #' @param g.hat
#' #' @param g.bar
#' #' @param n
#' #' @param d.star
#' #' @param t2
#' #'
#' #' @return
#' #'
#' #' @examples
#' postmean <- function (g.hat, g.bar, n, d.star, t2)
#' {
#'   (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
#' }
#'
#' #' Posterior var function from sva, used for original ComBat
#' #'
#' #' @param sum2
#' #' @param n
#' #' @param a
#' #' @param b
#' #'
#' #' @return
#' #'
#' #' @examples
#' postvar <- function (sum2, n, a, b)
#' {
#'   (0.5 * sum2 + b)/(n/2 + a - 1)
#' }
#'
#' #' Iterated solving function from sva, used for original ComBat
#' #'
#' #' @param sdat
#' #' @param g.hat
#' #' @param d.hat
#' #' @param g.bar
#' #' @param t2
#' #' @param a
#' #' @param b
#' #' @param conv
#' #'
#' #' @return
#' #'
#' #' @examples
#' it.sol <- function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 1e-04)
#' {
#'   n <- rowSums(!is.na(sdat))
#'   g.old <- g.hat
#'   d.old <- d.hat
#'   change <- 1
#'   count <- 0
#'   while (change > conv) {
#'     g.new <- postmean(g.hat, g.bar, n, d.old, t2)
#'     sum2 <- rowSums((sdat - g.new %*% t(rep(1, ncol(sdat))))^2,
#'                     na.rm = TRUE)
#'     d.new <- postvar(sum2, n, a, b)
#'     change <- max(abs(g.new - g.old)/g.old, abs(d.new -
#'                                                   d.old)/d.old)
#'     g.old <- g.new
#'     d.old <- d.new
#'     count <- count + 1
#'   }
#'   adjust <- rbind(g.new, d.new)
#'   rownames(adjust) <- c("g.star", "d.star")
#'   adjust
#' }
#'
#' #' Row variance function from sva
#' #'
#' #' @param x
#' #' @param na.rm
#' #'
#' #' @return
#' #'
#' #' @examples
#' rowVars <- function (x,na.rm = TRUE)
#' {
#'   sqr = function(x) x * x
#'   n = rowSums(!is.na(x))
#'   n[n <= 1] = NA
#'   return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
#' }
