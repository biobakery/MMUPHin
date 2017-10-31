#' Main function for batch effect adjustment
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
batch_adj <- function(physeq, batch, adj_model=NULL) {
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
  log_dat_adj <- ComBat_mod(log_dat, batch, mod)
  mat_otu_adj <- mat_otu
  mat_otu_adj[!no_adj_index, ] <- exp(log_dat_adj)
  mat_otu_adj[mat_otu == 0] <- 0

  physeq_return <- physeq
  otu_table(physeq_return) <- otu_table(mat_otu_adj, taxa_are_rows = T)
  return(physeq_return)
}

#' Running modified ComBat model
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
ComBat_mod <- function (dat, batch, mod = NULL, mean.only = FALSE,
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
  return(bayesdata)
}
