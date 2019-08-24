#' Construct a non-intercept design model matrix given a metadata data frame
#'
#' @param data metadata data frame
#'
#' @return design matrix
construct_design <- function(data) {
  # Returns NULL if data is NULL. This happens if the covariate data frame is
  # NULL (when no covariates are provided)
  if(is.null(data)) return(NULL)

  # Construct the matrix using all variables in data but excluding the intercept
  # term (because batch dummy variables encompass the intercept term)
  model.matrix(~ . - 1, data = data)
}

#' Centralize (by design matrix) and standardize (by pooled variance across all batches)
#' feature abundances for empirical Bayes fit
#' @param y vector of non-zero abundance of a single feature (if zero-inflated is true).
#' @param i.design design matrix for the feature; samples with zeros are taken out
#' (if zero-inflated is true).
#' @param n.batch number of batches in the data.
#'
#' @return a list with component: y.stand for vector of centralized and standardized feature
#' abundance, and stand.mean/varpooled for the location and scale factor (these are used
#' later to back transform the batch-shrinked feature abundance).
standardize_feature <- function(y,
                                i.design,
                                n.batch) {
  beta.hat <- solve(crossprod(i.design),
                    crossprod(i.design, y))
  grand.mean <- mean(i.design[, 1:n.batch] %*%
                       beta.hat[1:n.batch, ])

  ## change var.pooled for ref batch
  var.pooled <- var(y - (i.design %*% beta.hat)[, 1])
  stand.mean <- rep(grand.mean, length(y))
  if(ncol(i.design) > n.batch){
    stand.mean <- stand.mean +
      (i.design[, -(1:n.batch), drop = FALSE] %*%
         beta.hat[-(1:n.batch), ])[, 1]
  }
  y.stand <- (y - stand.mean) / sqrt(var.pooled)
  return(list(y.stand = y.stand, stand.mean = stand.mean, var.pooled = var.pooled))
}

aprior <- function(gamma.hat, na.rm = FALSE) {
  m <- mean(gamma.hat, na.rm = na.rm)
  s2 <- var(gamma.hat, na.rm = na.rm)
  (2*s2 + m^2) / s2
}

bprior <- function(gamma.hat, na.rm = FALSE){
  m <- mean(gamma.hat, na.rm = na.rm)
  s2 <- var(gamma.hat, na.rm = na.rm)
  (m*s2 + m^3) / s2
}

it.sol  <- function(sdat,
                    g.hat,
                    d.hat,
                    g.bar,
                    t2,
                    a,
                    b,
                    conv=.0001){
  n <- rowSums(!is.na(sdat))
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- rowSums((sdat - g.new %*% t(rep(1,ncol(sdat))))^2, na.rm=TRUE)
    sum2[sum2 == 0] <- NA
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new-g.old) / g.old, abs(d.new-d.old) / d.old, na.rm=TRUE)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  ## cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

postmean <- function(g.hat,g.bar,n,d.star,t2){
  (t2*n*g.hat + d.star*g.bar) / (t2*n + d.star)
}

postvar <- function(sum2,n,a,b){
  (.5*sum2 + b) / (n/2 + a - 1)
}
