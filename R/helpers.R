standardize_feature <- function(y,
                                i_design,
                                n.batch) {
  beta.hat <- solve(crossprod(i_design), crossprod(i_design, t(y)))
  grand.mean <- tcrossprod(t(beta.hat[1:n.batch, , drop = FALSE]), i_design[, 1:n.batch, drop = FALSE]) %>%
    mean()

  ## change var.pooled for ref batch
  var.pooled <- var((y - t(i_design %*% beta.hat))[1, ])
  stand.mean <- rep(grand.mean, ncol(y))
  if(ncol(i_design) > n.batch){
    stand.mean <- stand.mean +
      tcrossprod(beta.hat[-(1:n.batch), , drop = FALSE],
                 i_design[, -(1:n.batch), drop = FALSE])
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

