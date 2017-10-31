#' Title
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

#' Title
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

#' Title
#'
#' @param obj
#' @param mod
#' @param coef
#' @param B
#' @param szero
#' @param spos
#'
#' @return
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



#' Title
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
