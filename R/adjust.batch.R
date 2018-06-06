
#' Zero-inflated Empirical Bayes Adjustment of Batch Effect in Feature Abundance Data
#'
#' @param feature.count Feature x sample matrix of feature abundance (counts)
#' @param batch Variable indicating batch membership
#' @param formula.adj Formula indicating additional covariate to adjust for
#' in the batch correction model
#' @param data.adj Data frame for covarites to adjust for.
#' @param zero.inflation Flag for whether or not a zero-inflated model should be run.
#' Default to TRUE (zero-inflated model). If set to FALSE then vanilla ComBat
#' (with parametric adjustment) will be performed.
#' @param pseudo.count Pseudo count to add to the count table before log-transformation. Default
#' to 0.5.
#' @param diagnostics Flag for whether or not diagnostic figures should be generated.
#' Deafault to TRUE.
#' @param verbose Flag for whether or not verbose modelling information should be printed.
#' Default to yes.
#'
#' @return A feature x sample matrix of adjusted feature abundance.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
adjust.batch <- function(feature.count,
                         batch,
                         formula.adj = NULL,
                         data.adj = NULL,
                         zero.inflation = TRUE,
                         pseudo.count = 0.5,
                         diagnostics = TRUE,
                         verbose = TRUE) {

  ## Ensure data formatts are as expected
  feature.count <- as.matrix(feature.count)
  if(any(feature.count < 0, na.rm = TRUE))
    stop("Found negative values in the feature table!")
  if(any(is.na(feature.count)))
    stop("Found missing values in the feature table!")
  batch <- as.factor(batch)
  if(any(is.na(batch))) stop("Found missing values in the batch variable!")
  data.adj <- as.data.frame(data.adj)

  ## Data dimensions need to agree with each other
  if(ncol(feature.count) != length(batch)) {
    stop("Dimensions of feature table and batch variable do not agree!")
     } else if(!is.null(formula.adj)) {
    if(ncol(feature.count) != nrow(data.adj))
      stop("Dimensions of feature table and covariate table do not agree!")
     }

  ## If specified, construct covariate adjustment table
  mod <- NULL
  ind.sample <- rep(TRUE, ncol(feature.count))
  if(!is.null(formula.adj)) {
    mod <- model.matrix(formula.adj,
                        model.frame(~., data.adj, na.action = na.pass))
    if(any(is.na(mod))) warning("Found missing values in the covariate table; only fully observed records will be adjusted.")
    ind.sample <- !apply(is.na(mod), 1, any)
  }

  ## Construct batch adjustment model
  n.batch <- length(unique(batch[ind.sample]))
  if(n.batch < 2) stop("Batch variable has only one level!")
  if(verbose) message("Found ", n.batch, " batches")
  batchmod <- model.matrix(~ -1 + batch)

  ## Construct design matrix
  design <- cbind(batchmod, mod)
  ## Check for intercept, or covariates not meaningful after missing value removal, remove if present
  check <- apply(design[ind.sample, ], 2, function(x) all(x == 1) | all(x == 0))
  design <- design[, !check]
  batchmod <- design[, 1:n.batch]
  ## Number of covariates or covariate levels
  if(verbose) message("Adjusting for ", ncol(design) - n.batch, " covariate(s) or covariate level(s).")

  ## Transform data for ComBat fit
  lib.size <- apply(feature.count, 2, sum)
  if(all(lib.size <= 1)) {
    warning("Feature table appears to be on the relative abundance scale. ",
            "Setting pseudo.count to be min(feature.count)/2! ")
    pseudo.count <- min(setdiff(feature.count, 0)) / 2
  }
  log.data <- log(apply(feature.count + pseudo.count,
                        2,
                        function(x) x / sum(x)))

  ## Check if the design is confounded
  if(qr(design[ind.sample, ])$rank < ncol(design[ind.sample, ])) {
      if(qr(design[ind.sample, -c(1:n.batch)])$rank <
         ncol(design[ind.sample, -c(1:n.batch)])) {
        stop("The covariates are confounded!")
      } else {
        stop("At least one covariate is confounded with batch!")
      }
  }

  ## Identify data to adjust for
  ind.data <- matrix(TRUE, nrow(feature.count), ncol(feature.count)) # which non-missing feature table values are zero
  ind.data[, !ind.sample] <- FALSE
  ind.gamma <- matrix(TRUE, nrow(feature.count), n.batch) # which feature x batch pairs are adjustable
  ind.mod <- rep(TRUE, ncol(design) - n.batch) # covariates are always adjusted for
  if(zero.inflation) {
    ind.data[feature.count == 0] <- FALSE
    for(i.feature in 1:nrow(feature.count)) {
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
  ## Batch has to have more than one adjustable feature to make EB estimates
  ind.gamma[, apply(ind.gamma, 2, sum) < 2] <- FALSE # FIXME?
  ind.feature <- apply(ind.gamma, 1, any)
  if(verbose) message("(After filtering) adjusting for ", sum(ind.feature), " feature(s).")

  ## Standardize data across features
  if(verbose) message("Standardizing Data across features.")
  s.data <- log.data
  l.standFit <- list()
  for(i.feature in 1:nrow(s.data)) {
    if(ind.feature[i.feature]) {
      i.design <- design[ind.data[i.feature, ],
                         c(ind.gamma[i.feature, ], ind.mod),
                         drop = FALSE] # FIXME
      ## For debugging, this shouldn't happen
      if(nrow(i.design) <= 1 | ncol(i.design) <= 1) stop("Something wrong happened!" ) # FIXME
      standFit <- standardize.feature(y = s.data[i.feature, ind.data[i.feature, ]],
                                      i.design = i.design,
                                      n.batch = sum(ind.gamma[i.feature, ]))
      s.data[i.feature, ind.data[i.feature, ]] <- standFit$y.stand
      l.standFit[[i.feature]] <- standFit
    } else l.standFit[[i.feature]] <- NULL
  }

  ## Estimate per-batch mean/variance
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
      ## For debugging, this shouldn't happen
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

  ## Shrink per-batch mean/variance
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

  ## Adjust the data
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

  adj.data <- exp(adj.data)
  adj.data[feature.count == 0] <- 0
  adj.data <- apply(adj.data, 2, function(x) x / sum(x))
  adj.feature.count <- t(t(adj.data) * lib.size)
  dimnames(adj.feature.count) <- dimnames(feature.count)

  ## If required, generate diagnostic plots
  if(diagnostics) {
    df.plot <- data.frame(gamma.hat = gamma.hat %>% as.vector,
                          gamma.star = gamma.star %>%  as.vector)
    plot1 <- ggplot2::ggplot(df.plot, aes(x = gamma.hat, y = gamma.star)) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(intercept = 0, slope = 1) +
      ggplot2::theme_bw() +
      ggtitle("Shrinkage of batch mean parameters")

    df.ra <- map2_df(
      list(feature.count, adj.feature.count),
      c("Original", "Adjusted"),
      function(i.count.table, i.normalization) {
        i.count.table %>%
          apply(2, function(x) x / sum(x)) %>%
          t() %>%
          as.data.frame() %>%
          tibble::rownames_to_column("Sample") %>%
          tidyr::gather(key = "feature",
                        value = "relative abundance",
                        -Sample) %>%
          dplyr::mutate(Adjustment = i.normalization)
      }
    ) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(Adjustment = Adjustment %>%
                      factor(levels = c("Original", "Adjusted")))
    df.batch <- data.frame(batch = batch, Sample = colnames(feature.count),
                           stringsAsFactors = FALSE)
    df.plot <- df.ra %>%
      dplyr::left_join(df.batch, by = "Sample") %>%
      dplyr::group_by(feature) %>%
      dplyr::mutate(`Overall Mean` = mean(`relative abundance`[Adjustment == "Original"])) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Adjustment, feature, batch, `Overall Mean`) %>%
      dplyr::summarise(`Batch Mean` = mean(`relative abundance`)) %>%
      dplyr::ungroup()
    plot2 <- df.plot %>%
      ggplot2::ggplot(aes(x = `Overall Mean`,
                          y = `Batch Mean`)) +
      ggplot2::geom_point(aes(group = paste0(Adjustment, feature),
                              color = Adjustment),
                          position = position_dodge(width = max(df.plot$`Overall Mean`) / 100)) +
      ggplot2::geom_line(aes(group = paste0(Adjustment, feature),
                             color = Adjustment),
                         position = position_dodge(width = max(df.plot$`Overall Mean`) / 100)) +
      ggplot2::geom_abline(intercept = 0, slope = 1) +
      ggplot2::scale_color_manual(values = c("Original" = "black", "Adjusted" = "red")) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position=c(0, 1),
                     legend.justification=c(0, 1),
                     legend.direction="horizontal") +
      ggplot2::ggtitle("Batch mean relative abundance for all features")
    plot <- cowplot::plot_grid(plot1, plot2, nrow = 1)
    print(plot)
  }

  return(adj.feature.count)
}

