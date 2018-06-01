adjust.batch <- function(feature.count,
                         batch,
                         formula.adj = NULL,
                         data = NULL,
                         zero.inflation = TRUE,
                         pseudo.count = 0.5,
                         filtering = TRUE,
                         verbose = TRUE) {

  ## Ensure feature.count is a count table
  otu.count <- as.matrix(otu.count)
  if(any(otu.count < 0, na.rm = TRUE))
    stop("Found negative values in the feature count table!")

  ## Filter features
  features.count <- feature.count[, ]

  ## Construct covariate adjustment model matrix
  mod <- model.matrix(formula.adj, data)

  ## Filter samples with missing values
  if(ncol(feature.count) != length(batch) |
     ncol(feature.count) != nrow(mod)
     )


  ## Check for missing values
  if(any(is.na(feature.count)))
    stop("Found missing values in the feature table!")
  if(any(is.na(batch)))
    stop("Found missing values in the batch variable!")
  if(any(is.na(mod)))
    stop("Found missing values in the covariates!")



  ## Transform data for ComBat fit
  feature.count <- as.matrix(feature.count)
  log.data <- log(apply(feature.count + pseudo.count,
                        2,
                        function(x) x / sum(x, na.rm - TRUE)))

  ## Summarise batch variable
  batch <- as.factor(batch)
  n.batch <- nlevels(batch)
  message("Found ", n.batch, " batches.")
  if(n.batch < 2)
    stop("!")

  batchmod <- model.matrix(~ -1 + batch)

  #
  # batches <- list()
  # for (i in 1:n.batch) {
  #   batches[[i]] <- which(batch == levels(batch)[i])
  # } # list of samples in each batch
  # n.batches <- sapply(batches, length)
  # n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)

  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[,!check])

  ## Number of covariates or covariate levels
  message("Adjusting for ", ncol(design)-ncol(batchmod), ' covariate(s) or covariate level(s).')

  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded.")
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat.")
      }
  }

  ## Identify data features to adjust for
  ind_data <- matrix(TRUE, nrow(feature.count), ncol(feature.count)) # which feature table values are zero
  ind_gamma <- matrix(TRUE, nrow(feature.count), n.batch) # which batch X bug pairs are present
  ind_feature <- rep(TRUE, length = nrow(feature.count)) # which features have enough information for adjustment
  ind_mod <- rep(TRUE, length = ncol(design) - n.batch) # All covariates are included
  if(!zero.inflation) {
    for(i_feature in 1:nrow(feature.count)) {
      ind_data[i_feature, ] <- feature.count[i_feature, ] != 0
      ind_gamma[i_feature, ] <- apply(batchmod[ind_data[i_feature, ], , drop = FALSE] == 1, 2, any)
      ind_feature[i_feature] <- {
        sum(ind_data[i_feature, ]) > ncol(design) - sum(!ind_gamma[i_feature, ]) && # should have at least one more sample than covariate
          sum(ind_gamma[i_feature, ]) > 1 && # and at least two batches for correction
          qr(design[ind_data[i_feature, ],
                    c(ind_gamma[i_feature, ], ind_mod)]
          )$rank == ncol(design) - sum(!ind_gamma[i_feature, ]) # and the design matrix should be full rank
      }
      if(!ind_feature[i_feature]) ind_gamma[i_feature, ] <- FALSE
    }
  }
  message("Adjusting for ", sum(ind_feature), ' feature(s).')

  ## Standardize data across features
  message("Standardizing Data across features.")
  s.data <- log.data
  l_standFit <- list()
  for(i_feature in 1:nrow(s.data)) {
    if(ind_feature[i_feature]) {
      i_ind_sample <- ind_data[i_feature, ]
      standFit <- standardize_feature(y = s.data[i_feature, i_ind_sample, drop = FALSE],
                                      i_design = design[i_ind_sample, c(ind_gamma[i_feature, ], ind_mod), drop = FALSE],
                                      n.batch = sum(ind_gamma[i_feature, ]))
      s.data[i_feature, i_ind_sample] <- standFit$y.stand
      l_standFit[[i_feature]] <- standFit
    } else l_standFit[[i_feature]] <- NULL
  }

  ## Sanity check
  # (1:nrow(s.data))[ind_feature] %>%
  #   sapply(function(i_feature) mean(s.data[i_feature,
  #                                          ind_data[i_feature, ] &
  #                                            apply(batchmod[, ind_gamma[i_feature, ], drop = FALSE] == 1,
  #                                                  1,
  #                                                  any)
  #                                          ])
  #   )
  # (1:nrow(s.data))[ind_feature] %>%
  #   sapply(function(i_feature) sd(s.data[i_feature,
  #                                        ind_data[i_feature, ] &
  #                                          apply(batchmod[, ind_gamma[i_feature, ], drop = FALSE] == 1,
  #                                                1,
  #                                                any)
  #                                        ])
  #   )

  ## Get regression batch effect parameters
  message("Estimating EB prior parameters.")
  gamma.hat <-
    delta.hat <-
    gamma.star <-
    delta.star <-
    matrix(NA, nrow = nrow(s.data), ncol = n.batch)
  for(i_feature in 1:nrow(s.data)) {
    if(ind_feature[i_feature]) {
      i_ind_sample <- ind_data[i_feature, ]
      abd_batch <- s.data[i_feature, i_ind_sample] * batchmod[i_ind_sample, ind_gamma[i_feature, ], drop = FALSE]
      i_gamma <- apply(abd_batch, 2, function(x) mean(setdiff(x, 0)))
      i_delta <- apply(abd_batch, 2, function(x) var(setdiff(x, 0)))
      i_delta[is.na(i_delta)] <- 1
      gamma.hat[i_feature, ind_gamma[i_feature, ]] <- i_gamma
      delta.hat[i_feature, ind_gamma[i_feature, ]] <- i_delta
    }
  }

  t2 <- apply(gamma.hat, 2, var, na.rm = TRUE)
  gamma.bar <- apply(gamma.hat, 2, mean, na.rm = TRUE)
  a.prior <- apply(delta.hat, 2, aprior, na.rm = TRUE) # FIXME
  b.prior <- apply(delta.hat, 2, bprior, na.rm = TRUE) # FIXME

  if(any(apply(!is.na(gamma.hat), 2, sum) < 2) |
     any(apply(!is.na(delta.hat), 2, sum) < 2))
    stop("One batch has only one feature with valid parameter estimate!")

  ## Making adjustments
  message("Performing shrinkage adjustments on per-batch parameters.")
  results <- lapply(1:n.batch, function(i_batch) {
    i_s.data <- s.data
    i_s.data[!ind_data] <- NA # set all zeros to NA
    i_s.data[, !as.logical(batchmod[, i_batch])] <- NA # set all other batches to NA
    i_s.data[!ind_gamma[, i_batch], ] <- NA # set features not adjustable to NA
    temp <- it.sol(sdat = i_s.data,
                   g.hat = gamma.hat[, i_batch],
                   d.hat = delta.hat[, i_batch],
                   g.bar = gamma.bar[i_batch],
                   t2 = t2[i_batch],
                   a = a.prior[i_batch],
                   b = b.prior[i_batch])
    gamma.star <- temp[1, ]
    delta.star <- temp[2, ]
    list(gamma.star=gamma.star, delta.star=delta.star)
  })
  for (i_batch in 1:n.batch) {
    gamma.star[, i_batch] <- results[[i_batch]]$gamma.star
    delta.star[, i_batch] <- results[[i_batch]]$delta.star
  }

  # visualize shrinkage
  # tb_plot <- tibble::tibble(gamma = c(as.vector(gamma.hat), as.vector(gamma.star)),
  #                           delta = c(as.vector(delta.hat), as.vector(delta.star)),
  #                           feature_group = rep((1:length(gamma.hat)), 2),
  #                           adjustment = rep(c("before", "after"), each = length(gamma.hat)))
  # tb_plot %>% ggplot(aes(x = gamma, y = delta)) +
  #   geom_point(aes(color = adjustment)) +
  #   geom_line(aes(group = feature_group)) +
  #   theme_bw()

  adj.data <- s.data
  for (i_batch in 1:n.batch){
    i_ind_feature <- !is.na(gamma.star[, i_batch]) & !is.na(delta.star[, i_batch])
    for(i_feature in 1:nrow(adj.data)) {
      if(i_ind_feature[i_feature]) {
        i_ind_sample <- ind_data[i_feature, ] & as.logical(batchmod[, i_batch])
        adj.data[i_feature, i_ind_sample] <-
          (adj.data[i_feature, i_ind_sample] -
             gamma.star[i_feature, i_batch]) /
          sqrt(delta.star[i_feature, i_batch])
      }
    }
  }

  message("Performing batch adjustments.")
  for(i_feature in 1:nrow(adj.data)) {
    if(ind_feature[i_feature]) {
      ind_sample <- ind_data[i_feature, ]
      standFit <- l_standFit[[i_feature]]
      adj.data[i_feature, ind_sample] <-
        adj.data[i_feature, ind_sample] * sqrt(standFit$var.pooled) + standFit$stand.mean
    }
  }

  if(any(is.na(adj.data)))
    stop("There are missing values in the adjusted data!")

  adj.data <- exp(adj.data)
  adj.data.ra <- apply(adj.data, 2, function(x) x / sum(x))
  adj.data.ra[feature.count == 0] <- 0
  adj.data.count <- t(t(adj.data.ra) * apply(feature.count, 2, sum))

  return(adj.data.count)
}

