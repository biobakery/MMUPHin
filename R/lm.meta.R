#' Main function for batch effect adjustment
#'
#' @param feature.count Feature x sample matrix of feature abundance
#' @param batch Variable indicating batch membership
#' @param exposure Primary exposure of interest
#' @param formula.adj Formula indicating additional covariate to adjust for
#' in the batch correction model
#' @param data Data frame for covarites to adjust for.
#' @param zero.inflation Flag for whether or not a zero-inflated model should be run.
#' Default to TRUE (zero-inflated model). If set to FALSE then vanilla ComBat
#' (with parametric adjustment) will be performed.
#' @param pseudo.count Pseudo count to add to the count table before log-transformation. Default
#' to 0.5.
#' @param forest Flag for whether or not forest plots figures should be generated.
#' Deafault to TRUE.
#' @param directory Directory for generated forest plots
#' @param verbose Flag for whether or not verbose modelling information should be printed.
#' Default to yes.
#' @return a long-format data.frame of main exposure's individual and aggregated effect sizes
#' @importFrom magrittr %>%
#' @export
#'
lm.meta <- function(feature.count,
                    batch,
                    exposure,
                    formula.adj = NULL,
                    re.exposure = FALSE,
                    data = NULL,
                    zero.inflation = TRUE, # for Maaslin
                    pseudo.count = 0.5,
                    forestplots = TRUE,
                    directory = "./",
                    verbose = TRUE,
                    ...) {
  ## Ensure data formatts are as expected
  feature.count <- as.matrix(feature.count)
  if(is.null(rownames(feature.count)))
    rownames(feauture.count) <- paste0("Feature", 1:nrow(feature.count))
  if(any(feature.count < 0, na.rm = TRUE))
    stop("Found negative values in the feature table!")
  if(any(is.na(feature.count)))
    stop("Found missing values in the feature table!")
  batch <- as.factor(batch)
  if(any(is.na(batch))) stop("Found missing values in the batch variable!")
  data <- as.data.frame(data)

  ## Data dimensions need to agree with each other
  if(ncol(feature.count) != length(batch)) {
    stop("Dimensions of feature table and batch variable do not agree!")
  } else if(ncol(feature.count) != length(exposure)) {
    stop("Dimensions of feature table and exposure variable do not agree!")
  }else if(ncol(feature.count) != nrow(data))
      stop("Dimensions of feature table and covariate table do not agree!")

  ## Check for RE components in the formulae

  ## If specified, construct covariate adjustment table
  mod <- NULL
  if(!is.null(formula.adj)) {
    mod <- model.matrix(formula.adj,
                        model.frame(~., data, na.action = na.pass))
    check <- apply(mod, 2, function(x) all(x == 1, na.rm = TRUE) |
                     all(x == 0, na.rm = TRUE))
    mod <- mod[, !check]
  }

  ## Construct exposure model
  expomod <- model.matrix(~ exposure)[, -1, drop = FALSE]
  nlvl.expo <- ncol(expomod)
  names.expo <- colnames(expomod)

  ## Construct design matrix
  design <- cbind(expomod, mod)
  ind.sample <- rep(TRUE, ncol(feature.count))
  if(any(is.na(design))) warning("Found missing values in the covariate table; only fully observed records will be included.")
  ind.sample <- !apply(is.na(design), 1, any)
  if(qr(design[ind.sample, , drop = FALSE])$rank < ncol(design[ind.sample, , drop = FALSE]))
    stop("The covariates are confounded!")

  ## Identify batch variables
  n.batch <- length(unique(batch[ind.sample]))
  if(n.batch < 2) stop("Batch variable has only one level!")
  if(verbose) message("Found ", n.batch, " batches")
  lvl.batch <- unique(batch) %>% setdiff(NA)

  ## Transform sample counts
  lib.size <- apply(feature.count, 2, sum)
  if(all(lib.size <= 1)) {
    warning("Feature table appears to be on the relative abundance scale. ",
            "Setting pseudo.count to be min(feature.count)/2! ")
    pseudo.count <- min(setdiff(feature.count, 0)) / 2
  }
  log.data <- log(apply(feature.count + pseudo.count,
                        2,
                        function(x) x / sum(x)))

  ## Fit individual LM models
  if(verbose) message("Fitting individual regressions.")
  lm.results <- levels(batch) %>%
    purrr::map_df(function(i.batch) {
      i.log.data <- log.data[, batch == i.batch]
      i.design <- design[batch == i.batch, ]
      i.check <- apply(i.design, 2, function(x) all(x == 1, na.rm = TRUE) |
                       all(x == 0, na.rm = TRUE))
      i.design <- i.design[, !i.check]
      rownames(i.log.data) %>%
        purrr::map_df(function(i.feature) {
          i.result <- data.frame(feature = i.feature,
                                 exposure = names.expo,
                                 tau2 = NA,
                                 I2 = NA,
                                 stringsAsFactors = FALSE)
          i.data <- data.frame(y = i.log.data[i.feature, ],
                               i.design)
          i.lm.fit <- try(lm(y ~ ., data = i.data) %>%
                          summary %>% coef,
                        silent = TRUE)
          if(class(i.lm.fit) != "try-error") {
            i.lm.fit <- i.lm.fit %>%
              data.frame() %>%
              tibble::rownames_to_column(var = "exposure")
            i.result <- i.result %>%
              dplyr::left_join(i.lm.fit, by = "exposure") %>%
              dplyr::mutate(beta = Estimate,
                     sd = Std..Error,
                     p = Pr...t..)
          } else i.result <- i.result %>%
            dplyr::mutate(beta = NA, sd = NA, p = NA)
          return(i.result %>%
                  dplyr::select(feature, exposure, beta, sd, p, tau2, I2))
        }) %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(batch = i.batch)
    }) %>%
    dplyr::bind_rows()


  # Fit fixed/random effects models
  if(verbose) message("Fitting meta-analysis model.")
  meta.results <- lm.results %>%
    dplyr::filter(!is.na(beta),
                  !is.na(sd),
                  sd != 0) %>%
    dplyr::group_by(feature, exposure) %>%
    dplyr::filter(n() > 2) %>%
    ungroup()
  meta.fit <- unique(meta.results$feature) %>%
    purrr::map_df(function(i.feature) {
      unique(meta.results$exposure) %>%
        purrr::map_df(function(i.exposure) {
          i.result <- data.frame(feature = i.feature,
                                 exposure = i.exposure,
                                 stringsAsFactors = FALSE)
          i.data <- meta.results %>%
            dplyr::filter(feature %in% i.feature,
                          exposure %in% i.exposure)
          i.rma.fit <- try(metafor::rma.uni(yi = i.data$beta,
                                            sei = i.data$sd,
                                            method = method),
                           silent = TRUE)
          if(all(class(i.rma.fit) != "try-error")) {
            i.result <- i.result %>%
              mutate(beta = i.rma.fit$beta,
                     sd = i.rma.fit$se,
                     p = i.rma.fit$pval,
                     tau2 = i.rma.fit$tau2,
                     I2 = i.rma.fit$I2)
          } else {
            i.result <- i.result %>%
              mutate(beta = NA,
                     sd = NA,
                     p = NA,
                     tau2 = NA,
                     I2 = NA)
          }
          return(i.result)
        }) %>% dplyr::bind_rows()
    }) %>% dplyr::bind_rows() %>%
    dplyr::mutate(batch = "Meta-Analysis Model")
  results <- meta.results %>%
    dplyr::mutate(p.adj = NA) %>%
    dplyr::bind_rows(
      meta.fit %>%
        dplyr::mutate(p.adj = p.adjust(p, method = "fdr"))
    )

  if(verbose) message("After correction, found ", sum(results$p.adj, na.rm = TRUE),
                      " significant feature ~ exposure associations.")

  # if(forestplots & sum(results$p.adj, na.rm = TRUE) > 0) {
  #   pdf(paste0(directory, "forestplots.pdf", width = 8, height = 1.5*nlvl.expo))
  #   meta.fit <- for(i.feature in unique(results$feature)) {
  #     for(i.exposure in unique(results$exposure)) {
  #       i.results <- results %>%
  #         dplyr::filter(feature %in% i.feature,
  #                       exposure %in% i.exposure)
  #       if(any(i.results$p.adj < 0.05, na.rm = TRUE))
  #         i.results <- i.results %>%
  #           dplyr::mutate(batch = factor(batch,
  #                                        levels = c("Meta-Analysis Model",
  #                                                   setdiff(batch, "Meta-Analysis Model"))
  #           )
  #           )
  #         i.plot <- i.results %>%
  #           ggplot2::ggplot(ggplot2::aes(x = beta,
  #                                        y = batch,
  #                                        xmin = beta - 1.96*sd,
  #                                        xmax = beta + 1.96*sd)) +
  #           ggplot2::geom_errorbarh(alpha=0.5, color="black", height = 0.3) +
  #           ggplot2::geom_point() +
  #           ggplot2::theme_bw() +
  #           ggplot2::xlab("Effect Size") +
  #           ggplot2::ylab("Study") +
  #           ggplot2::ggtitle(paste0(i.feature, " ~ ", i.exposure))
  #         print(i.plot)
  #     }
  #   }
  #   dev.off()
  # }

  return(results)
}
