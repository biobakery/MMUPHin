standardize.feature <- function(y,
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

Maaslin2.wrapper <- function(taxa,
                             metadata,
                             variables,
                             covariates.random = NULL,
                             directory = "./",
                             normalization = "TSS",
                             transform = "AST",
                             analysis_method = "LM",
                             ...) {
  # Specify files
  output.file <- file.path(directory, "output.txt")

  # Create temporary feature/sample/covariate names to avoid
  # Weird scenarios
  taxa.rename <- taxa
  metadata.rename <- metadata[, c(variables, covariates.random), drop = FALSE]
  taxa.names.rename <- rename.Maaslin(rownames(taxa), prefix = "T")
  sample.names.rename <- rename.Maaslin(colnames(taxa), prefix = "S")
  variables.rename <- rename.Maaslin(variables, prefix = "X")
  covariates.random.rename <- rename.Maaslin(covariates.random, prefix = "RX")
  dimnames(taxa.rename) <- list(taxa.names.rename, sample.names.rename)
  dimnames(metadata.rename) <- list(sample.names.rename,
                                    c(variables.rename,
                                      covariates.random.rename))
  # covariates.random.rename <- rename.Maaslin(NULL, prefix = "RX")

  # Run Maaslin2
  log.Maaslin <- suppressWarnings(
    capture.output(
      res.rename <- Maaslin2::Maaslin2(input_data = taxa.rename,
                                       input_metadata = metadata.rename,
                                       output = directory,
                                       min_abundance = 0,
                                       min_prevalence = 0,
                                       normalization = normalization,
                                       transform = transform,
                                       analysis_method = analysis_method,
                                       max_significance = 1,
                                       random_effects = covariates.random.rename,
                                       fixed_effects = variables.rename,
                                       standardize = FALSE,
                                       ...)$results
    ))

  cat(paste(log.Maaslin, collapse = "\n"),
      file = file.path(directory, "Maaslin.log"))
  # Read Maaslin results
  table.taxa.rename <-
    data.frame(Feature = names(taxa.names.rename),
               Feature.rename = taxa.names.rename,
               stringsAsFactors = FALSE)
  res <- list()
  for(variable in variables) {
    i.result <- dplyr::filter(res.rename, metadata == variables.rename[variable])
    i.result <- dplyr::left_join(table.taxa.rename,
                                 i.result,
                                 by = c("Feature.rename" = "feature"))
    i.result$Variable <- variable
    i.result <- dplyr::select(i.result,
                              Variable,
                              Feature = Feature,
                              Value = value,
                              Coefficient = coef,
                              N = N,
                              N.not.0 = N.not.zero,
                              P.value = pval,
                              Q.value = qval)
    i.result$Standard.error <- get.se.Maaslin(i.result$Coefficient, i.result$P.value)
    res[[variable]] <- i.result
  }

  return(res)
}

rename.Maaslin <- function(old.names, prefix) {
  if(is.null(old.names)) return(NULL)
  new.names <- paste0(prefix, seq_along(old.names))
  names(new.names) <- old.names
  return(new.names)
}

get.se.Maaslin <- function(coefficient, p) {
  ifelse(p != 1,
         abs(coefficient / qnorm(p/2)),
         NA)
}

# Maaslin.wrapper <- function(taxa,
#                             metadata,
#                             variables,
#                             covariates.random = NULL,
#                             directory = "./",
#                             ...) {
#   # Specify files
#   conf.file <- file.path(directory, "data.conf")
#   data.file.tsv <- file.path(directory, "data.tsv")
#   output.file <- file.path(directory, "output.txt")
#
#   # Create temporary feature/sample/covariate names to avoid
#   # Weird scenarios
#   taxa.rename <- taxa
#   metadata.rename <- metadata[, c(variables, covariates.random), drop = FALSE]
#   taxa.names.rename <- rename.Maaslin(rownames(taxa), prefix = "T")
#   sample.names.rename <- rename.Maaslin(colnames(taxa), prefix = "S")
#   variables.rename <- rename.Maaslin(variables, prefix = "X")
#   covariates.random.rename <- rename.Maaslin(covariates.random, prefix = "RX")
#   dimnames(taxa.rename) <- list(taxa.names.rename, sample.names.rename)
#   dimnames(metadata.rename) <- list(sample.names.rename,
#                                     c(variables.rename,
#                                       covariates.random.rename))
#   # Write Maaslin files
#   write.config(t(taxa.rename),
#                c(variables.rename,
#                  covariates.random.rename),
#                conf.file)
#   write.data(t(taxa.rename), metadata.rename,
#              c(variables.rename, covariates.random.rename),
#              data.file.tsv)
#
#   # Run the Maaslin command
#   log.Maaslin <- suppressWarnings(capture.output(Maaslin::Maaslin(strInputTSV = data.file.tsv,
#                                                                   strInputConfig = conf.file,
#                                                                   strOutputDIR = directory,
#                                                                   strRandomCovariates = covariates.random.rename,
#                                                                   dMinAbd=0,
#                                                                   dMinSamp=0,
#                                                                   dSignificanceLevel=1,
#                                                                   ...)))
#   cat(paste(log.Maaslin, collapse = "\n"),
#       file = file.path(directory, "Maaslin.log"))
#   # Read Maaslin results
#   table.taxa.rename <-
#     data.frame(Feature = names(taxa.names.rename),
#                Feature.rename = taxa.names.rename,
#                stringsAsFactors = FALSE)
#   res.rename <- read.Maaslin(directory)
#   res <- list()
#   for(variable in variables) {
#     if(variables.rename[variable] %in% names(res.rename)) {
#       i.result.rename <- res.rename[[variables.rename[variable]]]
#       i.lvls <- unique(gsub(variables.rename[variable], "", i.result.rename$Value, fixed = TRUE))
#       i.table.taxa.rename <- expand.grid(Feature.rename = taxa.names.rename,
#                                          Value = paste0(variables.rename[variable],
#                                                         i.lvls),
#                                          stringsAsFactors = FALSE)
#       i.table.taxa.rename <- dplyr::left_join(i.table.taxa.rename, table.taxa.rename,
#                                               by = "Feature.rename")
#       i.result <- dplyr::left_join(i.table.taxa.rename,
#                                    i.result.rename,
#                                    by = c("Feature.rename" = "Feature",
#                                           "Value" = "Value"))
#       i.result$Value <- gsub(paste0("^", variables.rename[variable]),
#                              variable,
#                              i.result$Value)
#       i.result$Variable <- variable
#       i.result <- i.result[, c("Variable",
#                                "Feature",
#                                "Value",
#                                "Coefficient",
#                                "N",
#                                "N.not.0",
#                                "P.value",
#                                "Q.value")]
#
#     } else {
#       i.result <- data.frame(
#         Variable = variable,
#         Feature = names(taxa.names.rename),
#         Value = NA,
#         Coefficient = NA,
#         N = NA,
#         N.not.0 = NA,
#         P.value = NA,
#         Q.value = NA
#       )
#     }
#     i.result$Standard.error <- get.se.Maaslin(i.result$Coefficient, i.result$P.value)
#     res[[variable]] <- i.result
#   }
#
#   return(res)
# }

# write.config <- function(taxa, variables, filename){
#   template <-
#     "Matrix: Metadata
#   Read_PCL_Rows: %s-%s
#
#   Matrix: Abundance
#   Read_PCL_Rows: %s-"
#
#   cat(sprintf(template,
#               variables[1],
#               variables[length(variables)],
#               colnames(taxa)[1]),
#       file = filename)
# }

# write.data <- function(taxa, metadata, variables, filename){
#   meta <- as.matrix(metadata[, variables, drop = FALSE])
#   otu <- as.matrix(taxa)
#   header <- matrix(rownames(taxa), ncol = 1, dimnames = list(NULL, "ID"))
#
#   mat <- cbind(header, meta, otu)
#
#   write.table(mat, sep = "\t", row.names = F, col.names = T, quote = F, file = filename)
# }


