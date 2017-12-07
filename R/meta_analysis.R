#' Main function for batch effect adjustment
#'
#' @param physeq combined phyloseq object with otu table and metadata across studies
#' @param study study indicator variables
#' @param fit_model formula for model to fit with Maaslin
#' @param method fit method for rma
#' @param dir_output output directory for Maaslin fit results
#'
#' @return a phyloseq object, with batch-corrected abundance data
#' @export
#' @import phyloseq magrittr Maaslin metafor
meta_fit <- function(physeq, batch, fit_model=NULL, method, dir_output="./") {
  physeq_ra <- transform_sample_counts(physeq, function(x) x / sum(x))
  mat_otu <- otu_table(physeq)@.Data
  study <- as.factor(study)
  variables_fit <- as.character(fit_model)[[2]] %>% strsplit(" + ", fixed = T) %>%
    `[[`(1)

  # fit Maaslin model in individual datasets
  for(iStudy in levels(study)) {
    df_meta_forMsln <- sample_data(physeq_ra)[study %in% iStudy, variables_fit]
    df_meta_forMsln <- data.frame(samples = rownames(df_meta_forMsln), df_meta_forMsln)
    otu_table_forMsln <- otu_table(physeq_ra)@.Data[, study %in% iStudy]
    otu_table_forMsln <- otu_table_forMsln[apply(otu_table_forMsln != 0, 1, sum) > 0, ]
    data.frame(df_meta_forMsln,
               otu_table_forMsln %>% t, check.names = F) %>%
      write.table(file = paste0(dir_output, "for_Maaslin.tsv"), sep = '\t', quote = F,
                  row.names = F)
    Maaslin(strInputTSV = paste0(dir_output, "for_Maaslin.tsv"),
            strOutputDIR = paste0(dir_output, iStudy, "/"),
            dMinSamp = 0,
            fZeroInflated = FALSE, # no zero-inflation
            strModelSelection = "none", # no model selection
            strMethod = "lm", # linear model
            strTransform = "asinsqrt",
            iLastMetadata = 1 + length(variables_fit)
    )
  }

  # random effects model for each fitted variables
  l_fit <- list()
  for (variable in variables_fit) {
    df_fit <- lapply(levels(study), function(iStudy) {
      df_fit_tmp <- read.table(paste0(dir_output, iStudy, "/",
                                      iStudy, "-", variable, ".txt"),
                               header = T, row.names = 2, sep = '\t',
                               stringsAsFactors = F)
      df_fit_tmp <- df_fit_tmp %>% select(- Value, - N.not.0)
      df_fit_tmp$Variable <- df_fit_tmp$Variable
      df_fit_tmp$Coefficient <- -df_fit$Coefficient
      df_fit_tmp$Study <- iStudy
    }) %>% Reduce("rbind", .)
    rma_fit <- rma(df_fit$Coefficient, df_fit$Var, method = method)
    l_fit[[variable]] <- rma_fit
  }

  return(l_fit)
}
