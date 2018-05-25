#' #' Main function for batch effect adjustment
#' #'
#' #' @param physeq combined phyloseq object with otu table and metadata across studies
#' #' @param study study indicator variables
#' #' @param fit_model formula for model to fit with Maaslin
#' #' @param method fit method for rma
#' #' @param dir_output output directory for Maaslin fit results
#' #'
#' #' @return a phyloseq object, with batch-corrected abundance data
#' #' @export
#' #' @import phyloseq magrittr Maaslin metafor
#' meta_fit <- function(physeq, study, fit_model=NULL, method, dir_output="./") {
#'   physeq_ra <- transform_sample_counts(physeq, function(x) x / sum(x))
#'   mat_otu <- otu_table(physeq)@.Data
#'   study <- as.factor(study)
#'   df_meta <- sample_data(physeq)
#'   class(df_meta) <- "data.frame"
#'   variables_fit <- as.character(fit_model)[[2]] %>% strsplit(" + ", fixed = T) %>%
#'     `[[`(1)
#'   fit_model_tmp <- paste0("y ~ ", paste(variables_fit, collapse = "+")) %>%
#'     as.formula()
#'
#'   # fit Maaslin model in individual datasets
#'   # for(iStudy in levels(study)) {
#'   #   df_meta_forMsln <- sample_data(physeq_ra)[study %in% iStudy, variables_fit]
#'   #   df_meta_forMsln <- data.frame(samples = rownames(df_meta_forMsln), df_meta_forMsln)
#'   #   otu_table_forMsln <- otu_table(physeq_ra)@.Data[, study %in% iStudy]
#'   #   otu_table_forMsln <- otu_table_forMsln[apply(otu_table_forMsln != 0, 1, sum) > 0, ]
#'   #   data.frame(df_meta_forMsln,
#'   #              otu_table_forMsln %>% t, check.names = F) %>%
#'   #     write.table(file = paste0(dir_output, "for_Maaslin.tsv"), sep = '\t', quote = F,
#'   #                 row.names = F)
#'   #   Maaslin(strInputTSV = paste0(dir_output, "for_Maaslin.tsv"),
#'   #           strOutputDIR = paste0(dir_output, iStudy, "/"),
#'   #           dMinSamp = 0,
#'   #           fZeroInflated = FALSE, # no zero-inflation
#'   #           strModelSelection = "none", # no model selection
#'   #           strMethod = "lm", # linear model
#'   #           strTransform = "asinsqrt",
#'   #           iLastMetadata = 1 + length(variables_fit)
#'   #   )
#'   # }
#'   l_fit <- list()
#'   for(iStudy in levels(study)) {
#'     df_meta_forMsln <- sample_data(physeq_ra)[study %in% iStudy, variables_fit]
#'     df_meta_forMsln <- data.frame(samples = rownames(df_meta_forMsln), df_meta_forMsln)
#'     otu_table_forMsln <- otu_table(physeq_ra)@.Data[, study %in% iStudy]
#'     # otu_table_forMsln <- otu_table_forMsln[apply(otu_table_forMsln != 0, 1, sum) > 0, ]
#'     otu_table_forMsln <- otu_table_forMsln %>% sqrt() %>% asin()
#'     l_fit[[iStudy]] <- lapply(rownames(otu_table_forMsln), function(feature) {
#'       y <- otu_table_forMsln[feature, ]
#'       df_meta_forMsln$y <- y
#'       df_result_tmp <- try(lm(fit_model_tmp, data = df_meta_forMsln) %>%
#'         summary %>% coef, silent = TRUE)
#'       if(class(df_result_tmp) == "try-error") {
#'         df_result_tmp <- data.frame(Estimate = NA,
#'                                     `Std. Error` = NA,
#'                                     `t value` = NA,
#'                                     `Pr(>|t|)` = NA,
#'                                     variable = colnames(model.matrix(fit_model_tmp,
#'                                                                      data = df_meta %>% data.frame(y = 7))),
#'                                     feature = feature,
#'                                     study = iStudy,
#'                                     check.names = F)
#'       }
#'       else {
#'         df_result_tmp <- data.frame(df_result_tmp,
#'                                     variable = rownames(df_result_tmp),
#'                                     feature = feature,
#'                                     study = iStudy,
#'                                     check.names = F)
#'       }
#'       return(df_result_tmp)
#'     }) %>% Reduce("rbind", .)
#'   }
#'   df_fit <- Reduce("rbind", l_fit)
#'
#'   # random effects model for each fitted variables
#'   df_return <- data.frame(beta = NULL,
#'                           se = NULL,
#'                           zval = NULL,
#'                           pval = NULL,
#'                           tau2 = NULL,
#'                           variable = NULL,
#'                           feature = NULL,
#'                           method = NULL,
#'                           beta_all = NULL,
#'                           se_all = NULL)
#'   for (iVariable in unique(df_fit$variable)) {
#'     for (iFeature in taxa_names(physeq)) {
#'       df_fit_tmp <- df_fit %>% subset(variable == iVariable &
#'                                         feature == iFeature)
#'       if(any(df_fit_tmp$`Std. Error` == 0, na.rm = T)) {
#'         df_return <- data.frame(beta = NA,
#'                                 se = NA,
#'                                 zval = NA,
#'                                 pval = NA,
#'                                 tau2 = NA,
#'                                 variable = iVariable,
#'                                 feature = iFeature,
#'                                 method = method,
#'                                 beta_all = df_fit_tmp$`Estimate`,
#'                                 se_all = df_fit_tmp$`Std. Error`,
#'                                 study = df_fit_tmp$study) %>% rbind(df_return, .)
#'       } else {
#'         rma_fit_tmp <- try(rma(df_fit_tmp$`Estimate`,
#'                            sei = df_fit_tmp$`Std. Error`,
#'                            method = method), silent = T)
#'         if(class(rma_fit_tmp) == "try-error") {
#'           df_return <- data.frame(beta = NA,
#'                                   se = NA,
#'                                   zval = NA,
#'                                   pval = NA,
#'                                   tau2 = NA,
#'                                   variable = iVariable,
#'                                   feature = iFeature,
#'                                   method = method,
#'                                   beta_all = df_fit_tmp$`Estimate`,
#'                                   se_all = df_fit_tmp$`Std. Error`,
#'                                   study = df_fit_tmp$study
#'           ) %>% rbind(df_return, .)
#'         } else {
#'           df_return <- data.frame(beta = rma_fit_tmp$beta,
#'                                   se = rma_fit_tmp$se,
#'                                   zval = rma_fit_tmp$zval,
#'                                   pval = rma_fit_tmp$pval,
#'                                   tau2 = rma_fit_tmp$tau2,
#'                                   variable = iVariable,
#'                                   feature = iFeature,
#'                                   method = method,
#'                                   beta_all = df_fit_tmp$`Estimate`,
#'                                   se_all = df_fit_tmp$`Std. Error`,
#'                                   study = df_fit_tmp$study
#'           ) %>% rbind(df_return, .)
#'         }
#'       }
#'
#'     }
#'   }
#'   return(df_return)
#' }
