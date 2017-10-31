#' Wrapper for differential abundance anlysis using metagenomeSeq
#'
#' @param physeq
#' @param outcome_var
#' @param adj_var
#'
#' @return
#' @export
#'
#' @import magrittr phyloseq
#' @importFrom metagenomeSeq newMRexperiment
#' @importFrom metagenomeSeq newMRexperiment calcPosComponent
#' @examples
DA_single_metagenomeSeq <- function(physeq, model) {
  df_meta <- sample_data(physeq)
  class(df_meta) <- 'data.frame'

  mat_otu <- otu_table(physeq)@.Data
  depth <- apply(mat_otu, 2, sum)
  mat_otu_ra <- apply(mat_otu, 2, function(x) x / sum(x))

  fit_logNormal <- mat_otu %>% newMRexperiment(normFactors = depth) %>%
    fitFeatureModel2(model.matrix(model, data = df_meta))

  df_return <- data.frame(matrix(NA, ntaxa(physeq), 4))
  df_return <- data.frame(
    coef = fit_logNormal$fitZeroLogNormal$logFC,
    sd = fit_logNormal$fitZeroLogNormal$se,
    p = fit_logNormal$pvalues,
    p.adj = fit_logNormal$pvalues %>% p.adjust(method = 'BH')
  )
  colnames(df_return) <- c('coef', 'sd', 'p', 'p.adj')
  rownames(df_return) <- taxa_names(physeq)

  return(df_return)
}
