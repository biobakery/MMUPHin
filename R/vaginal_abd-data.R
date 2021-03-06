#' Species level feature abundance data of two public vaginal studies
#'
#' Species level relative abundance profiles of vaginal samples in
#' the two public studies provided in
#' \code{\link[curatedMetagenomicData]{curatedMetagenomicData}}.
#'
#' @docType data
#'
#' @usage data(vaginal_abd)
#'
#' @format A feature-by-sample \code{matrix} of species-level profiles
#'
#' @references Pasolli, Edoardo, Lucas Schiffer, Paolo Manghi, Audrey Renson, 
#' Valerie Obenchain, Duy Tin Truong, Francesco Beghini et al. "Accessible, 
#' curated metagenomic data through ExperimentHub." Nature methods 14, no. 11 
#' (2017): 1023.
#'
#' @source \code{\link[curatedMetagenomicData]{curatedMetagenomicData}}
#'
#' @examples
#' data(vaginal_abd)
#' # features included
#' rownames(vaginal_abd)
#' # These are relative abundances
#' apply(vaginal_abd, 2, sum)
#' # The following were used to generate the object
#' # library(curatedMetagenomicData)
#' # library(phyloseq)
#' # datasets <- curatedMetagenomicData(
#' #   "*metaphlan_bugs_list.vagina*",
#' #   dryrun = FALSE)
#' # Construct phyloseq object from the five datasets
#' # physeq <-
#'   # Aggregate the five studies into ExpressionSet
#' #   mergeData(datasets) %>%
#'   # Convert to phyloseq object
#' #   ExpressionSet2phyloseq() %>%
#'   # Subset features to species
#' #   subset_taxa(!is.na(Species) & is.na(Strain)) %>%
#'   # Normalize abundances to relative abundance scale
#' #   transform_sample_counts(function(x) x / sum(x)) %>%
#'   # Filter features to be of at least 1e-5 relative abundance in two samples
#' #   filter_taxa(kOverA(2, 1e-5), prune = TRUE)
#' # vaginal_abd <- otu_table(physeq)@.Data
"vaginal_abd"