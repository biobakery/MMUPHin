#' Sample metadata of five public CRC studies
#'
#' Metadata information of CRC and control patients in
#' the five public studies used in Thomas et al. (2019). These were accessed
#' through \code{\link[curatedMetagenomicData]{curatedMetagenomicData}}.
#'
#' @docType data
#'
#' @usage data(CRC_meta)
#'
#' @format A \code{data.frame} of per-sample metadata information
#'
#' @references Thomas, Andrew Maltez, Paolo Manghi, Francesco Asnicar, 
#' Edoardo Pasolli, Federica Armanini, Moreno Zolfo, Francesco Beghini et al. 
#' "Metagenomic analysis of colorectal cancer datasets identifies cross-cohort 
#' microbial diagnostic signatures and a link with choline degradation." Nature 
#' medicine 25, no. 4 (2019): 667.
#'
#' @source \code{\link[curatedMetagenomicData]{curatedMetagenomicData}}
#'
#' @examples
#' data(CRC_meta)
#' # has CRC and control samples across five studies
#' table(CRC_meta$studyID, CRC_meta$study_condition)
#' # The following were used to generate the object
#' # library(curatedMetagenomicData)
#' # library(phyloseq)
#' # library(genefilter)
#' # datasets <- curatedMetagenomicData(
#' #   c("FengQ_2015.metaphlan_bugs_list.stool"  ,
#' #     "HanniganGD_2017.metaphlan_bugs_list.stool",
#' #     "VogtmannE_2016.metaphlan_bugs_list.stool",
#' #     "YuJ_2015.metaphlan_bugs_list.stool",
#' #     "ZellerG_2014.metaphlan_bugs_list.stool"),
#' #   dryrun = FALSE)
#' # Construct phyloseq object from the five datasets
#' # physeq <-
#'     # Aggregate the five studies into ExpressionSet
#' #   mergeData(datasets) %>%
#'     # Convert to phyloseq object
#' #   ExpressionSet2phyloseq() %>%
#'     # Subset samples to only CRC and controls
#' #   subset_samples(study_condition %in% c("CRC", "control")) %>%
#'     # Subset features to species
#' #   subset_taxa(!is.na(Species) & is.na(Strain)) %>%
#'     # Normalize abundances to relative abundance scale
#' #   transform_sample_counts(function(x) x / sum(x)) %>%
#'     # Filter features to be of at least 1e-5 relative abundance in five 
#'     # samples
#' #   filter_taxa(kOverA(5, 1e-5), prune = TRUE)
#' # CRC_meta <- data.frame(sample_data(physeq))
#' # CRC_meta$studyID <- factor(CRC_meta$studyID)
"CRC_meta"