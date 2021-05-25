# These are notes that were originally in the vignette as commented code. They were removed because it made the vignette
# harder to read. Some were integrated into the main text of the vignette.

# CRC_abd is the feature (species) abundance matrix. Rows are features and 
# columns are samples.

# CRC_meta is the metadata data frame. Columns are samples.

# A total of five studies are included 

# The following were used to access and format the two objects
# library(curatedMetagenomicData)
# library(phyloseq)
# library(genefilter)
# datasets <- curatedMetagenomicData(
#   c("FengQ_2015.metaphlan_bugs_list.stool"  ,
#     "HanniganGD_2017.metaphlan_bugs_list.stool",
#     "VogtmannE_2016.metaphlan_bugs_list.stool",
#     "YuJ_2015.metaphlan_bugs_list.stool",
#     "ZellerG_2014.metaphlan_bugs_list.stool"),
#   dryrun = FALSE)
# # Construct phyloseq object from the five datasets
# physeq <-
#   # Aggregate the five studies into ExpressionSet
#   mergeData(datasets) %>%
#   # Convert to phyloseq object
#   ExpressionSet2phyloseq() %>%
#   # Subset samples to only CRC and controls
#   subset_samples(study_condition %in% c("CRC", "control")) %>%
#   # Subset features to species
#   subset_taxa(!is.na(Species) & is.na(Strain)) %>%
#   # Normalize abundances to relative abundance scale
#   transform_sample_counts(function(x) x / sum(x)) %>%
#   # Filter features to be of at least 1e-5 relative abundance in five samples
#   filter_taxa(kOverA(5, 1e-5), prune = TRUE)
# CRC_abd <- otu_table(physeq)@.Data
# CRC_meta <- data.frame(sample_data(physeq))
# CRC_meta$studyID <- factor(CRC_meta$studyID)

# The function call indicates for adjust_batch to correct for the effect
# of the batch variable, studyID, while controlling for the effect of the 
# disease variable, study_condition. Many additional options are available 
# through the control parameter, here we specify verbose=FALSE to avoid 
# excessive messages, although they can often be helpful in practice!

# Note that adjust_batch returns a list of more than one components, and 
# feature_abd_adj is the corrected feature abundance matrix. See 
# help(adjust_batch) for the meaning of other components.

# First subset both feature abundance table and metadata to only control samples

# discrete_discover takes as input sample-by-sample dissimilarity measurements 
# rather than abundance table. The former can be easily computed from the 
# latter with existing R packages.

# By default, fit_discrete evaluates cluster numbers 2-10
# By default, fit_discrete evaluates cluster numbers 2-10

# library(curatedMetagenomicData)
# library(phyloseq)
# datasets <- curatedMetagenomicData(
#   "*metaphlan_bugs_list.vagina*",
#   dryrun = FALSE)
# # Construct phyloseq object from the five datasets
# physeq <-
#   # Aggregate the five studies into ExpressionSet
#   mergeData(datasets) %>%
#   # Convert to phyloseq object
#   ExpressionSet2phyloseq() %>%
#   # Subset features to species
#   subset_taxa(!is.na(Species) & is.na(Strain)) %>%
#   # Normalize abundances to relative abundance scale
#   transform_sample_counts(function(x) x / sum(x)) %>%
#   # Filter features to be of at least 1e-5 relative abundance in two samples
#   filter_taxa(kOverA(2, 1e-5), prune = TRUE)
# vaginal_abd <- otu_table(physeq)@.Data
# vaginal_meta <- data.frame(sample_data(physeq))
# vaginal_meta$studyID <- factor(vaginal_meta$studyID)

