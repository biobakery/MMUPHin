library(phyloseq)
library(tidyverse)
library(curatedMetagenomicData)
studies_stool <- combined_metadata %>%
  filter(body_site %in% "stool",
         disease %in% c("healthy")) %>%
  `$`(dataset_name) %>%
  table()
studies_stool <- names(studies_stool)[studies_stool > 30]
subjects_stool <- combined_metadata$subjectID[combined_metadata$dataset_name %in% studies_stool]

all_stool <- curatedMetagenomicData("*metaphlan_bugs_list.stool", dryrun = FALSE)
all_stool <- mergeData(all_stool)
all_stool <- all_stool[, pData(all_stool)$disease %in% "healthy"]
all_stool <- all_stool[, pData(all_stool)$subjectID %in% subjects_stool]
all_stool.phyloseq <- ExpressionSet2phyloseq(all_stool)
all_stool.phyloseq <- all_stool.phyloseq %>%
  subset_taxa(is.na(Species) & !is.na(Genus))
all_stool.phyloseq <- all_stool.phylos
