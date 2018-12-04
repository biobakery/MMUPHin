rm(list = ls())
library(tidyverse)
load("../ibd_paper/data/phyloseq/genera_adj.RData")
source("../ibd_paper/functions/phyloseq_helpers.R")
source("../ibd_paper/functions/helpers.R")

mat_otu_adj <- otu_table2(phylo_genera_adj)
metadata <- sample_data2(phylo_genera_adj)
dir_output.tmp <- "debugging/11_05_Maaslin2/"
dir.create(dir_output.tmp, showWarnings = FALSE, recursive = TRUE)
metadata_test <- metadata %>%
  dplyr::mutate(B.cat_fill = fill_na(B.cat),
                L.cat_fill = fill_na(L.cat),
                E.cat_fill = fill_na(E.cat),
                race_fill = fill_na(race),
                gender_fill = fill_na(gender),
                age.cat_fill =
                  dplyr::case_when(age < 18 ~ "child",
                                   !is.na(age) ~ "adult",
                                   TRUE ~ NA_character_) %>% fill_na,
                age_at_diagnosis.cat_fill =
                  dplyr::case_when(age_at_diagnosis < 18 ~ "A1",
                                   age_at_diagnosis < 40 ~ "A2",
                                   !is.na(age_at_diagnosis) ~ "A3",
                                   TRUE ~ NA_character_) %>% fill_na,
                antibiotics_fill = fill_na(antibiotics),
                immunosuppressants_fill = fill_na(immunosuppressants),
                steroids_fill = fill_na(steroids),
                mesalamine_5ASA_fill = fill_na(mesalamine_5ASA),
                subject_accession_lgtdn = ifelse(dataset_name %in% c("Herfarth_CCFA_Microbiome_3B_combined",
                                                                    "HMP2",
                                                                    "Jansson_Lamendella_Crohns",
                                                                    "LSS-PRISM",
                                                                    "PROTECT"),
                                                 subject_accession,
                                                 "not_longitudinal"))
metadata_test <- metadata_test %>%
  dplyr::mutate(study_site = paste0(dataset_name, sample_type),
                study_site_disease = paste0(dataset_name, sample_type, disease),
                IBD = disease %>%
                  dplyr::recode("CD" = "IBD",
                                "UC" = "IBD",
                                "control" = "control",
                                .default = NA_character_))

test_tmp <- "IBD_vs_control"
metadata_tmp <- metadata_test %>%
  dplyr::filter(!is.na(IBD)) %>%
  dplyr::mutate(test_variable = IBD %>%
                  dplyr::recode("IBD" = "1",
                                "control" = "0",
                                .missing = NA_character_))
metadata_filter_tmp <-
  metadata_tmp %>%
  dplyr::group_by(study_site) %>%
  dplyr::summarise(available = dplyr::n_distinct(test_variable, na.rm = TRUE))
metadata_tmp <- metadata_tmp %>%
  dplyr::left_join(metadata_filter_tmp, by = "study_site") %>%
  dplyr::filter(available == 2)
rownames(metadata_tmp) <- metadata_tmp$sample_accession_16S
test <- MMUPHin::lm.meta(
  feature.count = mat_otu_adj[, metadata_tmp$sample_accession_16S],
  batch = "study_site",
  exposure = "test_variable",
  covariates = c("gender_fill", "race_fill", "age.cat_fill"),
  covariates.random = "subject_accession_lgtdn",
  data = metadata_tmp,
  directory = dir_output.tmp
)


metadata_tmp <- metadata_tmp %>%
  dplyr::filter(study_site == "CS-PRISMbiopsy")
rownames(metadata_tmp) <- metadata_tmp$sample_accession_16S

test <- Maaslin2::Maaslin2(input_data = mat_otu_adj[, rownames(metadata_tmp)],
                           input_metadata = metadata_tmp,
                           output = dir_output.tmp,
                           min_abundance = 0,
                           min_prevalence = 0,
                           normalization = "TSS",
                           transform = "AST",
                           analysis_method = "LM",
                           random_effects = "subject_accession",
                           fixed_effects = c("test_variable", "age.cat_fill"),
                           standardize = FALSE)


otu_test_tmp <- matrix(c(runif(4),
                         c(0, 0, 0.1, 1),
                         c(0, 0, 0, 1),
                         c(0, 0, 0, 0),
                         c(0, 0, 0, 0)),
                       ncol = 4,
                       byrow = TRUE)
meta_test_tmp <- data.frame(x1 = c("a", "a", "a", "b"),
                            x2 = c("a", "a", "b", "b"))


colnames(otu_test_tmp) <- rownames(meta_test_tmp) <- 1:4
test <- Maaslin2::Maaslin2(otu_test_tmp,
                           meta_test_tmp,
                           output = dir_output.tmp,
                           min_abundance = 0,
                           min_prevalence = 0,
                           normalization = "TSS",
                           transform = "AST",
                           analysis_method = "LM",
                           fixed_effects = c("x1"),
                           random_effects = "x2",
                           standardize = FALSE,plot_heatmap = FALSE, plot_scatter = FALSE)

# test
set.seed(1)
yi.test <- rnorm(5)
test1 <- metafor::rma.uni(yi = yi.test, sei = rep(1, 5))
test2 <- metafor::rma.uni(yi = yi.test, sei = rep(1, 5), mods = c(1, 1, 2, 2, 2))
