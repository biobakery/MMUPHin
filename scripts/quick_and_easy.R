library(tidyverse)
library(phyloseq)
l_biom <- c("BIDMC-FMT", "CS-PRISM", "LSS-PRISM",
            "Pouchitis", "RISK", "Jansson_Lamendella_Crohns") %>%
  purrr::map(function(study) {
    tmp <- import_biom(paste0("../ibd_meta_analysis/processed/",
                              study,
                              "/16s/all_samples_taxonomy_closed_reference.biom"))
    return(tmp)
  })
l_taxa <- l_biom %>%
  purrr::map(function(i_biom) {
    taxa_names(i_biom)
  })
l_biomRa <- l_biom %>%
  purrr::map(function(i_biom) i_biom %>%
               transform_sample_counts(function(x) x / sum(x)))
taxa_common <- l_taxa[1:5] %>% Reduce("union", .)

otu_table_all <- l_biom[1:5] %>%
  purrr::map(function(i_biom) {
    mat_tmp <- otu_table(i_biom)@.Data
    mat_na <- matrix(0, nrow = length(taxa_common),
                     ncol = nsamples(i_biom))
    dimnames(mat_na) <- list(taxa_common, sample_names(i_biom))
    mat_na[rownames(mat_tmp), ] <- mat_tmp
    return(mat_na)
  }) %>%
  Reduce("cbind", .)
tax_table_all <- taxa_common %>%
  sapply(function(taxa) {
    for(i in 1:5) {
      i_table <- tax_table(l_biom[[i]])@.Data
      if(taxa %in% rownames(i_table)) return(i_table[taxa, ])
    }
  }) %>% t
dimnames(tax_table_all) <- list(taxa_common, paste0("Rank", 1:7))

l_metadata <-  c("BIDMC-FMT", "CS-PRISM", "LSS-PRISM", "RISK") %>%
  purrr::map(function(study) {
    read_tsv(paste0("../ibd_meta_analysis/raw/",
                    study,
                    "/metadata/",
                    study,
                    "_common.txt"))
  })
metadata1 <- l_metadata %>% bind_rows()
load("../IBD_structure/data/clean/phylo_all_genus_with_tree.RData")
metadata2 <- sample_data(phylo_all_genus) %>%
  data.frame %>%
  filter(collection == "MSH")
mapping <- read_tsv("../ibd_meta_analysis/data/metadata_raivo/sample2project_common.txt")
metadata2 <- metadata2 %>%
  mutate(Project = "Pouchitis",
         OriginalID = DirkSampleID,
         DonorID = subject,
         Location = biopsy_location %>%
           recode("Terminalileum" = "Terminal Ileum",
                  .default = NA_character_),
         Diagnosis = Diagnosis %>%
           recode("control" = "Control"),
         Age = age
  ) %>%
  dplyr::select(Project, DonorID, OriginalID, Location, Diagnosis, Gender, Age)

metadata <- metadata1 %>%
  bind_rows(metadata2) %>%
  left_join(mapping, by = "OriginalID") %>%
  filter(!duplicated(GID), !is.na(GID)) %>%
  as.data.frame()
rownames(metadata) <- metadata$GID
samples <- intersect(rownames(metadata), colnames(otu_table_all))
qe_biom <- phyloseq(
  otu_table(otu_table_all[, samples], taxa_are_rows = TRUE),
  sample_data(metadata[samples,]),
  tax_table(tax_table_all)
)
save(qe_biom, file = "results/quick_and_easy/qe_biom.RData")
metadata_qe <- sample_data(qe_biom) %>% data.frame()
RISK_GID <- metadata_qe$GID[metadata_qe$Project.x == "RISK"]
RISK_remove <- RISK_GID[sample.int(n = length(RISK_GID),
                                   size = 682,
                                   replace = FALSE)]
metadata_qe <- metadata_qe %>%
  filter(!(GID %in% RISK_remove),
         !is.na(Location.x),
         !is.na(Diagnosis.x),
         Diagnosis.x %in% c("CD", "UC", "Control")) %>%
  mutate(Location = Location.x %>%
           recode("L1" = "Terminal Ileum",
                  "L2" = "Terminal Ileum",
                  "L3" = "Terminal Ileum",
                  "L1 + L4" = "Terminal Ileum"),
         AgeG = c("<=18", ">18")[(Age.x > 18) + 1],
         Diagnosis = Diagnosis.x,
         Gender = Gender.x,
         Project = Project.x)
rownames(metadata_qe) <- metadata_qe$GID
qe_biom <- phyloseq(
  otu_table(otu_table_all[, rownames(metadata_qe)], taxa_are_rows = TRUE),
  sample_data(metadata_qe),
  tax_table(tax_table_all)
)
mat_otu_qe <- otu_table(qe_biom)@.Data
samples_remain <- colnames(mat_otu_qe)[!apply(mat_otu_qe == 0, 2, all) ]
mat_otu_qe <- mat_otu_qe[, samples_remain]
metadata_qe <- metadata_qe[samples_remain, ]
qe_biom <- phyloseq(
  otu_table(mat_otu_qe, taxa_are_rows = TRUE),
  sample_data(metadata_qe),
  tax_table(tax_table_all)
)
save(qe_biom, file = "results/quick_and_easy/qe_biom.RData")
qe_biom_genus <- qe_biom %>%
  tax_glom(taxrank = "Rank6")
save(qe_biom_genus, file = "results/quick_and_easy/qe_biom.RData")
mat_otu_qe <- otu_table(qe_biom_genus)@.Data
save(qe_biom_genus, file = "data/clean/qe_genus.Rdata")
qe_biom_ra <- qe_biom_genus %>% transform_sample_counts(function(x) x / sum(x))
library(vegan)
vegan::adonis(distance(qe_biom_ra, "bray") ~
                Project + Diagnosis + Location,
              data = metadata_qe,
              permutations = 2)
plot_ordination(qe_biom_ra, ordinate(qe_biom_ra, method = "NMDS", distance = "bray"),
                color = "Project") +
  theme_bw()
mat_otu_adj <- MMUPHin::adjust.batch(
  feature.count = otu_table(qe_biom_genus)@.Data,
  batch = metadata_qe$Project,
  formula.adj = ~ Diagnosis + Location + AgeG,
  data.adj = metadata_qe,
  zero.inflation = TRUE,
  diagnostics = TRUE
)
vegan::adonis(distance(mat_otu_adj %>%
                         apply(2, function(x) x / sum(x)) %>%
                         otu_table(taxa_are_rows = TRUE), "bray") ~
                Project,
              data = metadata_qe,
              permutations = 2)


plot_ordination(qe_biom_ra, ordinate(otu_table(mat_otu_adj,
                                               taxa_are_rows = TRUE) %>%
                                       transform_sample_counts(function(x) x / sum(x)),
                                     method = "NMDS", distance = "bray"),
                color = "Project") +
  theme_bw()
results.tmp <- results.tmp %>%
  left_join(
    tax_table_all %>%
      data.frame() %>%
      rownames_to_column("feature"),
    by = "feature"
  )
write_csv(results.tmp, "~/Downloads/tmp.csv")
