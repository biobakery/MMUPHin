rm(list = ls())
library(MMUPHin)
library(tidyverse)

load("results/quick_and_easy/qe_biom.RData")
metadata <- phyloseq::sample_data(qe_biom_genus) %>% data.frame
metadata <- metadata %>%
  rownames_to_column("Sample.tmp") %>%
  mutate(SequencingRun_fill = ifelse(is.na(SequencingRun),
                                      Project,
                                      SequencingRun)
  ) %>%
  column_to_rownames("Sample.tmp")
phyloseq::sample_data(qe_biom_genus) <- metadata

qe_biom_adj <- qe_biom_genus
pdf("results/quick_and_easy/batch/per_study.pdf", width = 8, height = 4)
for(project in unique(metadata$Project)) {
  phylo_tmp <- phyloseq::subset_samples(qe_biom_genus, Project == project)
  feature.count_tmp <- phyloseq::otu_table(phylo_tmp)@.Data
  metadata_tmp <- phyloseq::sample_data(phylo_tmp) %>% data.frame
  if(length(unique(metadata_tmp$SequencingRun_fill)) > 1) {
    feature.count_tmp_adj <-
      MMUPHin::adjust.batch(feature.count_tmp,
                            batch = "SequencingRun_fill",
                            data = metadata_tmp)
    phyloseq::otu_table(qe_biom_adj)@.Data[, metadata$Project == project] <-
      feature.count_tmp_adj
  }
}
dev.off()

# Sanity check for adjust.batch
# debugonce(MMUPHin::adjust.batch)
# project <- "RISK"
# phylo_tmp <- phyloseq::subset_samples(qe_biom_genus, Project == project)
# feature.count_tmp <- phyloseq::otu_table(phylo_tmp)@.Data
# metadata_tmp <- phyloseq::sample_data(phylo_tmp) %>% data.frame
# feature.count_tmp_adj <-
#   MMUPHin::adjust.batch(feature.count_tmp,
#                         batch = "SequencingRun_fill",
#                         data = metadata_tmp)
# data.frame(gamma.hat) %>%
#   mutate(Center = apply(gamma.hat, 1, mean, na.rm = TRUE)) %>%
#   gather(key = "Batch", value = "Gamma.hat", X1:X11) %>%
#   ggplot(aes(x = Center, y = Gamma.hat)) +
#     geom_point() +
#     theme_bw()
# log.data.tmp <- log.data
# log.data.tmp[feature.count == 0] <- NA
# df.outlier <- log.data.tmp[apply(abs(gamma.hat), 1, max, na.rm = TRUE) > 5,
#                            ] %>%
#   t %>%
#   data.frame(check.names = FALSE)
# df.outlier %>%
#   cbind(data) %>%
#   gather(key = "taxa",
#          value = "relative abundance",
#          one_of(colnames(df.outlier))) %>%
#   ggplot(aes(x = SequencingRun_fill, y = `relative abundance`)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitter()) +
#   facet_grid(taxa~.) +
#   theme_bw()
# (1:nrow(gamma.hat))[apply(abs(gamma.hat), 1, max, na.rm = TRUE) > 20]
# gamma.hat[553, ]


pdf("results/quick_and_easy/batch/ordination_per_study.pdf",
    width = 8,
    height = 3.5)
for(project in unique(metadata$Project)) {
  phylo_tmp <- phyloseq::subset_samples(qe_biom_genus, Project == project)
  metadata_tmp <- phyloseq::sample_data(phylo_tmp) %>%
    data.frame()
  phylo_tmp_adj <- phyloseq::subset_samples(qe_biom_adj, Project == project)
  dist_tmp <- phyloseq::distance(phylo_tmp %>%
                                   phyloseq::transform_sample_counts(
                                     function(x) x / sum(x)
                                   ),
                                 method = "bray")
  dist_tmp_adj <- phyloseq::distance(phylo_tmp_adj %>%
                                       phyloseq::transform_sample_counts(
                                         function(x) x / sum(x)
                                       ), method = "bray")
  R2_tmp <- R2_tmp_adj <- NA
  if(length(unique(metadata_tmp$SequencingRun_fill)) > 1) {
    fit.pnova <- vegan::adonis(
      dist_tmp ~ SequencingRun_fill,
      data = metadata_tmp,
      permutations = 2
    )
    fit.pnova_adj <- vegan::adonis(
      dist_tmp_adj ~ SequencingRun_fill,
      data = metadata_tmp,
      permutations = 2
    )
    R2_tmp <- fit.pnova$aov.tab[1, 5]
    R2_tmp_adj <- fit.pnova_adj$aov.tab[1, 5]
  }
  p1_tmp <- phyloseq::ordinate(phylo_tmp,
                               method = "MDS",
                               distance = dist_tmp) %>%
    phyloseq::plot_ordination(phylo_tmp,
                              ordination = .,
                              type = "samples",
                              color = "SequencingRun_fill") +
    theme_bw() +
    ggtitle(" Before adjustment R2 = " %>% paste0(project, ., round(R2_tmp*100, digits = 2))) +
    guides(color = FALSE)
  p2_tmp <- phyloseq::ordinate(phylo_tmp_adj,
                               method = "MDS",
                               distance = dist_tmp_adj) %>%
    phyloseq::plot_ordination(phylo_tmp_adj,
                              ordination = .,
                              type = "samples",
                              color = "SequencingRun_fill") +
    theme_bw() +
    ggtitle("After adjustment R2 = " %>% paste0(round(R2_tmp_adj*100, digits = 2))) +
    guides(color = FALSE)
  print(cowplot::plot_grid(p1_tmp, p2_tmp, nrow = 1))
}
dev.off()

phyloseq::otu_table(qe_biom_adj) <-
  MMUPHin::adjust.batch(
    feature.count = phyloseq::otu_table(qe_biom_adj)@.Data,
    batch = "Project",
    covariates = c("Diagnosis", "Location", "AgeG"),
    data = metadata
  ) %>%
  phyloseq::otu_table(taxa_are_rows = TRUE)
dist_all <- phyloseq::distance(qe_biom_genus %>%
                                 phyloseq::transform_sample_counts(
                                   function(x) x / sum(x)
                                 ),
                               method = "bray")
dist_adj <- phyloseq::distance(qe_biom_adj %>%
                                 phyloseq::transform_sample_counts(
                                   function(x) x / sum(x)
                                 ), method = "bray")
fit.pnova <- vegan::adonis(
  dist_all ~ Project,
  data = metadata,
  permutations = 2
)
fit.pnova_adj <- vegan::adonis(
  dist_adj ~ Project,
  data = metadata,
  permutations = 2
)
R2 <- fit.pnova$aov.tab[1, 5]
R2_adj <- fit.pnova_adj$aov.tab[1, 5]
p1 <- phyloseq::ordinate(qe_biom_genus,
                         method = "MDS",
                         distance = dist_all) %>%
  phyloseq::plot_ordination(qe_biom_genus,
                            ordination = .,
                            type = "samples",
                            color = "Project") +
  theme_bw() +
  ggtitle("Before adjustment R2 = " %>% paste0(R2)) +
  guides(color = FALSE)
p2 <- phyloseq::ordinate(qe_biom_adj,
                         method = "MDS",
                         distance = dist_adj) %>%
  phyloseq::plot_ordination(qe_biom_adj,
                            ordination = .,
                            type = "samples",
                            color = "Project") +
  theme_bw() +
  ggtitle("After adjustment R2 = " %>% paste0(R2_adj)) +
  guides(color = FALSE)
pdf("results/quick_and_easy/batch/all_studies.pdf", width = 8, height = 4)
print(cowplot::plot_grid(p1, p2, nrow = 1))
dev.off()


# Meta_analysis -----------------------------------------------------------
qe_tmp <- phyloseq::subset_samples(qe_biom_adj,
                                   Project %in% c("CS-PRISM", "Pouchitis", "RISK") &
                                     Diagnosis %in% c("CD", "Control"))
results <- MMUPHin::lm.meta(phyloseq::otu_table(qe_tmp %>%
                                        phyloseq::transform_sample_counts(
                                          function(x) x / sum(x))
                                      )@.Data,
                            data = phyloseq::sample_data(qe_tmp) %>% data.frame(),
                            batch = "Project",
                            exposure = "Diagnosis",
                            covariates = "Location",
                            directory = "results/quick_and_easy/meta_analysis/")
results <- results %>%
  data.frame() %>%
  rownames_to_column("OTU_ID") %>%
  left_join(phyloseq::tax_table(qe_tmp)@.Data %>%
              data.frame %>%
              rownames_to_column("OTU_ID")) %>%
  arrange(DiagnosisControl.P.val)

# sanity check that Maaslin treats all zero features as NA
# phylo_tmp <- phyloseq::subset_samples(qe_tmp, Project == "CS-PRISM")
# feature.count.tmp <- phylo_tmp %>%
#   phyloseq::transform_sample_counts(function(x) x / sum(x)) %>%
#   phyloseq::otu_table() %>%
#   unclass %>%
#   as.matrix
# metadata.tmp <- phylo_tmp %>%
#   phyloseq::sample_data() %>%
#   data.frame
# feature.count.tmp.NA <- feature.count.tmp
# feature.count.tmp.NA[apply(feature.count.tmp.NA == 0, 1, all), ] <- NA
# result_zero <- MMUPHin:::Maaslin.wrapper(taxa = feature.count.tmp,
#                                          metadata = metadata.tmp,
#                                          variables = "Diagnosis",
#                                          directory = "results/quick_and_easy/test_zero_na/",
#                                          fZeroInflated = FALSE,
#                                          strModelSelection = "none")
# result_NA <- MMUPHin:::Maaslin.wrapper(taxa = feature.count.tmp.NA,
#                                        metadata = metadata.tmp,
#                                        variables = "Diagnosis",
#                                        directory = "results/quick_and_easy/test_zero_na/",
#                                        fZeroInflated = FALSE,
#                                        strModelSelection = "none")
# all(is.na(result_zero$Diagnosis$Coefficient) == is.na(result_NA$Diagnosis$Coefficient))
# plot(result_zero$Diagnosis$Coefficient, result_NA$Diagnosis$Coefficient)
