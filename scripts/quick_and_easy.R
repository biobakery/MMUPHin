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
load("results/quick_and_easy/qe_biom.RData")

feature.count <- otu_table(qe_biom_genus)@.Data
df.meta <- sample_data(qe_biom_genus) %>% data.frame()
feature.count.adj <- feature.count %>%
  MMUPHin::adjust.batch(batch = df.meta$Project,
                        formula.adj = ~ AgeG + Location + Diagnosis + Gender,
                        data = df.meta,
                        zero.inflation = TRUE,
                        pseudo.count = 0.5,
                        diagnostics = TRUE)
feature.count <- feature.count.adj
batch <- df.meta$Project %>% as.factor
lvl.batch <- levels(batch)
n.batch <- length(lvl.batch)

## Transform sample counts
lib.size <- apply(feature.count, 2, sum)
if(all(lib.size <= 1)) {
  warning("Feature table appears to be on the relative abundance scale. ",
          "Setting pseudo.count to be min(feature.count)/2! ")
  pseudo.count <- min(setdiff(feature.count, 0)) / 2
}
log.data <- log(apply(feature.count + pseudo.count,
                      2,
                      function(x) x / sum(x)))
log.data <- asin(sqrt(apply(feature.count, 2, function(x) x / sum(x))))

# Calculate PC for the training sets
if(verbose) message("Performing PCA in individual datasets...")
pca.all <- lapply(lvl.batch, function(lvl)
{
  pc <- log.data[, batch == lvl]
  dat.pca <- prcomp(t(pc))
  return(dat.pca)
})

# Specify the smallest number of PCs to include
ind.var.perc <- lapply(pca.all, function(dat.pca) {
  (cumsum(dat.pca$sdev^2) / sum(dat.pca$sdev^2)) > var.perc.cutoff
})
n.pc.top <- min(sapply(ind.var.perc, length))
ind.var.perc <- sapply(ind.var.perc, function(x) x[1:n.pc.top])
ind.var.perc <- apply(ind.var.perc, 1, all)
if(any(ind.var.perc)) {
  n.pc.top <- min((1:n.pc.top)[ind.var.perc])
  if(verbose) message("Smallest number of PCs passing the specified threshold (",
                      var.perc.cutoff,
                      ") is ",
                      n.pc.top,
                      ".")
} else {
  warning("No number of PCs passed the threshold!\n",
          "Setting number of PCs per dataset to largest possible",
          " value (",
          n.pc.top,
          ").")
}

# First n.pc.top PC loadings for each training dataset
data.loadings <- lapply(pca.all, function(x)
{
  loadings <- x$rotation[,1:n.pc.top]
  return(loadings)
})
mat.data.loading <- Reduce("cbind", data.loadings)
colnames(mat.data.loading) <- unlist(lapply(lvl.batch, function(lvl) {
  paste0(lvl, ", PC", 1:n.pc.top)
}))

# Calculate correlation matrix of loadings across datasets
if(verbose) message("Calculating correlation between PCs across datasets...")
cor.matrix <- abs(t(mat.data.loading) %*% mat.data.loading)
if(!is.null(cor.cutoff)) {
  cor.matrix[cor.matrix < cor.cutoff] <- 0
}

# Create igraph graph object
pc.graph <- igraph::graph_from_adjacency_matrix(cor.matrix,
                                                mode = "undirected",
                                                weighted = TRUE,
                                                diag = FALSE)
# Perform graph community detection
if(verbose) message("Performing network clustering...")
pc.cluster <- igraph::cluster_fast_greedy(pc.graph)
size.communities <- igraph::sizes(pc.cluster)
if(all(size.communities < 2)) {
  warning("No clusters found in the PC network!\n",
          "Consider lowering the value of cor.cutoff.")
  return(NULL)
}
ind.consensus.loading <- size.communities > 1
membership.loading <- igraph::membership(pc.cluster)

# Generate consensus loadings
if(verbose) message("Calculating consensus loadings...")
mat.cons.loading <- sapply((1:max(membership.loading))[ind.consensus.loading],
                           function(i) {
                             loading.tmp <- mat.data.loading[, membership.loading == i]
                             for(j in (2:ncol(loading.tmp))) {
                               if (loading.tmp[, 1] %*% loading.tmp[, j] < 0)
                                 loading.tmp[, j] <- -loading.tmp[, j]
                             }
                             apply(loading.tmp,
                                   1,
                                   mean)
                           })
# Internal validation
mat.cor.vali <- sapply(data.loadings, function(loadings) {
  apply(abs(cor(mat.cons.loading, loadings)), 1, max)
})
qe_biom_adj <- qe_biom_genus
qe_biom_phylum <- qe_biom_adj %>%
  tax_glom(taxrank = "Rank2")
otu_table(qe_biom_adj) <- otu_table(feature.count.adj, taxa_are_rows = TRUE)
dist.tmp <- distance(qe_biom_adj %>% transform_sample_counts(function(x) x / sum(x)), method = "bray")
sample_data(qe_biom_adj) <- df.meta
vegan::adonis2(dist.tmp~`Score 1`, data = df.meta, permutations = 2)
vegan::adonis2(dist.tmp~Project, data = df.meta, permutations = 2)
ordinate.tmp <- ordinate(qe_biom_adj, distance = dist.tmp, method = "MDS")
plot_ordination(qe_biom_adj, ordinate.tmp, color = "Prevotella") +
  facet_grid(.~Project)
mat.cons.loading <- mat.cons.loading[, c(1, 11, 9)]
colnames(mat.cons.loading) <- paste0("Loading ", 1:3)
mat.cons.loading <- mat.cons.loading %>% apply(2, function(x) x / sqrt(sum(x^2)))
df.meta.score <- t(log.data) %*% mat.cons.loading
colnames(df.meta.score) <- paste0("Score", 1:3)
df.meta <- cbind(df.meta, df.meta.score)
mat.tax <- tax_table(qe_biom_adj)@.Data
rownames(mat.tax)[mat.tax[, 6] == "g__Prevotella"]
df.meta$Prevotella <- otu_table(qe_biom_adj %>% transform_sample_counts(function(x)x/sum(x)))@.Data["568118", ]
df.meta %>%
  ggplot(aes(x = Diagnosis, y = `Score 2`)) +
  geom_boxplot() +
  facet_grid(.~Project)
df.meta
tb.loading <- mat.cons.loading %>%
  data.frame(check.names = FALSE) %>%
  cbind(tax_table(qe_biom_genus)@.Data %>% data.frame()) %>%
  as.tibble()
tb.loading <- tb.loading %>%
  arrange(desc(abs(`1`))) %>%
  slice(1:10) %>%
  as.tibble()
tb.loading <- tb.loading %>%
  arrange(desc(abs(`Loading 1`))) %>%
  select(`Loading 1`, `Loading 2`, `Loading 3`, Rank2, Rank5, Rank6)
write_csv(tb.loading, path = "results/quick_and_easy/loadings.csv")


pc.graph.tmp <- pc.graph %>%
  delete.vertices(
    V(pc.graph)[membership(pc.cluster) %in%
                  (names(sizes(pc.cluster))[sizes(pc.cluster) == 1] %>% as.numeric)]
  )
list.membership.tmp <- membership(pc.cluster)[names(V(pc.graph.tmp))]
list.membership.tmp <- unique(list.membership.tmp) %>%
  lapply(function(x) names(list.membership.tmp)[list.membership.tmp == x])
plot.igraph(pc.graph.tmp, mark.groups = list.membership.tmp,
            vertex.size = degree(pc.graph.tmp) / max(degree(pc.graph.tmp)) * 15,
            edge.width = E(pc.graph.tmp)$weight * 10)
pc.cluster.tmp <- igraph::cluster_fast_greedy(pc.graph.tmp)
plot(pc.cluster.tmp, pc.graph.tmp)
# Visualization
pdf(paste0(directory, "network_communities.pdf"),
    width = 10,
    height = 10)
plot(pc.cluster, pc.graph)
dev.off()
if(diagnostics) {
  df.cor.vali <- data.frame(mat.cor.vali)
  colnames(df.cor.vali) <- lvl.batch
  df.cor.vali$community <- factor(1:nrow(df.cor.vali))
  df.plot <- tidyr::gather(df.cor.vali,
                           key = dataset,
                           value = correlation,
                           -community)
  plot.diagnostic <- ggplot2::ggplot(df.plot,
                                     ggplot2::aes(x = community,
                                                  y = correlation)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_point(position = ggplot2::position_jitter()) +
    ggplot2::geom_text(ggplot2::aes(label = dataset), size = 3) +
    ggplot2::theme_bw()
  if(!is.null(cor.cutoff))
    plot.diagnostic <- plot.diagnostic +
    ggplot2::geom_hline(yintercept = cor.cutoff,
                        color = "red")
  ggsave(filename = paste0(directory, "diagnostic.pdf"),
         plot.diagnostic,
         width = min(20, sum(size.communities > 1)),
         height = 4)
}

list(
  consensus.loading = mat.cons.loading,
  membership = lapply((1:length(size.communities))[size.communities > 1],
                      function(i) colnames(mat.data.loading)[igraph::membership(pc.cluster) == i]),
  mat.cor = mat.cor.vali
)

continuous.fit <- MMUPHin::continuous.train(
  feature.count = feature.count,
  batch = batch,
  var.perc.cutoff = 0.8,
  cor.cutoff = 0.5,
  directory = "results/quick_and_easy/"
)
mat.cons.loading <- continuous.fit$consensus.loading
df.cons.loading <-  mat.cons.loading %>%
  abs() %>%
  as.data.frame() %>%
  rownames_to_column("taxa") %>%
  left_join(
    tax_table(qe_biom_genus)@.Data %>%
      as.data.frame() %>%
      rownames_to_column("taxa")
  )
write.csv(df.cons.loading, file = "results/quick_and_easy/loadings.csv")
mat_otu_qe <- otu_table(qe_biom_genus)@.Data
save(qe_biom_genus, file = "data/clean/qe_genus.Rdata")


rm(list = ls())
load("results/quick_and_easy/qe_biom.RData")
qe_biom_ra <- qe_biom_genus %>% transform_sample_counts(function(x) x / sum(x))
mat_otu <- otu_table(qe_biom_genus)@.Data
tb_meta <- sample_data(qe_biom_genus)
class(tb_meta) <- "data.frame"
mat_otu_adj <- mat_otu %>% adjust.batch(batch = tb_meta$Project,
                                        formula.adj = ~ AgeG + Location + Diagnosis,
                                        data = tb_meta,
                                        zero.inflation = TRUE,
                                        pseudo.count = 0.5,
                                        diagnostics = TRUE)
qe_biom_adj <- qe_biom_genus
otu_table(qe_biom_adj) <- otu_table(mat_otu_adj, taxa_are_rows = TRUE)
qe_biom_adj_ra <- qe_biom_adj %>% transform_sample_counts(function(x) x / sum(x))
vegan::adonis(distance(qe_biom_ra, "bray") ~
                Project,
              data = tb_meta,
              permutations = 2)
vegan::adonis(distance(qe_biom_adj_ra, "bray") ~
                Project,
              data = tb_meta,
              permutations = 2)

p1 <- plot_ordination(
  qe_biom_ra,
  ordinate(qe_biom_ra, method = "MDS", distance = "bray"),
  color = "Project"
) +
  theme_bw() +
  ggtitle(expression(paste("Before adjustment, ",
                           R^{2},
                           '=',
                           11.36,
                           "%"))) +
  guides(color=guide_legend(nrow = 3,byrow = TRUE, title = "Study")) +
  theme(legend.position = c(0, 0),
        legend.direction = "horizontal",
        legend.justification = c(0, 0),
        legend.background = element_rect(color = "black"))
p2 <- plot_ordination(
  qe_biom_adj_ra,
  ordinate(qe_biom_adj_ra, method = "MDS", distance = "bray"),
  color = "Project"
) +
  theme_bw() +
  ggtitle(expression(paste("After adjustment, ",
                           R^{2},
                           '=',
                           6.67,
                           "%"))) +
  theme(legend.position = "none")

library(cowplot)
plot_grid(p1, p2, nrow = 1) %>%
  ggsave(file = "~/Dropbox (Huttenhower Lab)/IBD_structure/posters/poster_ISMB_2018/batch.pdf",
         width = 10,
         height = 5)


# boxplots ----------------------------------------------------------------
tb_plot <- otu_table(qe_biom_adj_ra)@.Data %>%
  data.frame %>%
  rownames_to_column("taxa") %>%
  gather(key = GID, value = `Relative Abundance`, -taxa) %>%
  left_join(tb_meta, by = "GID") %>%
  mutate(Diagnosis = Diagnosis %>% factor(levels = c("Control", "CD", "UC")))
p1 <- tb_plot %>%
  filter(Diagnosis != "UC", taxa == "372348",
         Project %in% c("CS-PRISM", "Pouchitis", "RISK")) %>%
  ggplot(aes(x = Project, y = `Relative Abundance`, color = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_bw() +
  xlab("Study") +
  scale_color_manual(values = c("Control" = "blue", "CD" = "red")) +
  guides(color=guide_legend(byrow = TRUE, title = "Diagnosis")) +
  theme(legend.position = c(0, 1),
        legend.direction = "horizontal",
        legend.justification = c(0, 1),
        legend.background = element_rect(color = "black"))
ggsave(p1, file = "~/Dropbox (Huttenhower Lab)/IBD_structure/posters/poster_ISMB_2018/boxplot1.pdf",
       width = 4.5, height = 2.5)
lm1 <- tb_plot %>%
  filter(Diagnosis != "UC", taxa == "372348",
         Project %in% c("CS-PRISM")) %>%
  lm(log(`Relative Abundance` + 1e-5) ~ Diagnosis, data = .) %>%
  summary %>%
  `$`("coef")
lm2 <- tb_plot %>%
  filter(Diagnosis != "UC", taxa == "372348",
         Project %in% c("Pouchitis")) %>%
  lm(log(`Relative Abundance` + 1e-5) ~ Diagnosis, data = .) %>%
  summary %>%
  `$`("coef")
lm3 <- tb_plot %>%
  filter(Diagnosis != "UC", taxa == "372348",
         Project %in% c("RISK")) %>%
  lm(log(`Relative Abundance` + 1e-5) ~ Diagnosis, data = .) %>%
  summary %>%
  `$`("coef")
yis <- c(lm1[2, 1], lm2[2, 1], lm3[2, 1])
seis <- c(lm1[2, 2], lm2[2, 2], lm3[2, 2])
library(metafor)
rma.fit <- rma.uni(yi = yis,
                   sei = seis,
                   slab =  c("CS-PRISM", "Pouchitis", "RISK"))
pdf("~/Dropbox (Huttenhower Lab)/IBD_structure/posters/poster_ISMB_2018/meta1.pdf",
    width = 5.5,
    height = 5.5)
forest(rma.fit, xlab = "Effect Size (CD vs. Control)")
dev.off()

p1 <- tb_plot %>%
  filter(Diagnosis != "UC", taxa == "342427",
         Project %in% c("CS-PRISM", "Pouchitis", "RISK")) %>%
  ggplot(aes(x = Project, y = `Relative Abundance`, color = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  theme_bw() +
  xlab("Study") +
  scale_color_manual(values = c("Control" = "blue", "CD" = "red")) +
  guides(color=guide_legend(byrow = TRUE, title = "Diagnosis")) +
  theme(legend.position = c(0, 1),
        legend.direction = "horizontal",
        legend.justification = c(0, 1),
        legend.background = element_rect(color = "black"))
ggsave(p1, file = "~/Dropbox (Huttenhower Lab)/IBD_structure/posters/poster_ISMB_2018/boxplot2.pdf",
       width = 4.5, height = 2.5)
lm1 <- tb_plot %>%
  filter(Diagnosis != "UC", taxa == "342427",
         Project %in% c("CS-PRISM")) %>%
  lm(log(`Relative Abundance` + 1e-5) ~ Diagnosis, data = .) %>%
  summary %>%
  `$`("coef")
lm2 <- tb_plot %>%
  filter(Diagnosis != "UC", taxa == "342427",
         Project %in% c("Pouchitis")) %>%
  lm(log(`Relative Abundance` + 1e-5) ~ Diagnosis, data = .) %>%
  summary %>%
  `$`("coef")
lm3 <- tb_plot %>%
  filter(Diagnosis != "UC", taxa == "342427",
         Project %in% c("RISK")) %>%
  lm(log(`Relative Abundance` + 1e-5) ~ Diagnosis, data = .) %>%
  summary %>%
  `$`("coef")
yis <- c(lm1[2, 1], lm2[2, 1], lm3[2, 1])
seis <- c(lm1[2, 2], lm2[2, 2], lm3[2, 2])
library(metafor)
rma.fit <- rma.uni(yi = yis,
                   sei = seis,
                   slab =  c("CS-PRISM", "Pouchitis", "RISK"))
pdf("~/Dropbox (Huttenhower Lab)/IBD_structure/posters/poster_ISMB_2018/meta2.pdf",
    width = 5.5,
    height = 5.5)
forest(rma.fit, xlab = "Effect Size (CD vs. Control)")
dev.off()
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
(tb.loading %>%
    arrange(desc(abs(`Loading 1`))) %>%
    slice(1:10) %>%
  ggplot(aes(y = `Loading 1`, x = paste(Rank2 %>% gsub("p__", "", ., fixed = TRUE),
                                        Rank5 %>% gsub("f__", "", ., fixed = TRUE),
                                        Rank6 %>% gsub("g__", "", ., fixed = TRUE), sep = " | "))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  xlab("Genus")
) %>% ggsave(filename = "results/quick_and_easy/loading1.pdf", width = 7, height = 5)
(tb.loading %>%
    arrange(desc(abs(`Loading 2`))) %>%
    slice(1:10) %>%
    ggplot(aes(y = `Loading 2`, x = paste(Rank2, Rank5, Rank6, sep = ","))) +
    geom_bar(stat = "identity") +
    coord_flip()
) %>% ggsave(filename = "results/quick_and_easy/loading2.pdf", width = 10, height = 10)
(tb.loading %>%
    arrange(desc(abs(`Loading 3`))) %>%
    slice(1:10) %>%
    ggplot(aes(y = `Loading 3`, x = paste(Rank2, Rank5, Rank6, sep = ","))) +
    geom_bar(stat = "identity") +
    coord_flip()
) %>% ggsave(filename = "results/quick_and_easy/loading3.pdf", width = 10, height = 10)
df.meta <- df.meta %>%
  cbind(t(log.data) %*% (mat.cons.loading %>%
          apply(2, function(x) x / sqrt(sum(x^2)))))

df.cluster <- lvl.batch %>%
  lapply(function(babtch) {
    physeq_tmp <- qe_biom_adj %>%
      subset_samples(Project == batch)
    dist_tmp <- distance(physeq_tmp, "bray")
    fpc::prediction.strength(dist_tmp,
                             Gmax = 10,
                             M = 30,
                             clustermethod = claraCBI)
  }
  )
df.cluster <- df.cluster %>%
  lapply(function(result_ps) {
    result_ps <- data.frame(
      K = 2:10,
      stat = result_ps$mean.pred[-1],
      se = sapply(result_ps$predcorr[2:10], sd),
      measure = "Pred.Str"
    )
    return(result_ps)
  })
df.cluster <- 1:length(df.cluster) %>%
  lapply(function(i) {
    data.frame(df.cluster[[i]], Study = lvl.batch[i])
  })
df.cluster %>%
  Reduce("rbind", .) %>%
  filter(K <= 7) %>%
  ggplot(aes(x = K, y = stat)) +
  geom_point() +
  geom_line() +
  facet_grid(.~Study) +
  geom_errorbar(aes(ymin = stat - se, ymax = stat + se)) +
  theme_bw() +
  xlab("Number of clusters") +
  ylab("Prediction strength")

