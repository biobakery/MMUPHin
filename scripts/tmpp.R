# dir.create("tmp")
# load("../IBD_structure/data/clean/phylo_all_genus_with_tree.RData")
# phylo_all_genus <- phylo_all_genus %>%
#   subset_samples(biopsy_location %in% c("Rectum", "Terminalileum", "stool"))
#
# # panel A -----------------------------------------------------------------
# df_meta <- sample_data(phylo_all_genus)
# class(df_meta) <- "data.frame"
# # df_meta$sample_type
# phylo_all_genus_adj <- phylo_all_genus %>%
#   batch_adj(batch = df_meta$collection,
#             adj_model = ~ Diagnosis + biopsy_location)
# dist_orig <- phylo_all_genus %>%
#   normalize(method = "TSS") %>%
#   distance(method = "bray")
# dist_adj <- phylo_all_genus_adj %>%
#   normalize(method = "TSS") %>%
#   distance(method = "bray")
#
# l_fit <- c("collection", "Diagnosis", "antibiotics", "sample_type") %>%
#   lapply(function(var) {
#     model <- as.formula(paste0("~", var))
#     fit1 <- model %>% as.character %>% c('dist_orig', .) %>% paste(collapse = '') %>%
#       as.formula %>% vegan::adonis(data = df_meta, permutations = 2)
#     fit2 <- model %>% as.character %>% c('dist_adj', .) %>% paste(collapse = '') %>%
#       as.formula %>% vegan::adonis(data = df_meta, permutations = 2)
#     data.frame(
#       R2 = c(fit1$aov.tab$R2[1], fit2$aov.tab$R2[1]),
#       variable = var,
#       data = c("original", "adjusted")
#     )
#   })
# df_plot <- l_fit %>% Reduce("rbind", .) %>%
#   mutate(Variable = variable %>%
#            recode_factor(collection = "Study",
#                          Diagnosis = "Disease Subtype",
#                          antibiotics = "Antibiotics",
#                          sample_type = "Body Site"),
#          Adjustment = data %>% factor(levels = c(
#            "original", "adjusted"
#          )))
# pA <- df_plot %>%
#   ggplot(aes(x = Adjustment, y = R2)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~Variable, scales = "free_y", nrow = 1) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 30, hjust = 1))
# # panel B -----------------------------------------------------------------
# load("../IBD_structure/data/clean/phylo_all_genus_with_tree.RData")
# phylo_CD <- phylo_all_genus
# otu_filtered <- phylo_CD %>%
#   transform_sample_counts(function(x) x / sum(x)) %>%
#   filter_taxa(filterfun(kOverA(5, 1e-04)))
# phylo_CD <- prune_taxa(taxa = otu_filtered, x = phylo_CD)
#
# df_meta <- sample_data(phylo_CD)
# df_meta$Disease <- c("IBD", "control")[(df_meta$Diagnosis == "control") + 1] %>%
#   factor(levels = c("control", "IBD"))
# sample_data(phylo_CD) <- df_meta
# class(df_meta) <- "data.frame"
# sample_data(phylo_CD) <- df_meta
# study <- paste0(df_meta$collection, ", ", df_meta$sample_type)
# fit <- meta_fit(physeq = phylo_CD, fit_model = ~ Disease + antibiotics + Gender + age,
#                 study = study, method = "REML",
#                 dir_output = "tmp/")
# fit_contrast <- fit %>%
#   filter(variable == "DiseaseIBD") %>%
#   arrange(order(pval))
# df_tax <- tax_table(phylo_CD)@.Data %>% as.data.frame
# fit_contrast <- df_tax %>%
#   rownames_to_column("feature") %>%
#   right_join(fit_contrast, by = "feature")
# fit_contrast <- fit_contrast %>%
#   filter(!is.na(beta))
# fit_contrast_meta <- fit_contrast %>%
#   filter(!duplicated(feature)) %>%
#   arrange(pval)
# fit_contrast_meta
# df_plot1 <- fit_contrast %>%
#   filter(Rank6 == "g__Roseburia")
# rma_fit_tmp <- rma(yi = df_plot1$beta_all, sei = df_plot1$se_all)
# pB <- ggforest(rma_fit_tmp, bug = "Rosburia")
# df_plot2 <- fit_contrast %>%
#   filter(Rank6 == "g__Veillonella")
# rma_fit_tmp <- rma(yi = df_plot2$beta_all, sei = df_plot2$se_all)
# pC <- ggforest(rma_fit_tmp, bug = "Veillonella")
#
# ggforest = function(x, bug){
#   require("ggplot2")
#   # Function to convert REM results in `rma`-format into a data.frame
#   rma2df = function(x){
#     rbind(
#       data.frame(Study = c("MSH", "PRISM", "RISK biopsy", "RISK stool"), LogFC = x$yi,
#                  CILB=x$yi - 2*sqrt(x$vi),
#                  CIUB=x$yi + 2*sqrt(x$vi),
#                  p = x$pval,
#                  stringsAsFactors = FALSE),
#       data.frame(Study = "RE Model", LogFC = x$b, CILB=x$ci.lb, CIUB=x$ci.ub,
#                  p = x$pval,
#                  stringsAsFactors = FALSE)
#     ) %>% mutate(Study = factor(Study, levels = c("RE Model",
#                                                   "RISK stool",
#                                                   "RISK biopsy",
#                                                   "PRISM",
#                                                   "MSH")),
#                  aggregate = ifelse(Study == "RE Model", "all", "individual"))
#   }
#   remresdf = rma2df(x)
#   remresdf <- transform(remresdf, interval = CIUB - CILB)
#   remresdf <- transform(remresdf, RelConf = 1/interval)
#   p = ggplot(remresdf,
#              aes(LogFC, Study, xmax=CIUB, xmin=CILB)) +
#     # coord_cartesian(xlim=c(-2, 2)) +
#     scale_alpha_discrete(range = c(0.2, 1)) +
#     geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) +
#     geom_errorbarh(alpha=0.5, color="black", height = 0.3) +
#     geom_point(aes(size = RelConf, shape = aggregate)) +
#     geom_point(data = subset(remresdf, Study=="RE Model"), aes(shape = aggregate), size=7) +
#     scale_size(range = c(2, 5), guide=FALSE) +
#     theme_bw() +
#     theme(text = element_text(size=14)) +
#     xlab("Differnetial Abundance (IBD vs. control)") +
#     ggtitle(bug) +
#     scale_shape_manual(values = c(15, 16), guide = F) +
#     theme(plot.title=element_text(face="italic"))
#   return(p)
# }
#
#
# # generate plot -----------------------------------------------------------
# library(cowplot)
# p <- plot_grid(pA, pB, pC, labels = c("A", "B", ""), nrow = 1, rel_widths = c(1.5, 1, 1))
# ggsave(p, file = "tmp/figure.pdf", height = 5, width = 16)
