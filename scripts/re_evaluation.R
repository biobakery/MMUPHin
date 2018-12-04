# nCores <- 4
# simparamslabels = c("ZIF", "RE", "metadataType","nSubjects", "nPerSubject",
#                     "nMicrobes", "spikeMicrobes", "nMetadata", "spikeMetadata",
#                     "effectSize", "rep")
#
# # test for no random effect
# simList <- trigger_sparseDOSSA_Simulator(noZeroInflate=FALSE,
#                                          RandomEffect=FALSE,
#                                          metadataType="UVA",
#                                          nSubjects=250,
#                                          nPerSubject=1,
#                                          nMicrobes=100,
#                                          spikeMicrobes=0.2,
#                                          nMetadata=1,
#                                          spikeMetadata=1,
#                                          effectSize=c(1, 2, 5),
#                                          nIterations=25,
#                                          rSeed=1234,
#                                          nCores=4)
#
# simparams <- names(simList)
# no_cores <- nCores
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# df_result_FE <- foreach(i = simparams,
#                         .packages = c("metafor", "phyloseq", "magrittr")) %dopar% {
#
#                           params = strsplit(i, '_')[[1]]
#                           names(params) <- simparamslabels
#
#                           # Extract Relevant Parameters
#                           metadataType <-  as.character(params["metadataType"]) # Type of Metadata
#                           nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
#                           nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
#                           nSamples<-round(nSubjects*nPerSubject) # Number of Samples
#                           nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
#                           spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Percentage of Spiked-in Microbes
#                           nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
#                           spikeMetadata<-as.numeric(params["spikeMetadata"])  # Percentage of Spiked-in Metadata
#                           effectSize<-as.character(params["effectSize"]) # Effect Size
#
#
#                           physeq_tmp <- phyloseq(simList[[i]]$features %>% as.matrix %>% t %>%
#                                                    otu_table(taxa_are_rows = T),
#                                                  simList[[i]]$metadata %>% sample_data())
#
#                           df_result <- meta_fit(physeq_tmp,
#                                                 study = rep(1:5, each = 50),
#                                                 fit_model = ~ Metadata1_TP,
#                                                 method = "FE") %>%
#                             data.frame(metadataType = metadataType,
#                                        spikeMicrobes = spikeMicrobes,
#                                        effectSize = effectSize,
#                                        check.names = T)
#                           return(df_result)
#                         } %>% Reduce("rbind", .)
# df_result_RE <- foreach(i = simparams,
#                         .packages = c("metafor", "phyloseq", "magrittr")) %dopar% {
#
#                           params = strsplit(i, '_')[[1]]
#                           names(params) <- simparamslabels
#
#                           # Extract Relevant Parameters
#                           metadataType <-  as.character(params["metadataType"]) # Type of Metadata
#                           nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
#                           nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
#                           nSamples<-round(nSubjects*nPerSubject) # Number of Samples
#                           nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
#                           spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Percentage of Spiked-in Microbes
#                           nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
#                           spikeMetadata<-as.numeric(params["spikeMetadata"])  # Percentage of Spiked-in Metadata
#                           effectSize<-as.character(params["effectSize"]) # Effect Size
#
#
#                           physeq_tmp <- phyloseq(simList[[i]]$features %>% as.matrix %>% t %>%
#                                                    otu_table(taxa_are_rows = T),
#                                                  simList[[i]]$metadata %>% sample_data())
#
#                           df_result <- meta_fit(physeq_tmp,
#                                                 study = rep(1:5, each = 50),
#                                                 fit_model = ~ Metadata1_TP,
#                                                 method = "REML") %>%
#                             data.frame(metadataType = metadataType,
#                                        spikeMicrobes = spikeMicrobes,
#                                        effectSize = effectSize,
#                                        check.names = T)
#                           return(df_result)
#                         } %>% Reduce("rbind", .)
# stopCluster(cl)
# df_result <- rbind(df_result_FE, df_result_RE) %>%
#   mutate(Alternative = grepl("TP", feature, fixed = T),
#          rej = pval < 0.05)
# save(df_result, file = "result/sim_noRE.RData")
# plot <- df_result %>%
#   group_by(effectSize, method) %>%
#   summarise(FPR = sum(rej & !Alternative, na.rm = T) / sum(!Alternative, na.rm = T),
#             sd_FPR = sqrt(FPR / sum(!Alternative, na.rm = T)),
#             power = sum(rej & Alternative, na.rm = T) / sum(Alternative, na.rm = T),
#             sd_power = sqrt(power / sum(Alternative, na.rm = T))) %>%
#   gather(key = "metric", value = "value", FPR, power) %>%
#   gather(key = "metric2", value = "sd", sd_FPR, sd_power) %>%
#   filter((metric == "FPR" & metric2 == "sd_FPR") | (metric == "power" & metric2 == "sd_power")) %>%
#   ggplot(aes(x = effectSize %>% as.character %>% as.numeric, y = value, color = method)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.3) +
#   facet_grid(metric ~ method, scales = "free_y") +
#   xlab("effect size") +
#   theme_bw()
# ggsave(plot, file = "result/noRE.pdf", width = 6, height = 4)
#
# # test for random effect
# simList <- trigger_sparseDOSSA_Simulator(noZeroInflate=FALSE,
#                                          RandomEffect=TRUE,
#                                          metadataType="UVA",
#                                          nSubjects=5,
#                                          nPerSubject=50,
#                                          nMicrobes=100,
#                                          spikeMicrobes=c(0.2),
#                                          nMetadata=1,
#                                          spikeMetadata=1,
#                                          effectSize=c(1, 2, 5),
#                                          nIterations=25,
#                                          rSeed=1234,
#                                          nCores=4)
#
# simparams <- names(simList)
# no_cores <- nCores
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# df_result_FE <- foreach(i = simparams,
#                         .packages = c("metafor", "phyloseq", "magrittr")) %dopar% {
#
#                           params = strsplit(i, '_')[[1]]
#                           names(params) <- simparamslabels
#
#                           # Extract Relevant Parameters
#                           metadataType <-  as.character(params["metadataType"]) # Type of Metadata
#                           nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
#                           nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
#                           nSamples<-round(nSubjects*nPerSubject) # Number of Samples
#                           nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
#                           spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Percentage of Spiked-in Microbes
#                           nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
#                           spikeMetadata<-as.numeric(params["spikeMetadata"])  # Percentage of Spiked-in Metadata
#                           effectSize<-as.character(params["effectSize"]) # Effect Size
#
#
#                           physeq_tmp <- phyloseq(simList[[i]]$features %>% as.matrix %>% t %>%
#                                                    otu_table(taxa_are_rows = T),
#                                                  simList[[i]]$metadata %>% sample_data())
#
#                           df_result <- meta_fit(physeq_tmp,
#                                                 study = rep(1:5, each = 50),
#                                                 fit_model = ~ Metadata1_TP,
#                                                 method = "FE") %>%
#                             data.frame(metadataType = metadataType,
#                                        spikeMicrobes = spikeMicrobes,
#                                        effectSize = effectSize,
#                                        check.names = T)
#                           return(df_result)
#                         } %>% Reduce("rbind", .)
# df_result_RE <- foreach(i = simparams,
#                         .packages = c("metafor", "phyloseq", "magrittr")) %dopar% {
#
#                           params = strsplit(i, '_')[[1]]
#                           names(params) <- simparamslabels
#
#                           # Extract Relevant Parameters
#                           metadataType <-  as.character(params["metadataType"]) # Type of Metadata
#                           nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
#                           nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
#                           nSamples<-round(nSubjects*nPerSubject) # Number of Samples
#                           nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
#                           spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Percentage of Spiked-in Microbes
#                           nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
#                           spikeMetadata<-as.numeric(params["spikeMetadata"])  # Percentage of Spiked-in Metadata
#                           effectSize<-as.character(params["effectSize"]) # Effect Size
#
#
#                           physeq_tmp <- phyloseq(simList[[i]]$features %>% as.matrix %>% t %>%
#                                                    otu_table(taxa_are_rows = T),
#                                                  simList[[i]]$metadata %>% sample_data())
#
#                           df_result <- meta_fit(physeq_tmp,
#                                                 study = rep(1:5, each = 50),
#                                                 fit_model = ~ Metadata1_TP,
#                                                 method = "REML") %>%
#                             data.frame(metadataType = metadataType,
#                                        spikeMicrobes = spikeMicrobes,
#                                        effectSize = effectSize,
#                                        check.names = T)
#                           return(df_result)
#                         } %>% Reduce("rbind", .)
# stopCluster(cl)
# df_result <- rbind(df_result_FE, df_result_RE) %>%
#   mutate(Alternative = grepl("TP", feature, fixed = T),
#          rej = pval < 0.05)
# save(df_result, file = "result/sim_RE.RData")
# plot <- df_result %>%
#   group_by(effectSize, method) %>%
#   summarise(FPR = sum(rej & !Alternative, na.rm = T) / sum(!Alternative, na.rm = T),
#             sd_FPR = sqrt(FPR / sum(!Alternative, na.rm = T)),
#             power = sum(rej & Alternative, na.rm = T) / sum(Alternative, na.rm = T),
#             sd_power = sqrt(power / sum(Alternative, na.rm = T))) %>%
#   gather(key = "metric", value = "value", FPR, power) %>%
#   gather(key = "metric2", value = "sd", sd_FPR, sd_power) %>%
#   filter((metric == "FPR" & metric2 == "sd_FPR") | (metric == "power" & metric2 == "sd_power")) %>%
#   ggplot(aes(x = effectSize %>% as.character %>% as.numeric, y = value, color = method)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.3) +
#   facet_grid(metric ~ method, scales = "free_y") +
#   xlab("effect size") +
#   theme_bw()
# ggsave(plot, file = "result/RE.pdf", width = 6, height = 4)
# df_PM <- data.frame(
#   effectSize = c(1, 2, 5),
#   R2 = sapply(c(1, 2, 5), function(size) {
#     i <- paste0("ZeroInflate_RandomEffect_UVA_5_50_100_0.2_1_1_", size, "_ 1")
#     physeq_tmp <- phyloseq(simList[[i]]$features %>% as.matrix %>% t %>%
#                              otu_table(taxa_are_rows = T),
#                            simList[[i]]$metadata %>%
#                              data.frame(study = rep(1:5, each = 50)) %>%
#                              sample_data())
#     return(PERMANOVA(physeq_tmp, model = ~ study)$aov.tab[1, 5])
#   })
# )
# df_PM <- data.frame(
#   effectSize = c(1, 2, 5),
#   R2 = sapply(c(1, 2, 5), function(size) {
#     i <- paste0("ZeroInflate_noRandomEffect_UVA_250_1_100_0.2_1_1_", size, "_ 1")
#     physeq_tmp <- phyloseq(simList[[i]]$features %>% as.matrix %>% t %>%
#                              otu_table(taxa_are_rows = T),
#                            simList[[i]]$metadata %>%
#                              data.frame(study = rep(1:5, each = 50)) %>%
#                              sample_data())
#     return(PERMANOVA(physeq_tmp, model = ~ study)$aov.tab[1, 5])
#   })
# )
# plot <- df_PM %>% ggplot(aes(x = effectSize, y = R2*100)) +
#   geom_point() +
#   geom_line() +
#   theme_bw()
# ggsave(plot, file = "result/RE_PM.pdf", width = 6, height = 4)
#
#
# plot <- df_PM %>% ggplot(aes(x = effectSize, y = R2*100)) +
#   geom_point() +
#   geom_line() +
#   theme_bw()
# ggsave(plot, file = "result/noRE_PM.pdf", width = 6, height = 4)
