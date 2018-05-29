rm(list = ls())
library(sparseDOSSA)
library(tidyverse)
library(vegan)
library(phyloseq)
library(MMUPHin)


# simulation setup -------------------------------------------------------
nSubjects <- 100
nPerSubject <- 1
nSamples<- round(nSubjects*nPerSubject) # Number of Samples
nMicrobes <- 50
spikeMicrobes <- 0.2
nMetadata <- 2
noZeroInflate <- FALSE
spikeCount <- 1
Metadatafrozenidx <- 1
effectSize <- c(0.1, 0.5, 1, 5, 10, 50, 100)
rep <- 1:20

df_simSetup <- tidyr::crossing(
  nSubjects,
  nPerSubject,
  nSamples,
  nMicrobes,
  spikeMicrobes,
  nMetadata,
  noZeroInflate,
  spikeCount,
  Metadatafrozenidx,
  effectSize,
  rep
)


# metadata ----------------------------------------------------------------
l_metadata_null <- rep %>%
  purrr::map(function(i_rep) {
    rbind(rbinom(n = nSamples, size = 1, prob = 0.5),
          rbinom(n = nSamples, size = 1, prob = 0.5))
  })
l_tb_metadata <- l_metadata_null %>%
  purrr::map(function(metadata_null) {
    metadata_null %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Sample = paste0("Sample", 1:nSamples))
  })


# simulation --------------------------------------------------------------
l_simResults <- (1:nrow(df_simSetup)) %>%
  purrr::map(function(i_iter) {
    i_simSetup <- df_simSetup[i_iter, ]
    otu_count <- sparseDOSSA::sparseDOSSA(number_features = i_simSetup$nMicrobes,
                                          number_samples = i_simSetup$nSamples,
                                          UserMetadata = l_metadata_null[[i_simSetup$rep]],
                                          Metadatafrozenidx = i_simSetup$Metadatafrozenidx,
                                          datasetCount = 1,
                                          spikeCount = i_simSetup$spikeCount %>% as.character,
                                          spikeStrength = i_simSetup$effectSize %>% as.character,
                                          noZeroInflate = i_simSetup$noZeroInflate,
                                          percent_spiked = i_simSetup$spikeMicrobes,
                                          write_table = FALSE,
                                          verbose = FALSE) %>%
      MMUPHin::extract_sparseDOSSA(add_libsize_var = T) %>%
      `$`("features")
    return(list("otu_count" = otu_count,
                "i_iter" = i_iter))
  })
save(l_simResults, file = "results/test_sparseDOSSA/simResults.RData")

# Summarise differential abundance/variance -------------------------------
tb_results <- l_simResults %>%
  purrr::map(function(i_simResult) {
    otu_count <- i_simResult$otu_count
    i_iter <- i_simResult$i_iter
    i_simSetup <- df_simSetup[i_iter, ]
    effectSize <- i_simSetup$effectSize
    rep <- i_simSetup$rep
    tb_metadata <- l_tb_metadata[[rep]]

    otu_ra <- otu_count %>%
      apply(2, function(x) x / sum(x))
    tb_otu <- otu_ra %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample") %>%
      tidyr::gather(key = "feature",
                    value = "relative abundance",
                    -Sample)
    tb_combined <- tb_otu %>%
      dplyr::left_join(tb_metadata) %>%
      dplyr::bind_cols(
        i_simSetup %>%
          slice(rep(1, nrow(tb_otu)))
      ) %>%
      mutate(`i Iteration` = i_iter)
    return(tb_combined)
  }) %>%
  dplyr::bind_rows()
df_plots <- tb_results %>%
  dplyr::filter(feature %>% stringr::str_detect("_TP")) %>%
  dplyr::group_by(feature, `effectSize`, V1, `i Iteration`, rep) %>%
  dplyr::summarise(`Mean Abundance` = `relative abundance` %>% setdiff(0) %>% mean,
                   `sd Abundance` = `relative abundance` %>% setdiff(0) %>% sd,
                   `Mean Abundance with Zero` = `relative abundance` %>% mean,
                   `Percent Zero` = mean(`relative abundance` == 0)) %>%
  dplyr::ungroup() %>%
  tidyr::gather(key = Statistic, value = Value,
                `Mean Abundance`,
                `sd Abundance`,
                `Mean Abundance with Zero`,
                `Percent Zero`) %>%
  tidyr::unite("Group / Statistic", V1, Statistic) %>%
  tidyr::spread(key = `Group / Statistic`, value = Value) %>%
  dplyr::mutate(
    `Mean Difference` = `1_Mean Abundance` - `0_Mean Abundance`,
    `Overall Mean` = (`1_Mean Abundance` + `0_Mean Abundance`) / 2,
    `sd ratio` = `1_sd Abundance` / `0_sd Abundance`,
    `Overall sd` = sqrt((`1_sd Abundance`^2 + `0_sd Abundance`^2)/2),
    `Mean Difference with Zero` = `1_Mean Abundance with Zero` - `0_Mean Abundance with Zero`,
    `Overall Percent Zero` = (`1_Percent Zero` + `0_Percent Zero`) / 2
    )
plot1 <- df_plots %>%
  tidyr::gather(key = `including zero?`, value = `Mean Difference`,
         `Mean Difference`,
         `Mean Difference with Zero`) %>%
  dplyr::mutate(
    `including zero?` = `including zero?` %>%
      recode(
        "Mean Difference" = "no",
        "Mean Difference with Zero" = "yes"
      )
  ) %>%
  ggplot2::ggplot(aes(x = `Overall Mean`, y = abs(`Mean Difference`))) +
  ggplot2::geom_point(aes(color = `including zero?`, alpha = `including zero?`)) +
  ggplot2::scale_color_manual(values = c("no" = "black",
                                         "yes" = "grey")) +
  ggplot2::scale_alpha_manual(values = c("no" = 1,
                                         "yes" = 0.8)) +
  ggplot2::facet_grid(. ~ `effectSize`) +
  theme_bw() +
  theme(legend.position=c(0,1),
        legend.justification=c(0, 1))
plot2 <- df_plots %>%
  ggplot2::ggplot(aes(x = `Overall Mean`, y = log10(`sd ratio`))) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(. ~ `effectSize`) +
  theme_bw()
plot3 <- df_plots %>%
  ggplot2::ggplot(aes(x = `Overall Mean`, y = `Overall Percent Zero`)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(. ~ `effectSize`) +
  theme_bw()
library(cowplot)
ggsave("results/test_sparseDOSSA/effectSize_vs_differentialAbundance.pdf",
       cowplot::plot_grid(plot1, plot2, plot3, ncol = 1),
       width = 20,
       height = 12)

# Summarise PERMANOVA R^2 statistic -------------------------------
tb_results <- l_simResults %>%
  purrr::map(function(i_simResult) {
    otu_count <- i_simResult$otu_count
    i_iter <- i_simResult$i_iter
    i_simSetup <- df_simSetup[i_iter, ]
    effectSize <- i_simSetup$effectSize
    rep <- i_simSetup$rep
    tb_metadata <- l_tb_metadata[[rep]]

    otu_ra <- otu_count %>%
      apply(2, function(x) x / sum(x))
    distance <- otu_ra %>%
      phyloseq::otu_table(taxa_are_rows = T) %>%
      phyloseq::phyloseq() %>%
      phyloseq::distance(method = "bray")
    adonis_fit <- vegan::adonis(distance ~ V1 + V2,
                                data = tb_metadata,
                                permutations = 2)

    adonis_fit$aov.tab %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Variable") %>%
      dplyr::bind_cols(i_simSetup %>%
                         slice(rep(1, nrow(adonis_fit$aov.tab))))
  }) %>%
  dplyr::bind_rows()
df_plots <- tb_results %>%
  filter(Variable %in% c("V1", "V2"))
plot <- df_plots %>%
  ggplot(aes(x = Variable, y = R2)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(.~effectSize) +
  geom_point(position = position_jitter()) +
  theme_bw()
ggsave("results/test_sparseDOSSA/effectSize_vs_R2.pdf",
       plot,
       width = 12,
       height = 4)
