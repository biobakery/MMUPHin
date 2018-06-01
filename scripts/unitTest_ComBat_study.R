rm(list = ls())
source("R/sparseDOSSA_helpers.R")
library(sva)
library(sparseDOSSA)
library(tidyverse)

# simulation setup -------------------------------------------------------
nSubjects <- 200
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
save(l_metadata_null, file = "results/unitTest_ComBat_study/metadata.RData")
save(l_simResults, file = "results/unitTest_ComBat_study/simResults.RData")

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
    otu_count_ComBat <- adjust_batch(dat = otu_count,
                                     batch = tb_metadata$V1,
                                     mod = cbind(tb_metadata$V2),
                                     noZeroInflate = TRUE,
                                     pseudoCount = 0.5)
    otu_count_MMUPHin <- adjust_batch(dat = otu_count,
                                      batch = tb_metadata$V1,
                                      mod = cbind(tb_metadata$V2),
                                      noZeroInflate = FALSE,
                                      pseudoCount = 0.5)
    tb_otu <- otu_ra %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample") %>%
      tidyr::gather(key = "feature",
                    value = "relative abundance",
                    -Sample) %>%
      mutate(Normalization = "Original") %>%
      bind_rows(otu_count_ComBat %>%
                  apply(2, function(x) x / sum(x)) %>%
                  t() %>%
                  as.data.frame() %>%
                  tibble::rownames_to_column("Sample") %>%
                  tidyr::gather(key = "feature",
                                value = "relative abundance",
                                -Sample) %>%
                  mutate(Normalization = "ComBat")) %>%
      bind_rows(otu_count_MMUPHin %>%
                  apply(2, function(x) x / sum(x)) %>%
                  t() %>%
                  as.data.frame() %>%
                  tibble::rownames_to_column("Sample") %>%
                  tidyr::gather(key = "feature",
                                value = "relative abundance",
                                -Sample) %>%
                  mutate(Normalization = "MMUPHin"))
    tb_combined <- tb_otu %>%
      dplyr::left_join(tb_metadata) %>%
      dplyr::bind_cols(
        i_simSetup %>%
          slice(rep(1, nrow(tb_otu)))
      ) %>%
      mutate(`i Iteration` = i_iter)
    # return(tb_combined)
  }) %>%
  dplyr::bind_rows()
save(tb_results, file = "results/unitTest_ComBat_study/tb_combinedAbundance.RData")
df_plots <- tb_results %>%
  dplyr::filter(feature %>% stringr::str_detect("_TP")) %>%
  dplyr::group_by(feature, `effectSize`, `i Iteration`, rep) %>%
  dplyr::mutate(`Mean Abundance Overall` =
                  `relative abundance`[Normalization == "Original"] %>%
                  setdiff(0) %>%
                  mean) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(feature, `effectSize`, V1, `i Iteration`, rep, Normalization, `Mean Abundance Overall`) %>%
  dplyr::summarise(`Mean Abundance` = `relative abundance` %>% setdiff(0) %>% mean) %>%
  ungroup()
plot <- df_plots %>%
  dplyr::mutate(difference = abs(`Mean Abundance` - `Mean Abundance Overall`)) %>%
  dplyr::group_by(feature, `effectSize`, V1, `i Iteration`, rep) %>%
  dplyr::mutate(Original = difference[Normalization == "Original"]) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Normalization != "Original") %>%
  ggplot2::ggplot(aes(x = `Original`, y = `difference`)) +
  ggplot2::geom_point(aes(color = Normalization)) +
  # ggplot2::geom_line(aes(group = paste0(feature, V1, `i Iteration`))) +
  ggplot2::geom_abline(slope = 1, intercept = 0) +
  ggplot2::scale_color_manual(values = c("ComBat" = "red",
                                         "MMUPHin" = "blue")) +
  ggplot2::facet_grid(.~effectSize) +
  ggplot2::theme_bw()

ggplot2::ggsave(file = "results/unitTest_ComBat_study/compare_meanDiff.pdf",
         plot,
         width = 20,
         height = 6)

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
    otu_count_ComBat <- adjust_batch(dat = otu_count,
                                     batch = tb_metadata$V1,
                                     mod = cbind(tb_metadata$V2),
                                     noZeroInflate = TRUE,
                                     pseudoCount = 0.5)
    otu_count_MMUPHin <- adjust_batch(dat = otu_count,
                                      batch = tb_metadata$V1,
                                      mod = cbind(tb_metadata$V2),
                                      noZeroInflate = FALSE,
                                      pseudoCount = 0.5)
    l_distance <- list(
      "Original" = otu_count,
      "ComBat" = otu_count_ComBat,
      "MMUPHin" = otu_count_MMUPHin
    ) %>%
      purrr::map(
        function(otu) {
          otu %>%
            phyloseq::otu_table(taxa_are_rows = TRUE) %>%
            phyloseq::phyloseq() %>%
            phyloseq::distance(method = "bray")
        }
      )

    tb_R2 <- names(l_distance) %>%
      purrr::map_df(
        function(i_normalization) {
          distance <- l_distance[[i_normalization]]
          adonis_fit1 <- vegan::adonis(distance ~ V1,
                                      data = tb_metadata,
                                      permutations = 2)
          adonis_fit2 <- vegan::adonis(distance ~ V2,
                                       data = tb_metadata,
                                       permutations = 2)
          rbind(adonis_fit1$aov.tab[1, ],
                adonis_fit2$aov.tab[1, ]) %>%
            as.data.frame() %>%
            tibble::rownames_to_column("Variable") %>%
            dplyr::mutate(Normalization = i_normalization)
        }
      ) %>%
      dplyr::bind_rows()

    tb_combined <- tb_R2 %>%
      dplyr::bind_cols(
        i_simSetup %>%
          slice(rep(1, nrow(tb_R2)))
      ) %>%
      mutate(`i Iteration` = i_iter)
    # return(tb_combined)
  }) %>%
  dplyr::bind_rows()
plot <- tb_results %>%
  mutate(Normalization = Normalization %>%
           factor(levels = c("Original", "ComBat", "MMUPHin"))) %>%
  filter(Variable %in% c("V1", "V2")) %>%
  ggplot(aes(x = Variable, y = R2)) +
  geom_boxplot(aes(color = Normalization)) +
  facet_grid(.~effectSize) +
  scale_color_manual(values = c("Original" = "black",
                                "ComBat" = "red",
                                "MMUPHin" = "blue")) +
  theme_bw()
ggsave(file = "results/unitTest_ComBat_study/compare_R2.pdf",
       plot,
       width = 10,
       height = 4)
