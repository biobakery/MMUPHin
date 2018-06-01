rm(list = ls())
library(sparseDOSSA)
library(tidyverse)
library(vegan)
library(phyloseq)
library(MMUPHin)


# simulation setup -------------------------------------------------------
nSubjects <- c(50, 100, 200, 400)
nPerSubject <- 1
# nSamples<- round(nSubjects*nPerSubject) # Number of Samples
nMicrobes <- 50
spikeMicrobes <- 0.2
nMetadata <- 2
noZeroInflate <- FALSE
spikeCount <- 1
Metadatafrozenidx <- 1
effectSize <- c(0.1, 1, 10, 100)
rep <- 1:10

tb_simSetup <- tidyr::crossing(
  nSubjects,
  nPerSubject,
  nMicrobes,
  spikeMicrobes,
  nMetadata,
  noZeroInflate,
  spikeCount,
  Metadatafrozenidx,
  effectSize,
  rep
) %>%
  dplyr::mutate(
    nSamples = nSubjects * nPerSubject,
    `i setup` = 1:nrow(.)
  )


# metadata ----------------------------------------------------------------
l_metadata_null <- tb_simSetup$`i setup` %>%
  purrr::map(function(i_iter) {
    rbind(rnorm(n = tb_simSetup[i_iter, ]$nSamples),
          rnorm(n = tb_simSetup[i_iter, ]$nSamples))
  })
l_tb_metadata <- l_metadata_null %>%
  purrr::map(function(metadata_null) {
    metadata_null %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Sample = paste0("Sample", 1:nrow(.)))
  })


# simulation --------------------------------------------------------------
l_simResults <- (1:nrow(tb_simSetup)) %>%
  purrr::map(function(i_iter) {
    i_simSetup <- tb_simSetup[i_iter, ]
    otu_count <- sparseDOSSA::sparseDOSSA(number_features = i_simSetup$nMicrobes,
                                          number_samples = i_simSetup$nSamples,
                                          UserMetadata = l_metadata_null[[i_iter]],
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
save(l_metadata_null, file = "results/unitTest_sparseDOSSA_sampleSize/metadata.RData")
save(l_simResults, file = "results/unitTest_sparseDOSSA_sampleSize/simResults.RData")

# R2 -------------------------------
tb_R2 <- l_simResults %>%
  purrr::map_df(function(i_simResult) {
    otu_count <- i_simResult$otu_count
    i_iter <- i_simResult$i_iter
    i_simSetup <- tb_simSetup[i_iter, ]
    tb_metadata <- l_tb_metadata[[i_iter]]

    otu_ra <- otu_count %>%
      apply(2, function(x) x / sum(x))
    distance <- otu_ra %>%
      phyloseq::otu_table(taxa_are_rows = TRUE) %>%
      phyloseq::phyloseq() %>%
      phyloseq::distance(method = "bray")

    c("V1", "V2") %>%
      purrr::map_df(function(Variable) {
        tb_selected <- tibble::tibble(V = pull(tb_metadata, Variable))
        adonis_fit <- vegan::adonis(distance ~ V, data = tb_selected, permutations = 2)$aov.tab
        adonis_fit %>%
          dplyr::slice(1) %>%
          dplyr::mutate(Variable = Variable)
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(`i setup` = i_iter)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(tb_simSetup, by = 'i setup')
tb_R2 <- tb_R2 %>%
  mutate(`Sample Size` = nSamples %>% factor(levels = c(50, 100, 200, 400)))
plot <- tb_R2 %>%
  ggplot(aes(x = `Sample Size`, y = R2)) +
  geom_boxplot() +
  facet_grid(Variable ~ effectSize, scales = "free_y") +
  theme_bw()
ggsave(file = "results/unitTest_sparseDOSSA_sampleSize/R2.pdf",
       plot,
       width = 10,
       height = 8)

# is this the property of sparseDOSSA, or PERMANOVA? ----------------------
nSubjects <- c(50, 100, 200, 400)
nPerSubject <- 1
effectSize <- c(0.1, 0.5, 1, 5, 10)
rep <- 1:10

tb_simSetup <- tidyr::crossing(
  nSubjects,
  nPerSubject,
  effectSize,
  rep
) %>%
  dplyr::mutate(
    nSamples = nSubjects * nPerSubject,
    `i setup` = 1:nrow(.)
  )


# metadata ----------------------------------------------------------------
l_metadata_null <- tb_simSetup$`i setup` %>%
  purrr::map(function(i_iter) {
    rbind(rnorm(n = tb_simSetup[i_iter, ]$nSamples),
          rnorm(n = tb_simSetup[i_iter, ]$nSamples),
          rnorm(n = tb_simSetup[i_iter, ]$nSamples))
  })
l_tb_metadata <- l_metadata_null %>%
  purrr::map(function(metadata_null) {
    metadata_null %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Sample = paste0("Sample", 1:nrow(.)))
  })

# R2 -------------------------------
tb_R2 <- (1:length(l_tb_metadata)) %>%
  purrr::map_df(function(i_iter) {
    tb_metadata <- l_tb_metadata[[i_iter]]
    i_simSetup <- tb_simSetup[i_iter, ]
    distance <- dist(tb_metadata[, 1:2] %>%
                       mutate(V1 = V1 * i_simSetup$effectSize))

    c("V1", "V3") %>%
      purrr::map_df(function(Variable) {
        tb_selected <- tibble::tibble(V = pull(tb_metadata, Variable))
        adonis_fit <- vegan::adonis(distance ~ V, data = tb_selected, permutations = 2)$aov.tab
        adonis_fit %>%
          dplyr::slice(1) %>%
          dplyr::mutate(Variable = Variable)
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(`i setup` = i_iter)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(tb_simSetup, by = 'i setup')
tb_R2 <- tb_R2 %>%
  mutate(`Sample Size` = nSamples %>% factor(levels = c(50, 100, 200, 400)))
plot <- tb_R2 %>%
  ggplot(aes(x = `Sample Size`, y = R2)) +
  geom_boxplot() +
  facet_grid(Variable ~ effectSize, scales = "free_y") +
  theme_bw()
ggsave(file = "results/unitTest_sparseDOSSA_sampleSize/R2_normalData.pdf",
       plot,
       width = 10,
       height = 8)


