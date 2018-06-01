rm(list = ls())
source("R/sparseDOSSA_helpers.R")
source("R/adjust_batch.R")
source("R/helpers.R")
library(sparseDOSSA)
library(tidyverse)

# simulation setup -------------------------------------------------------
nSubjects <- c(10, 20)
nPerSubject <- c(5, 10, 20)
# nSamples <- round(nSubjects*nPerSubject) # Number of Samples
nMicrobes <- 50
spikeMicrobes <- 0.2
nMetadata <- 3
noZeroInflate <- FALSE
Metadatafrozenidx <- "1,2"
spikeCount <- 2
effectSize <- c(0.5, 5, 50)
rep <- 1:20

tb_simSetup <- tidyr::crossing(
  nSubjects,
  nPerSubject,
  nMicrobes,
  spikeMicrobes,
  noZeroInflate,
  Metadatafrozenidx,
  spikeCount,
  nMetadata,
  effectSize,
  rep
) %>%
  dplyr::mutate(
    nSamples = nSubjects * nPerSubject,
    `i setup` = 1:nrow(.)
  )


# metadata ----------------------------------------------------------------
l_metadata_null <- (1:nrow(tb_simSetup)) %>%
  purrr::map(function(i_iter) {
    i_simSetup <- tb_simSetup[i_iter, ]
    rbind(sample.int(n = i_simSetup$nSubjects,
                     size = i_simSetup$nSamples,
                     replace = TRUE),
          rnorm(n = i_simSetup$nSamples),
          rnorm(n = i_simSetup$nSamples))
  })

# simulation --------------------------------------------------------------
l_simResults <- (1:nrow(tb_simSetup)) %>%
  purrr::map(function(i_iter) {
    i_simSetup <- tb_simSetup[i_iter, ]
    otu_count <- sparseDOSSA::sparseDOSSA(number_features = i_simSetup$nMicrobes,
                                          number_samples = i_simSetup$nSamples,
                                          UserMetadata = l_metadata_null[[i_iter]],
                                          Metadatafrozenidx = i_simSetup$Metadatafrozenidx %>%
                                            stringr::str_split(",") %>%
                                            `[[`(1) %>%
                                            as.integer,
                                          datasetCount = 1,
                                          spikeCount = i_simSetup$spikeCount %>% as.character,
                                          spikeStrength = i_simSetup$effectSize %>% as.character,
                                          noZeroInflate = i_simSetup$noZeroInflate,
                                          percent_spiked = i_simSetup$spikeMicrobes,
                                          write_table = FALSE,
                                          verbose = FALSE) %>%
      extract_sparseDOSSA(add_libsize_var = T) %>%
      `$`("features")
    return(list("otu_count" = otu_count,
                "i_iter" = i_iter))
  })
save(l_metadata_null, file = "results/unitTest_ComBat_batch/metadata.RData")
save(l_simResults, file = "results/unitTest_ComBat_batch/simResults.RData")



# examine sparseDOSSA simulation makes sense ------------------------------
load("results/unitTest_ComBat_batch/metadata.RData")
load("results/unitTest_ComBat_batch/simResults.RData")
l_tb_metadata <- l_metadata_null %>%
  purrr::map(function(metadata_null) {
    metadata_null %>%
      t() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Sample = paste0("Sample", 1:ncol(metadata_null)))
  })
tb_results <- l_simResults %>%
  purrr::map_df(function(i_simResult) {
    otu_count <- i_simResult$otu_count
    i_iter <- i_simResult$i_iter
    i_simSetup <- tb_simSetup[i_iter, ]
    tb_metadata <- l_tb_metadata[[i_iter]]

    otu_ra <- otu_count %>%
      apply(2, function(x) x / sum(x))
    tb_otu <- otu_ra %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample") %>%
      tidyr::gather(key = "feature",
                    value = "relative abundance",
                    -Sample) %>%
      dplyr::left_join(tb_metadata, by = "Sample") %>%
      dplyr::mutate(`i setup` = i_iter)
  }) %>%
  dplyr::bind_rows %>%
  dplyr::left_join(tb_simSetup, by = 'i setup')

tb_summarise <- tb_results %>%
  dplyr::filter(stringr::str_detect(feature, "_TP$")) %>%
  dplyr::group_by(`feature`,
                  `i setup`,
                  `V1`) %>%
  dplyr::summarise(`Group Mean` =
                     `relative abundance` %>%
                     setdiff(0) %>%
                     mean(na.rm = TRUE)) %>%
  dplyr::group_by(`i setup`,
                  `feature`) %>%
  dplyr::summarise(`Overall Mean` = mean(`Group Mean`, na.rm = TRUE),
                   `sd Mean` = sd(`Group Mean`, na.rm = TRUE))
tb_summarise <- tb_summarise %>%
  dplyr::left_join(tb_simSetup, by = "i setup") %>%
  dplyr::mutate(`size per batch` = nPerSubject %>%
                  factor(levels = c(5, 10, 20)),
                `# batches`= nSubjects %>%
                  factor(levels = c(10, 20)))
plot <- tb_summarise %>%
  ggplot2::ggplot(aes(x = `size per batch`, y = `sd Mean`)) +
  ggplot2::geom_boxplot(aes(color = `# batches`)) +
  ggplot2::facet_grid(.~effectSize) +
  theme_bw()

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

    c("V1", "V2", "V3") %>%
      purrr::map_df(function(Variable) {
        tb_selected <- tibble::tibble(V = pull(tb_metadata, Variable))
        if(Variable == "V1") tb_selected <- tb_selected %>% dplyr::mutate(V = factor(V))
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
  mutate(`size per batch` = nPerSubject %>% factor(levels = c(5, 10, 20)),
         `# batches` = nSubjects %>% factor(levels = c(10, 20)))
plot <- tb_R2 %>%
  ggplot(aes(x = `# batches`, y = R2)) +
  geom_boxplot(aes(color = `size per batch`)) +
  facet_grid(Variable ~ effectSize, scales = "free_y") +
  theme_bw()
ggsave(file = "results/unitTest_ComBat_batch/R2.pdf",
       plot,
       width = 10,
       height = 8)
