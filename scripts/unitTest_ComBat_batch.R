rm(list = ls())
source("R/sparseDOSSA_helpers.R")
# source("R/adjust.batch.R")
# source("R/helpers.R")
library(sparseDOSSA)
library(tidyverse)

# simulation setup -------------------------------------------------------
nSamples <- 100
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
  nSamples,
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
    nSubjects = nSamples / nPerSubject,
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
l_results <- l_simResults %>%
  purrr::map(function(i_simResult) {
    otu_count <- i_simResult$otu_count
    i_iter <- i_simResult$i_iter
    i_simSetup <- tb_simSetup[i_iter, ]
    tb_metadata <- l_tb_metadata[[i_iter]]

    list(i_iter = i_iter,
         Original = otu_count,
         ComBat = otu_count %>%
           MMUPHin::adjust.batch(batch = tb_metadata$V1,
                                 formula.adj = ~ V2 + V3,
                                 data.adj = tb_metadata,
                                 zero.inflation = FALSE,
                                 pseudo.count = 0.5,
                                 diagnostics = FALSE,
                                 verbose = FALSE),
         MMUPHin = otu_count %>%
           MMUPHin::adjust.batch(batch = tb_metadata$V1,
                                 formula.adj = ~ V2 + V3,
                                 data.adj = tb_metadata,
                                 zero.inflation = TRUE,
                                 pseudo.count = 0.5,
                                 diagnostics = FALSE,
                                 verbose = FALSE)
    )
  })

tb_R2 <- l_results[141] %>%
  purrr::map_df(function(i_result) {
    i_iter <- i_result$i_iter
    i_simSetup <- tb_simSetup[i_iter, ]
    tb_metadata <- l_tb_metadata[[i_iter]]


    c("Original", "ComBat", "MMUPHin") %>%
      purrr::map_df(function(Normalization) {
        distance <- i_result[[Normalization]] %>%
          apply(2, function(x) x / sum(x)) %>%
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
          dplyr::mutate(Normalization = Normalization)
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(`i setup` = i_iter)
})
tb_results <- tb_R2 %>%
  left_join(tb_simSetup, by = "i setup") %>%
  mutate(Normalization = Normalization %>% factor(levels = c("Original", "ComBat", "MMUPHin")))
tb_results <- tb_results %>%
  mutate(`size per batch` = nPerSubject %>% factor(),
         `# batches` = nSubjects %>% factor())
plot <- tb_results %>%
  ggplot(aes(x = `# batches`, y = R2)) +
  geom_boxplot(aes(color = Normalization)) +
  facet_grid(Variable ~ effectSize, scales = "free_y") +
  scale_color_manual(values = c("Original" = "black",
                                "ComBat" = "red",
                                "MMUPHin" = "blue")) +
  theme_bw()
ggsave(plot, file = "results/unitTest_ComBat_batch/R2.pdf", width = 10, height = 8)

tb_ra <- l_results[141] %>%
  purrr::map_df(function(i_result) {
    otu_count <- i_result$otu_count
    i_iter <- i_result$i_iter
    i_simSetup <- tb_simSetup[i_iter, ]
    tb_metadata <- l_tb_metadata[[i_iter]]

    c("Original", "ComBat", "MMUPHin") %>%
      purrr::map_df(function(Normalization) {
        i_result[[Normalization]] %>%
          apply(2, function(x) x / sum(x)) %>%
          t() %>%
          as.data.frame() %>%
          tibble::rownames_to_column("Sample") %>%
          tidyr::gather(key = "feature",
                        value = "relative abundance",
                        -Sample) %>%
          dplyr::left_join(tb_metadata, by = "Sample") %>%
          dplyr::mutate(Normalization = Normalization)
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(`i setup` = i_iter)
  })
tb_results <- tb_ra %>%
  left_join(tb_simSetup, by = "i setup") %>%
  # dplyr::filter(feature %>% stringr::str_detect("_TP")) %>%
  dplyr::group_by(feature, `i setup`) %>%
  dplyr::mutate(`Mean Abundance Overall` =
                  `relative abundance`[Normalization == "Original"] %>%
                  # setdiff(0) %>%
                  mean) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(feature, `effectSize`, V1, `i setup`,
                  rep, Normalization, `Mean Abundance Overall`,
                  nSubjects) %>%
  dplyr::summarise(`Mean Abundance` = `relative abundance` %>%
                     # setdiff(0) %>%
                     mean) %>%
  ungroup()
plot <- tb_results %>%
  dplyr::mutate(difference = abs(`Mean Abundance` - `Mean Abundance Overall`)) %>%
  dplyr::filter(Normalization != "Original") %>%
  ggplot2::ggplot(aes(x = `Mean Abundance Overall`, y = `difference`)) +
  ggplot2::geom_point(aes(color = Normalization)) +
  # ggplot2::geom_line(aes(group = paste0(feature, V1, `i Iteration`))) +
  ggplot2::geom_abline(slope = 1, intercept = 0) +
  ggplot2::scale_color_manual(values = c("ComBat" = "red",
                                         "MMUPHin" = "blue")) +
  ggplot2::facet_grid(nSubjects~effectSize) +
  ggplot2::theme_bw()
ggsave(plot, file = "results/unitTest_ComBat_batch/meanDiff.pdf", width = 10, height = 8)
