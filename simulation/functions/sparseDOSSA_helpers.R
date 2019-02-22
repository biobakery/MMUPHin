create_spikein.mt <- function(number_features,
                              percent_spiked,
                              effectSize,
                              seed,
                              same_features = FALSE) {
  set.seed(seed)
  nFeatureSpiked <- floor(number_features * percent_spiked)
  if(nFeatureSpiked == 0)
    stop("No features are spiked in with current configuration!")
  if(same_features) features_spike <- sample.int(n = number_features, size = nFeatureSpiked)
  l_spikein.mt <- lapply(1:length(effectSize), function(i) {
    if(!same_features) features_spike <- sample.int(n = number_features, size = nFeatureSpiked)
    data.frame(feature = features_spike,
               metadata = i,
               strength = effectSize[i],
               stringsAsFactors = FALSE, row.names = NULL)
  })
  return(Reduce("rbind", l_spikein.mt))
}

# this is to make exposure confounded with batch
create_spikein.mt_lm.meta <- function(number_features,
                                      percent_spiked,
                                      effectSize,
                                      seed) {
  set.seed(seed)
  nFeatureSpiked <- floor(number_features * percent_spiked)
  if(nFeatureSpiked == 0)
    stop("No features are spiked in with current configuration!")
  effect_batch <- effectSize[stringr::str_detect(names(effectSize), "batch")]
  effect_exposure <- effectSize[stringr::str_detect(names(effectSize), "exposure")]
  if(!setequal(names(effectSize), c(names(effect_batch),
                                    names(effect_exposure))))
    stop("Something's wrong!")
  l_spikein.mt_batch <- lapply(1:length(effect_batch), function(i) {
    features_spike <- sample.int(n = number_features, size = nFeatureSpiked)
    data.frame(feature = features_spike,
               metadata = i,
               strength = effectSize[i],
               stringsAsFactors = FALSE, row.names = NULL)
  })
  features_batch <- l_spikein.mt_batch %>%
    purrr::map("feature") %>%
    unlist() %>%
    unique()
  l_spikein.mt_exposure <- lapply(1:length(effect_exposure), function(i) {
    features_spike <- sample(x = features_batch, size = nFeatureSpiked, replace = FALSE)
    data.frame(feature = features_spike,
               metadata = i,
               strength = effectSize[i],
               stringsAsFactors = FALSE, row.names = NULL)
  })
  return(Reduce("rbind", c(l_spikein.mt_batch, l_spikein.mt_exposure)))
}

extract_sparseDOSSA <- function(sparseDOSSA_fit) {

  # metadata + feature data
  sparsedossa_results <- as.data.frame(sparseDOSSA_fit$OTU_count)
  rownames(sparsedossa_results) <- sparsedossa_results$X1
  nMetadata <- sum(grepl("Metadata", sparsedossa_results$X1, fixed = TRUE))
  nMicrobes <- sum(grepl("Feature_spike", sparsedossa_results$X1, fixed = TRUE))
  nSamples <- ncol(sparsedossa_results) - 1
  sparsedossa_results <- sparsedossa_results[-1, -1]
  colnames(sparsedossa_results) <- paste('Sample', 1:nSamples, sep='')
  data <- as.matrix(sparsedossa_results[-c((nMetadata+1):(2*nMicrobes+nMetadata)), ])
  data <- data.matrix(data)
  class(data) <- "numeric"

  # Spiked-in features and metadata
  # truth <- c(unlist(sparseDOSSA_fit$truth))
  # truth <- truth[!stringi::stri_detect_fixed(truth,":")]
  # significant_features <- truth[grepl("Feature", truth, fixed = TRUE)]
  # significant_metadata <- truth[-(1+nMetadata)] %>%
  #   stringr::str_subset("Metadata") %>%
  #   stringr::str_replace("_Level_.+", "") %>%
  #   unique

  # Extract Metadata
  # metadata <- as.data.frame(t(data[(1:nMetadata), ]))

  # Rename True Positive Metadata - Same Format at Mcmurdie and Holmes (2014)
  # which.TP <- colnames(metadata) %in% significant_metadata
  # meta_newname <- paste0(colnames(metadata)[which.TP], "_TP")
  # colnames(metadata)[which.TP] <- meta_newname

  # Extract Features
  features <- as.data.frame(t(data[-c(1:nMetadata),]))

  # Rename Features and True Positive Features - Same Format at Mcmurdie and Holmes (2014)
  # wh.TP <- colnames(features) %in% significant_features
  colnames(features) <- paste("Feature", 1:nMicrobes, sep = "")
  # newname <- paste0(colnames(features)[wh.TP], "_TP")
  # colnames(features)[wh.TP] <- newname

  # feature table has rows as features
  features <- t(features)

  # Return as list
  return(list(features=features))
}
