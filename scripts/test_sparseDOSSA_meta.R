rm(list = ls())
library(sparseDOSSA)
library(tidyverse)
library(vegan)
library(phyloseq)
library(MMUPHin)

nSubjects <- 100
nPerSubject <- 1
nSamples<- round(nSubjects*nPerSubject) # Number of Samples
nMicrobes <- 200
spikeMicrobes <- 0.2
nMetadata <- 3
noZeroInflate <- TRUE
spikeCount <- 2
Metadatafrozenidx <- 1:2

c(1, 10, 100, 1000, 10000) %>%
  map(function(effectSize) {

  })

for (effectSize inc(1, 10, 100, 1000, 10000)) {
  otu_count <- sparseDOSSA::sparseDOSSA(number_features = nMicrobes,
                                   number_samples = nSamples,
                                   UserMetadata = metadata_null,
                                   Metadatafrozenidx = Metadatafrozenidx,
                                   datasetCount = 1,
                                   spikeCount = spikeCount %>% as.character,
                                   spikeStrength = effectSize %>% as.character,
                                   noZeroInflate = noZeroInflate,
                                   percent_spiked = spikeMicrobes) %>%
    extract_sparseDOSSA(add_libsize_var = T) %>%
    `$`("features")
  perc_zero <- mean(otu_count == 0)
  distance <- otu %>%
    otu_table(taxa_are_rows = T) %>%
    phyloseq %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    phyloseq::distance(method = "bray")
  adonis_fit <- adonis(distance ~ X1 + X2 + X3,
                       data = metadata_null %>% t %>% data.frame,
                       permutations = 2)

}


metadata_null <- rbind(rbinom(n = nSamples, size = 1, prob = 0.5),
                       rbinom(n = nSamples, size = 1, prob = 0.5),
                       rbinom(n = nSamples, size = 1, prob = 0.5))
metadata1 <- metadata_null
metadata2 <- metadata_null
metadata1[1, ] <- metadata1[1, ] * 10
metadata2[2, ] <- metadata2[2, ] * 10
# metadata2 <- metadata_null
# metadata1 <- cbind(rbinom(n = nSamples, prob = 0.5, size = 1),
#                    rbinom(n = nSamples, prob = 0.5, size = 1) * 5.5) %>% t
# metadata2 <- cbind(rbinom(n = nSamples, prob = 0.5, size = 1) * 5.5,
#                    rbinom(n = nSamples, prob = 0.5, size = 1) * 0.5) %>% t
# metadata1 <- cbind(rnorm(n = nSamples),
#                    rnorm(n = nSamples) * 10) %>% t
# metadata2 <- cbind(rbinom(n = nSamples, prob = 0.5, size = 1) * 5.5,
#                    rbinom(n = nSamples, prob = 0.5, size = 1) * 0.5) %>% t


otu1 <- sparseDOSSA::sparseDOSSA(number_features = nMicrobes,
                                 number_samples = nSamples,
                                 UserMetadata = metadata1,
                                 Metadatafrozenidx = 1:2,
                                 datasetCount = 1,
                                 spikeCount = spikeCount %>% as.character,
                                 spikeStrength = effectSize %>% as.character,
                                 noZeroInflate = noZeroInflate,
                                 percent_spiked = spikeMicrobes) %>%
  extract_sparseDOSSA(add_libsize_var = T) %>%
  `$`("features")
otu1_ra <- otu1 %>% apply(2, function(x) x / sum(x))
median <- otu1_ra %>% apply(1, median)
sd <- otu1_ra %>% apply(1, sd)
zero <- apply(otu1 == 0, 1, mean)
plot(median, sd)
plot(median, zero)
otu2 <- sparseDOSSA::sparseDOSSA(number_features = nMicrobes,
                                 number_samples = nSamples,
                                 UserMetadata = metadata2,
                                 Metadatafrozenidx = 1,
                                 datasetCount = 1,
                                 spikeCount = spikeCount %>% as.character,
                                 spikeStrength = effectSize %>% as.character,
                                 noZeroInflate=noZeroInflate,
                                 percent_spiked=0.2) %>%
  extract_sparseDOSSA(add_libsize_var = T) %>%
  `$`("features")

library(phyloseq)
library(vegan)
distance1 <- otu1 %>%
  otu_table(taxa_are_rows = T) %>%
  phyloseq %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  phyloseq::distance(method = "bray")
adonis(distance1 ~ X1 + X2 + X3,
       data = metadata_null %>% t %>% data.frame,
       permutations = 2)
distance2 <- otu2 %>%
  otu_table(taxa_are_rows = T) %>%
  phyloseq %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  phyloseq::distance(method = "bray")
adonis(distance2 ~ X1 + X2, data = metadata2 %>% t %>% data.frame, permutations = 2)
if (is.null(DD) | inherits(DD, "try-error")) {
  tryAgain = TRUE
  infiniteloopcounter = infiniteloopcounter + 1
} else {
  tryAgain = FALSE
}
}
if (infiniteloopcounter >= 5) {
  stop("Consistent error found during simulation. Need to investigate cause.")
}

# Gather sparseDOSSA outputs
sparsedossa_results <- as.data.frame(DD$OTU_count)
rownames(sparsedossa_results)<-sparsedossa_results$X1
sparsedossa_results<-sparsedossa_results[-1,-1]
colnames(sparsedossa_results)<-paste('Sample', 1:ncol(sparsedossa_results), sep='')
data<-as.matrix(sparsedossa_results[-c((nMetadata+1):(2*nMicrobes+nMetadata)),])
data<-data.matrix(data)
class(data) <- "numeric"
truth<-c(unlist(DD$truth))
truth<-truth[!stri_detect_fixed(truth,":")]
truth<-truth[(5+nMetadata):length(truth)]
truth<-as.data.frame(truth)
significant_features<-as.vector(truth[seq(1, (as.numeric(spikeCount)+1)*(nMicrobes*spikeMicrobes), (as.numeric(spikeCount)+1)),])
