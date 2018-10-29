rm(list = ls())
library(MMUPHin)
load("debugging/10_12_Dmitry/sample.fffff.RData")
load("debugging/10_12_Dmitry/data.RData")
names(data.list)
feature.count <- data.list$feature.count
data <- data.list$data
# sum(sample.fffff$Cohort %in% c("RISK", "LSS"))
# nrow(sample.fffff)
# taxa.test <- c(runif(157), rep(0, 1666 - 157))
# feature.count <- rbind(feature.count, taxa.test)
debugonce(MMUPHin::lm.meta)
meta.fit <- MMUPHin::lm.meta(feature.count = feature.count,
                             batch = "Cohort",
                             exposure = "disease",
                             #covariates = “Antibiotics”,
                             covariates.random = "subject_new",
                             data = sample.fffff,
                             directory = "debugging/10_12_Dmitry/")
set.seed(1)
nsamples <- 50
test <- rbind(runif(nsamples), runif(nsamples), rep(0, nsamples))
rownames(test) <- paste0("taxa", 1:nrow(test))
metadata <- data.frame(
  fixed.cont = rnorm(nsamples),
  fixed.cat = c("Class1", "Class2")[rbinom(nsamples, size = 1, prob = 0.5) + 1],
  random = rep(1:5, each = 10)
)
rownames(metadata) <- colnames(test) <- paste0("sample", 1:nsamples)
tmp <- Maaslin2::Maaslin2(input_data = test[, , drop = FALSE],
                   input_metadata = metadata,
                   output = "debugging/10_12_Dmitry/",
                   min_abundance = 0, min_prevalence = 0,
                   max_significance = 1, random_effects = "random",
                   fixed_effects = "fixed.disc", standardize = FALSE)
debugonce(MMUPHin:::Maaslin2.wrapper)
tmp <- MMUPHin:::Maaslin2.wrapper(taxa = test[, , drop = FALSE],
                                  metadata = metadata,
                                  directory = "debugging/10_12_Dmitry/",
                                  min_abundance = 0, min_prevalence = 0,
                                  max_significance = 1, covariates.random = "random",
                                  variables = "fixed.cat", standardize = FALSE)
