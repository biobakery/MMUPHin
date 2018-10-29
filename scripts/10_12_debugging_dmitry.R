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
meta.fit.tmp <- lm.meta(feature.count = feature.count,
                        batch = "Cohort",
                        exposure = "disease",
                        #covariates = “Antibiotics”,
                        covariates.random = "subjectID",
                        data = data,
                        method= "REML", # method=“FE” change for fixed effects
                        directory = "debugging/10_12_Dmitry/")
test <- Maaslin2::Maaslin2(
  input_data = feature.count[, sample.fffff$Cohort == "LSS"],
  input_metadata = sample.fffff %>% subset(Cohort == "LSS"),
  output = "debugging/10_12_Dmitry/",
  min_abundance = 0,
  min_prevalence = 0,
  random_effects = "subject_new",
  fixed_effects = "disease",
  standardize = FALSE
)
set.seed(1)
nsamples <- 50
test <- rbind(runif(nsamples), runif(nsamples), runif(nsamples), rep(0, nsamples))
rownames(test) <- paste0("taxa", 1:nrow(test))
metadata <- data.frame(
  fixed.cont = rnorm(nsamples),
  fixed.cat = c("Class1", "Class2", "Class3")[rbinom(nsamples, size = 2, prob = 0.5) + 1],
  random = rep(1:5, each = 10)
)
rownames(metadata) <- colnames(test) <- paste0("sample", 1:nsamples)
tmp <- Maaslin2::Maaslin2(input_data = test[, , drop = FALSE],
                   input_metadata = metadata,
                   output = "debugging/10_12_Dmitry/",
                   min_abundance = 0, min_prevalence = 0.1,
                   max_significance = 1, random_effects = "random",
                   fixed_effects = "fixed.cat", standardize = FALSE,
                   plot_heatmap = FALSE, plot_scatter = FALSE)
# debugonce(MMUPHin:::Maaslin2.wrapper)
tmp <- MMUPHin:::Maaslin2.wrapper(taxa = test[, , drop = FALSE],
                                  metadata = metadata,
                                  directory = "debugging/10_12_Dmitry/",
                                  covariates.random = "random",
                                  variables = c("fixed.cat", "fixed.cont"))
