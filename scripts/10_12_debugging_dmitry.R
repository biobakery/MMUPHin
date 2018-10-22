library(MMUPHin)
load("debugging/10_12_Dmitry/data.RData")
names(data.list)
feature.count <- data.list$feature.count
data <- data.list$data
meta.fit <- MMUPHin::lm.meta(feature.count = feature.count,
                             batch = "Cohort",
                             exposure = "disease",
                             #covariates = “Antibiotics”,
                             covariates.random = "subject_new",
                             data = data,
                             directory = "debugging/10_12_Dmitry/")
