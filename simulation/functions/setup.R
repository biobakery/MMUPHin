lib.toLoad <-
  c("magrittr", # for easy operators
    "ggplot2",
    "phyloseq", # for phyloseq objects，
    "foreach", "doParallel"# for paralleling
  )

setup <- function() {

  # Load global packages
  for(i.lib in lib.toLoad) library(i.lib, character.only = TRUE)

  # Plotting options
  theme_set(theme_bw())
}
