control_adjust_batch <- list(
  # logical. Whether or not a zero-inflated model should be
  # run. Default to TRUE (zero-inflated model). If set to FALSE then ComBat
  # (with parametric adjustment) will be performed.
  zero_inflation = TRUE,
  # numeric. Pseudo count to add feature_abd before
  # log-transformation. Automatically set to half of minimal non-zero value if
  # not specified.
  pseudo_count = NULL,
  # character (or null). Name of diagnostic plot file.
  diagnostic_plot = "adjust_batch_diagnostic.pdf",
  conv = 1e-4,
  maxit = 1000,
  # logical. Whether or not verbose modeling information is printed.
  verbose = TRUE
)

control_lm_meta <- list(
  # character. Normalization parameter for Maaslin2.
  normalization = "TSS",
  # character. Transformation parameter for Maaslin2.
  transform = "AST",
  # character. Analysis method parameter for Maaslin2.
  analysis_method = "LM",
  # character. Method parameter for rma.
  rma_method = "REML",
  rma_conv = 1e-4,
  rma_maxit = 1000,
  # character. Output directory for Maaslin2 output and forest plots.
  output = "MMUPHin_lm_meta/",
  # character (or null). Suffix of forest plot file.
  forest_plot = "forest.pdf",
  # logical. Whether or not verbose modeling information is printed.
  verbose = TRUE
)

control_discrete_discover <- list(
  # integer. Maximum number of clusters to perform/evaluate on.
  k_max = 10,
  # function. Clustering function.
  cluster_function = fpc::claraCBI,
  # character. Needs to be either "centroid" or "knn" 
  classify_method = "centroid",
  # integer. Number of randomized iterations to run for prediction strength.
  M = 30,
  # intger. Number of nearest neighbors when using method knn
  nnk = 1,
  diagnostic_plot = "discrete_diagnostic.pdf",
  verbose = TRUE
)

control_continuous_discover <- list(
  normalization = "TSS",
  # transformation parameter.
  transform = "AST",
  # pseudo count to add to the count table before log-transformation. 
  pseudo_count = NULL,
  # percentage variance explained cutoff to choose the top PCs.
  var_perc_cutoff = 0.8,
  # correlation cutoff to construct edges for the PC network.
  cos_cutoff = 0.5,
  # function used to perform network community structure discovery.
  cluster_function = igraph::cluster_optimal,
  # output file for clustered PC network.
  network_plot = "clustered_network.pdf",
  # cluster size cutoff (for cluster to be included in the visualized network
  plot_size_cutoff = 2,
  # output file for diagnostic plot.
  diagnostic_plot = "continuous_diagnostic.pdf",
  verbose = TRUE
)
