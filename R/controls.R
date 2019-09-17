control_adjust_batch <- list(
  # logical. Whether or not a zero-inflated model should be
  # run. Default to TRUE (zero-inflated model). If set to FALSE then ComBat
  # (with parametric adjustment) will be performed.
  zero_inflation = TRUE,
  # numeric. Pseudo count to add feature_abd before
  # log-transformation. Automatically set to half of minimal non-zero value if
  # not specified.
  pseudo_count = NULL,
  # logical. Whether or not diagnostic plots are generated.
  diagnostics = "adjust_batch_diagnostics.pdf",
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
  # logical. Whether or not forest plots are generated.
  forest_plots = TRUE,
  # character. Output directory for Maaslin2 output and forest plots.
  output = "./MMUPHin_lm.meta/",
  rma_conv = 1e-4,
  rma_maxit = 1000,
  # logical. Whether or not verbose modeling information is printed.
  verbose = TRUE
)
