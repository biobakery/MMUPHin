#' #' Simulate microbial read counts data with
#' #' negative binomial distribution
#' #'
#' #' @param read_depth Average read depth
#' #' @param number_feature Number of simulated features
#' #' @param number_samples Number of simulated samples
#' #' @param UserMetadata A data frame with user-specified metadata
#' #' @param spikeStrength Magnitude of association between metadata and microbial features
#' #' @param noZeroInflate Whether zero-inflation should be simulated
#' #' @param prob_spiked Percentage of feature associated with metadata
#' #' @param verbose
#' #'
#' #' @return A phyloseq object with simulated read counts and
#' #' user-specified metadata
#' #'
#' #' @import magrittr phyloseq
#' #' @export
#' #'
#' simulate_NB <- function(read_depth = 10000,
#'                         number_feature,
#'                         number_samples,
#'                         UserMetadata,
#'                         spikeStrength = 1,
#'                         noZeroInflate = FALSE,
#'                         prob_spiked,
#'                         verbose = FALSE)
#' {
#'   if(number_samples != nrow(UserMetadata))
#'     stop("number_samples does not match that in UserMetadata!")
#'   rownames(UserMetadata) <- paste0('Sample', 1:number_samples)
#'   number_metadata <- ncol(UserMetadata)
#'
#'   # simulate null feature count table
#'   # simulate per-feature mean parameter, with log-normal distribution
#'   mean_feature <- exp(rnorm(number_feature, mean = 4.5, sd = 2))
#'   # simulate NB count table
#'   mat_count <- sapply(mean_feature,
#'                       function(mu) rnbinom(n = number_samples,
#'                                            size = 0.3,
#'                                            mu = mu)) %>% t
#'   # simulate zero inflation
#'   if(!noZeroInflate) {
#'     prob_ZI <- log(mean_feature) * (-0.7) + 3.3
#'     prob_ZI <- exp(prob_ZI) / (1 + exp(prob_ZI))
#'
#'     mat_count <- sapply(1:number_feature, function(k) {
#'       abd <- mat_count[k, ]
#'       abd[abd < quantile(abd, prob_ZI[k])] <- 0
#'       return(abd)
#'     }) %>% t
#'   }
#'   mat_count <- mat_count %>% apply(2, function(x) x / sum(x))
#'   # simulate library size
#'   mat_count <- (mat_count %>% t * read_depth *
#'                   exp(rnorm(number_samples, mean = 0, sd = 2.6))) %>%
#'     t
#'
#'   # Design matrix
#'   mat_design <- lapply(1:number_metadata, function(j) {
#'     model.matrix(~., data = data.frame(UserMetadata[, j]))[, -1]
#'   }) %>% Reduce('cbind', .)
#'
#'   # randomized effect size across the number of variables in the design
#'   # matrix and all the features
#'   mat_effectsize <- matrix(runif(ncol(mat_design) * number_feature,
#'                                  min = -2,
#'                                  max = 2),
#'                            nrow = ncol(mat_design),
#'                            ncol = number_feature)
#'
#'   # Deciding which feature/metadata pair to spike in, modified by the spike in
#'   # strength
#'   mat_spike <- matrix(rbinom(n = number_metadata * number_feature,
#'                              size = 1, prob = prob_spiked),
#'                       nrow = number_metadata,
#'                       ncol = number_feature)
#'   mat_spike <- mat_spike * spikeStrength
#'   mat_spike_expand <- lapply(1:number_metadata, function(j) {
#'     rep <- 1
#'     if(is.factor(UserMetadata[, j])) rep <- nlevels(UserMetadata[, j]) - 1
#'     matrix(mat_spike[j, ], nrow = number_feature, ncol = rep) %>% t %>% return
#'   }) %>% Reduce('rbind', .)
#'
#'   # modifier matrix
#'   mat_modify <- t(mat_design %*% (mat_effectsize * mat_spike_expand))
#'
#'   # final otu table
#'   mat_count <- mat_count * exp(mat_modify)
#'   colnames(mat_count) <- rownames(UserMetadata)
#'
#'   return(phyloseq(otu_table(mat_count, taxa_are_rows = T),
#'                   sample_data(UserMetadata)))
#' }
#'
#' # # Estimate distribution parameters from RISK
#' # load('data/clean/phylo_all_genus_with_tree.RData')
#' # phylo_ref <- phylo_all_genus %>%
#' #   subset_samples(collection %in% 'RISK') %>%
#' #   subset_samples(biopsy_location %in% 'stool')
#' # otu_filtered <- phylo_ref %>%
#' #   transform_sample_counts(function(x) x / sum(x)) %>%
#' #   filter_taxa(filterfun(kOverA(5, 1e-04)))
#' # phylo_ref <- prune_taxa(taxa = otu_filtered, x = phylo_ref)
#' # mat_otu_ref <- otu_table(phylo_ref)@.Data
#' # mat_otu_ref_nonzero <- mat_otu_ref
#' # mat_otu_ref_nonzero[mat_otu_ref_nonzero == 0] <- NA
#' # # estimate the feature-mean distribution parameters (modeled as log-normal)
#' # apply(mat_otu_ref_nonzero, 1, mean, na.rm = T) %>% log %>% mean #(~4.5)
#' # apply(mat_otu_ref_nonzero, 1, mean, na.rm = T) %>% log %>% sd #(~2)
#' # # estimate the size parameter for NB distribution
#' # means <- apply(mat_otu_ref_nonzero, 1, mean, na.rm = T)
#' # vars <- apply(mat_otu_ref_nonzero, 1, var, na.rm = T)
#' # fit.lmMeansVars <- lm(log(vars/means) ~ log(means)) #(size ~ 0.3)
#' # # estimate the zero-inflation parameter
#' # perczero <- apply(mat_otu_ref == 0, 1, mean)
#' # fit.lmZI <- lm(log(perczero / (1 - perczero))[perczero != 0] ~
#' #                  log(means)[perczero != 0]) # (intercept~3.3, slop~-0.7)
#' # # estimate distribution of library size (log normal)
#' # librarySizes <- apply(mat_otu_ref, 1, sum)
#' # sd(librarySizes %>% log) # (~2.6)
