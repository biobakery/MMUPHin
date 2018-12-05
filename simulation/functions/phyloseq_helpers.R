# These functions are to replace the otu_table, sample_data, and tax_table
# functions in phyloseq, so that the returned objects are matrix, data.frame, and
# matrix, respectively
otu_table2 <- function(phylo) {
  if(!(class(phylo) == "phyloseq"))
    stop("Must have phyloseq object as input!")
  phyloseq::otu_table(phylo)@.Data
}
sample_data2 <- function(phylo) {
  if(!(class(phylo) == "phyloseq"))
    stop("Must have phyloseq object as input!")
  data.frame(phyloseq::sample_data(phylo),
             check.names = FALSE,
             stringsAsFactors = FALSE)
}
tax_table2 <- function(phylo) {
  phyloseq::tax_table(phylo)@.Data
}

# create long format tibble object given phyloseq, aggregating otu table,
# metadata, and tax table together
phyloseq_to_tb <- function(phylo) {
  mat_otu <- otu_table2(phylo)
  df_metadata <- sample_data2(phylo)
  mat_tax <- tax_table2(phylo)
  # No conflicting column names between metadata, tax table, 
  # and "featuer", "rownames", and "abundance"
  check_commonNames <- c(
    intersect(c("feature", "rownames", "abundance"), colnames(df_metadata)) %>% 
      length() %>% is_greater_than(1),
    intersect(c("feature", "rownames", "abundance"), colnames(mat_tax)) %>% 
      length() %>% is_greater_than(1),
    intersect(colnames(df_metadata), colnames(mat_tax)) %>% 
      length() %>% is_greater_than(1)
  )
  if(any(check_commonNames))
    stop("There are overlapping column names, cannot proceed!")
  
  t(mat_otu) %>% 
    tibble::as_tibble(rownames = "rownames") %>% 
    tidyr::gather(key = feature,
                  value = abundance,
                  - rownames) %>% 
    dplyr::left_join(
      mat_tax %>% 
        tibble::as_tibble(rownames = "feature"),
      by = "feature"
    ) %>% 
    dplyr::left_join(
      df_metadata %>% 
        tibble::as_tibble(rownames = "rownames"),
      by = "rownames"
    )
}

# This function is to trim the phyloseq object so that won't run into empty taxa/samples
# which is prone to happenning whenever a phyloseq object is subsetted in any way!
prune_taxaSamples <- function(phylo, 
                              flist_taxa = kOverA2(k = 1, A = 0), # default non-empty pruning
                              flist_samples = function(x) sum(x > 0) > 0, # default non-empty pruning
                              max.iter = 3
                              ) {
  i.iter <- 1
  taxa.ind <- apply(otu_table2(phylo), 1, flist_taxa)
  samples.ind <- apply(otu_table2(phylo), 2, flist_samples)
  if (phyloseq::ntaxa(phylo) != length(taxa.ind) | 
      phyloseq::nsamples(phylo) != length(samples.ind)) {
    stop("Logic error in applying function(s). Logical result not same length as ntaxa(physeq)")
  }
  while(TRUE) {
    if(i.iter > max.iter) 
      stop("Something went wrong! Max iteration exceeded!")
    phylo <- phyloseq::prune_taxa(taxa.ind, phylo)
    samples.ind <- apply(otu_table2(phylo), 2, flist_samples)
    phylo <- phyloseq::prune_samples(samples.ind, phylo)
    taxa.ind <- apply(otu_table2(phylo), 1, flist_taxa)
    if(all(taxa.ind)) return(phylo)
    i.iter <- i.iter + 1
  }
}

# Helper function for filtering taxa
# k = 0 = no filtering
# A = 0, k = 1 = non-empty filtering
# A = 0, k = x>1 = presence filtering
# A = y>0, k = x>1 = abundance filtering
kOverA2 <- function(k = 1, A = 0) {
  if(A >= 1) stop("Filtering is for relative abundance!")
  function(x) {
    if(any(is.na(x))) stop("Missing values in the data!")
    sum(tss(x) > A) >= k
  } 
}

# Helper function for getting relative abundance
# Useful when possible samples are all zero
tss <- function(x) {
  if(all(x == 0)) return(x)
  return(x / sum(x))
}

makeBetterTaxaNames <- function(family, genus, species) {
  # family <- !!family
  # genus <- !!genus
  # species <- !!species
  dplyr::case_when(
    species != "s__" ~ paste(genus %>% 
                               stringr::str_replace_all(stringr::fixed("g__"), ""),
                             species %>% 
                               stringr::str_replace_all(stringr::fixed("s__"), "")),
    genus != "g__" ~ paste(genus %>% 
                             stringr::str_replace_all(stringr::fixed("g__"), ""),
                           "unclassified"),
    family != "f__" ~ paste0(family %>% 
                              stringr::str_replace_all(stringr::fixed("f__"), ""),
                            "(f) unclassified"),
    TRUE ~ "unclassified at family"
  )
}

# This funciton is to aggregate a list of (uniformly prepared) phyloseq objects into
# one single object
combine_phyloseq <- function(l_phylo) {
  studies <- names(l_phylo)
  if(is.null(studies)) stop("Names of l_phylo needed to provide study names!")
  
  # Check taxonamy tables are consistent; create overall taxonamy table
  l_tax <- l_phylo %>% 
    purrr::map(tax_table2)
  ranks1 <- l_tax %>% 
    purrr::map_dbl(~ is.na(.x) %>% 
                     not %>% 
                     apply(2, all) %>% 
                     sum)
  ranks2 <- l_tax %>% 
    purrr::map_dbl(~ is.na(.x) %>% 
                     not %>% 
                     apply(2, any) %>% 
                     sum)
  if(!all(ranks1 == ranks2)) 
    stop("Check the tax tables - some columns have sporadic missing values!")
  if(dplyr::n_distinct(ranks1) > 1)
    stop("Not all phyloseq objects have the same taxonomy level!")

  l_tax_all <- l_tax %>% 
    purrr::map(apply, MARGIN=1, FUN=paste, collapse = "|") 
  if(l_tax_all %>% 
     purrr::map_lgl(~anyDuplicated(.x) > 0) %>% 
     any) stop("Taxonamy needs to be uniquely mappable to features to aggregate!")
  tax_all <- l_tax_all %>% 
    unlist() %>% 
    unique()
  mat_tax_all <- tax_all %>% 
    sapply(strsplit, split = "|", fixed = TRUE) %>% 
    data.frame(check.names = FALSE) %>% 
    as.matrix() %>% 
    t()
  colnames(mat_tax_all) <- colnames(l_tax[[1]])
  
  # Create overall OTU table
  mat_otu_all <- studies %>% 
    purrr::map(function(study) {
      otu_tmp <- otu_table2(l_phylo[[study]])
      rownames(otu_tmp) <- l_tax[[study]] %>% apply(1, paste, collapse = "|")
      otu_new <- matrix(0, nrow = nrow(mat_tax_all), ncol = ncol(otu_tmp))
      dimnames(otu_new) <- list(rownames(mat_tax_all), colnames(otu_tmp))
      otu_new[rownames(otu_tmp), colnames(otu_tmp)] <- otu_tmp
      return(otu_new)
    }) %>% 
    purrr::reduce(cbind)
  
  # Check metadata tables are consistent; create overall metadata table
  l_metadata <- l_phylo %>% 
    purrr::map(sample_data2)
  for(i in 1:(length(l_tax) - 1)) {
    if(!all(colnames(l_metadata[[1]]) == colnames(l_metadata[[i + 1]])))
      stop("Not all phyloseq objects have the same metadata names!")
    if(!all(sapply(l_metadata[[1]], class) == sapply(l_metadata[[i + 1]], class)))
      stop("Not all phyloseq objects have the same metadata types!")
  }
  df_metadata_all <- l_metadata %>% 
    purrr::reduce(rbind)
  
  # Dimension names of the aggregated tables need to agree
  if(!all(colnames(mat_otu_all) == rownames(df_metadata_all))) {
    stop("Sample names in the aggregated OTU table and metadata table do not agree!")
  }
  if(!all(rownames(mat_otu_all) == rownames(mat_tax_all))) {
    stop("Feature names in the aggregated OTU table and metadata table do not agree!")
  }
  
  phyloseq(
    phyloseq::otu_table(mat_otu_all, taxa_are_rows = TRUE),
    phyloseq::sample_data(df_metadata_all),
    phyloseq::tax_table(mat_tax_all)
  )
}
