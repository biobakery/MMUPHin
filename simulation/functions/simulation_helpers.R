create_metadataMatrix <- function(df_metadata,
                                  metadata_type,
                                  fDummyData = TRUE,
                                  fDummyDirection = TRUE,
                                  scale = FALSE) {
  if(any(colnames(df_metadata) != names(metadata_type)))
    stop("Variable names in df_metadata and metadata_type do not agree!")
  if(any(sapply(df_metadata, class) != metadata_type))
    stop("Variable classes in df_metadata and metadata_type do not agree!")
  l_mat_metadata <- lapply(names(metadata_type), function(variable) {
    if(!fDummyData | metadata_type[variable] != "factor") {
      mat_tmp <- rbind(df_metadata[, variable])
      rownames(mat_tmp) <- variable
      return(mat_tmp)
    } else {
      lvls <- levels(df_metadata[, variable])
      nlvls <- nlevels(df_metadata[, variable])
      mat_tmp <- t(sapply(1:nlvls, function(ilvl) {
        (df_metadata[, variable] == lvls[ilvl]) * 1
      }))
      rownames(mat_tmp) <- paste0(variable, "_", 1:nlvls)
      if(fDummyDirection) {
        mat_tmp_opposite <- t(sapply(1:nlvls, function(ilvl) {
          (df_metadata[, variable] == lvls[ilvl]) * (-1)
        }))
        rownames(mat_tmp_opposite) <- paste0(variable, "_", 1:nlvls, "_opposite")
        mat_tmp <- rbind(mat_tmp, mat_tmp_opposite)
      }
      return(mat_tmp)
    }
  })
  mat_metadata <- Reduce("rbind", l_mat_metadata)
  if(scale) {
    mean_metadata <- apply(mat_metadata, 1, mean)
    sd_metadata <- apply(mat_metadata, 1, sd)
    sd_metadata[sd_metadata == 0] <- 1
    mat_metadata <- (mat_metadata - mean_metadata) / sd_metadata
  }
  return(mat_metadata)
}

# this is to help create metadata matrix that simulates "random effect"
# of exposure variables across batches, by differentially scaling the
# exposure within batches
create_metadataMatrix_RE <- function(df_metadata,
                                     metadata_type,
                                     batch = "batch", #name of the batch variable
                                     exposure = "exposure", #name of the exposure variable
                                     rel.mean,
                                     fDummyData = TRUE,
                                     fDummyDirection = TRUE) {
  if(any(colnames(df_metadata) != names(metadata_type)))
    stop("Variable names in df_metadata and metadata_type do not agree!")
  if(any(sapply(df_metadata, class) != metadata_type))
    stop("Variable classes in df_metadata and metadata_type do not agree!")
  if(!all(c(batch, exposure) %in% colnames(df_metadata)))
    stop("Batch and exposure variable names not in data!")
  if(metadata_type[batch] != "factor")
    stop("Batch variable must be factor!")

  batch <- df_metadata[, batch, drop = TRUE]
  if(length(rel.mean) != nlevels(batch))
    stop("Length of relative effects is different from number of batches!")
  l_mat_metadata <- lapply(names(metadata_type), function(variable) {
    if(!fDummyData | metadata_type[variable] != "factor") {
      mat_tmp <- rbind(df_metadata[, variable, drop = TRUE])
      rownames(mat_tmp) <- variable
      if(variable == exposure) mat_tmp <- mat_tmp %>%
        apply(1, scale_RE, batch = batch, rel.mean = rel.mean) %>%
        t()
      return(mat_tmp)
    } else {
      lvls <- levels(df_metadata[, variable])
      nlvls <- nlevels(df_metadata[, variable])
      mat_tmp <- t(sapply(1:nlvls, function(ilvl) {
        (df_metadata[, variable] == lvls[ilvl]) * 1
      }))
      rownames(mat_tmp) <- paste0(variable, "_", 1:nlvls)
      if(fDummyDirection) {
        mat_tmp_opposite <- t(sapply(1:nlvls, function(ilvl) {
          (df_metadata[, variable] == lvls[ilvl]) * (-1)
        }))
        rownames(mat_tmp_opposite) <- paste0(variable, "_", 1:nlvls, "_opposite")
        mat_tmp <- rbind(mat_tmp, mat_tmp_opposite)
      }
      if(variable == exposure) mat_tmp <- mat_tmp %>%
        apply(1, scale_RE, batch = batch, rel.mean = rel.mean) %>%
        t()
      return(mat_tmp)
    }
  })
  mat_metadata <- Reduce("rbind", l_mat_metadata)

  return(mat_metadata)
}

scale_RE <- function(x, batch, rel.mean) {
  x.scaled <- x - mean(x)
  lvl.batch <- levels(batch)
  for(i in 1:nlevels(batch)) {
    x.scaled[batch == lvl.batch[i]] <- x.scaled[batch == lvl.batch[i]] + rel.mean[i]
  }
  return(x.scaled + mean(x))
}

create_effectSize <- function(effectSize,
                              df_metadata,
                              metadata_type,
                              fDummyData = TRUE,
                              fDummyDirection = TRUE) {
  if(any(names(effectSize) != names(metadata_type)))
    stop("Variable names in effectSize and metadata_type do not agree!")
  if(any(colnames(df_metadata) != names(metadata_type)))
    stop("Variable names in df_metadata and metadata_type do not agree!")
  if(any(sapply(df_metadata, class) != metadata_type))
    stop("Variable classes in df_metadata and metadata_type do not agree!")
  l_effectSize <- lapply(names(metadata_type), function(variable) {
    if(!fDummyData | metadata_type[variable] != "factor") {
      effectSize_tmp <- effectSize[variable]
      names(effectSize_tmp) <- variable
      return(effectSize_tmp)
    } else {
      nlvls <- nlevels(df_metadata[, variable])
      effectSize_tmp <- rep(effectSize[variable], nlvls)
      names(effectSize_tmp) <- paste0(variable, "_", 1:nlvls)
      if(fDummyDirection) {
        effectSize_tmp_opposite <- rep(effectSize[variable], nlvls)
        names(effectSize_tmp) <- paste0(variable, "_", 1:nlvls, "_opposite")
        effectSize_tmp <- c(effectSize_tmp, effectSize_tmp_opposite)
      }
      return(effectSize_tmp)
    }
  })
  return(unlist(l_effectSize))
}
