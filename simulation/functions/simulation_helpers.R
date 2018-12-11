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
