create_metadataMatrix <- function(df_metadata,
                                  metadata_type,
                                  scale = FALSE,
                                  fDummyData = FALSE) {
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
      mat_tmp <- t(sapply(2:nlvls, function(ilvl) {
        (df_metadata[, variable] == lvls[ilvl]) * 1
      }))
      rownames(mat_tmp) <- paste0(variable, "_", 2:nlvls)
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

create_effectSize <- function(effectSize, df_metadata, metadata_type, fDummyData = FALSE) {
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
      effectSize_tmp <- rep(effectSize[variable], nlvls - 1)
      names(effectSize_tmp) <- paste0(variable, "_", 2:nlvls)
      return(effectSize_tmp)
    }
  })
  return(unlist(l_effectSize))
}

create_fronzenIndex <- function(df_metadata, metadata_type, fDummyData = FALSE) {
  if(any(colnames(df_metadata) != names(metadata_type)))
    stop("Variable names in df_metadata and metadata_type do not agree!")
  if(any(sapply(df_metadata, class) != metadata_type))
    stop("Variable classes in df_metadata and metadata_type do not agree!")
  ns <- sapply(names(metadata_type), function(variable) {
    if(!fDummyData | metadata_type[variable] != "factor") {
      return(1)
    } else {
      nlvls <- nlevels(df_metadata[, variable])
      return(nlvls - 1)
    }
  })
  return(1:sum(ns))
}
