sparseDOSSA_Wrapper <- function(simparams, simparamslabels,
                              noZeroInflate,
                              ncore ){
  foreach(i = simparams, .packages = c("sparseDOSSA", "MASS", "stringi"),
          .export = c("generateMetadata"), .errorhandling = 'remove') %dopar% {

            params = strsplit(i, '_')[[1]]
            names(params) <- simparamslabels

            # Extract Relevant Parameters
            metadataType = as.character(params["metadataType"]) # Type of Metadata
            nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
            nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
            nSamples<-round(nSubjects*nPerSubject) # Number of Samples
            nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
            spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Percentage of Spiked-in Microbes
            nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
            spikeMetadata<-as.numeric(params["spikeMetadata"])  # Percentage of Spiked-in Metadata
            effectSize<-as.character(params["effectSize"]) # Effect Size


            # Initialize
            DD = NULL

            # Error control using `try`,
            tryAgain = TRUE
            infiniteloopcounter = 1
            while (tryAgain & infiniteloopcounter < 5) {

              # Generate Metadata
              FF<-generateMetadata(metadataType=metadataType,
                                   nSubjects=nSubjects,
                                   nPerSubject=nPerSubject,
                                   nMetadata=nMetadata,
                                   spikeMetadata=spikeMetadata)

              # Extract Relevant Information
              UserMetadata<-FF$UserMetadata;
              Metadatafrozenidx<-FF$Metadatafrozenidx;
              significant_metadata<-FF$significant_metadata;
              spikeCount<-FF$spikeCount

              # Generate sparseDOSSA synthetic abundance
              DD<-sparseDOSSA::sparseDOSSA(number_features = nMicrobes,
                                           number_samples = nSamples,
                                           UserMetadata = UserMetadata,
                                           Metadatafrozenidx = Metadatafrozenidx,
                                           datasetCount = 1,
                                           spikeCount = spikeCount,
                                           spikeStrength = effectSize,
                                           noZeroInflate=noZeroInflate,
                                           percent_spiked=spikeMicrobes)
              if (is.null(DD) | inherits(DD, "try-error")) {
                tryAgain = TRUE
                infiniteloopcounter = infiniteloopcounter + 1
              } else {
                tryAgain = FALSE
              }
            }
            if (infiniteloopcounter >= 5) {
              stop("Consistent error found during simulation. Need to investigate cause.")
            }

            # Gather sparseDOSSA outputs
            sparsedossa_results <- as.data.frame(DD$OTU_count)
            rownames(sparsedossa_results)<-sparsedossa_results$X1
            sparsedossa_results<-sparsedossa_results[-1,-1]
            colnames(sparsedossa_results)<-paste('Sample', 1:ncol(sparsedossa_results), sep='')
            data<-as.matrix(sparsedossa_results[-c((nMetadata+1):(2*nMicrobes+nMetadata)),])
            data<-data.matrix(data)
            class(data) <- "numeric"
            truth<-c(unlist(DD$truth))
            truth<-truth[!stri_detect_fixed(truth,":")]
            truth<-truth[(5+nMetadata):length(truth)]
            truth<-as.data.frame(truth)
            significant_features<-as.vector(truth[seq(1, (as.numeric(spikeCount)+1)*(nMicrobes*spikeMicrobes), (as.numeric(spikeCount)+1)),])

            # Separate Metadata and Taxa - Unified Use of Maaslin and metagenomeSeq

            # Extract Metadata
            if (metadataType %in% c('UVA', 'UVB')){
              metadata<-as.data.frame(data[1,])
              colnames(metadata)<-rownames(data)[1]
            } else{
              metadata<-as.data.frame(t(data[(1:nMetadata),]))
            }

            # Rename the Metadata and True Positive Features - Same Format at Mcmurdie and Holmes (2014)
            which.TP = colnames(metadata) %in% significant_metadata
            meta_newname = paste0(colnames(metadata)[which.TP], "_TP")
            colnames(metadata)[which.TP] <- meta_newname

            # Extract Features
            features<-as.data.frame(t(data[-c(1:nMetadata),]))

            # Rename the Features and True Positive Features - Same Format at Mcmurdie and Holmes (2014)
            wh.TP = colnames(features) %in% significant_features
            colnames(features)<-paste("Feature", 1:nMicrobes, sep = "")
            newname = paste0(colnames(features)[wh.TP], "_TP")
            colnames(features)[wh.TP] <- newname;

            # Add Sample ID
            ID<-rep(paste('Subject', 1:nSubjects, sep=''), each = nPerSubject)

            # Add Normalization or Scaling Factor (Sequencing Depth)
            libSize<-rowSums(features)

            # Return as list
            return(list(metadata=metadata, features=features, ID=ID, libSize=libSize))
          }
}


#' Title
#'
#' @param sparseDOSSA_fit
#' @param add_libsize_var
#'
#' @return
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
extract_sparseDOSSA <- function(sparseDOSSA_fit, add_libsize_var = FALSE) {

  # metadata + feature data
  sparsedossa_results <- as.data.frame(sparseDOSSA_fit$OTU_count)
  rownames(sparsedossa_results) <- sparsedossa_results$X1
  nMetadata <- sum(grepl("Metadata", sparsedossa_results$X1, fixed = TRUE))
  nMicrobes <- sum(grepl("Feature_spike", sparsedossa_results$X1, fixed = TRUE))
  nSamples <- ncol(sparsedossa_results) - 1
  sparsedossa_results <- sparsedossa_results[-1, -1]
  colnames(sparsedossa_results) <- paste('Sample', 1:nSamples, sep='')
  data <- as.matrix(sparsedossa_results[-c((nMetadata+1):(2*nMicrobes+nMetadata)), ])
  data <- data.matrix(data)
  class(data) <- "numeric"

  # Spiked-in features and metadata
  truth <- c(unlist(sparseDOSSA_fit$truth))
  truth <- truth[!stringi::stri_detect_fixed(truth,":")]
  significant_features <- truth[grepl("Feature", truth, fixed = TRUE)]
  significant_metadata <- truth[-(1+nMetadata)] %>%
    stringr::str_subset("Metadata") %>%
    stringr::str_replace("_Level_.+", "") %>%
    unique

  # Extract Metadata
  metadata <- as.data.frame(t(data[(1:nMetadata), ]))

  # Rename True Positive Metadata - Same Format at Mcmurdie and Holmes (2014)
  which.TP <- colnames(metadata) %in% significant_metadata
  meta_newname <- paste0(colnames(metadata)[which.TP], "_TP")
  colnames(metadata)[which.TP] <- meta_newname

  # Extract Features
  features <- as.data.frame(t(data[-c(1:nMetadata),]))

  # Rename Features and True Positive Features - Same Format at Mcmurdie and Holmes (2014)
  wh.TP <- colnames(features) %in% significant_features
  colnames(features) <- paste("Feature", 1:nMicrobes, sep = "")
  newname <- paste0(colnames(features)[wh.TP], "_TP")
  colnames(features)[wh.TP] <- newname

  # If specified, further randomize library size
  if(add_libsize_var) {
    features <- features * exp(rnorm(n = nSamples, mean = 0.5, sd = 0.5))
  }
  # feature table has rows as features
  features <- t(features)
  libSize <- colSums(features)

  # Return as list
  return(list(metadata=metadata, features=features, libSize=libSize))
}
