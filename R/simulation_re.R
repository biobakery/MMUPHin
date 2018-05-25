library(sparseDOSSA)
library(doParallel)
library(parallel)

##############################
## Synthetic Data Generation #
##############################

# Generate Replicated Simulated Datasets For A Combination of Parameters
trigger_sparseDOSSA_Simulator<-function(noZeroInflate=FALSE,
                                        RandomEffect=FALSE,
                                        metadataType,
                                        nSubjects,
                                        nPerSubject,
                                        nMicrobes,
                                        spikeMicrobes,
                                        nMetadata,
                                        spikeMetadata,
                                        effectSize,
                                        nIterations=100,
                                        rSeed=1234,
                                        nCores=4){

  # Create Replicates
  reps = 1:nIterations

  ########################
  # Catch Obvious Errors #
  ########################

  # Character Values
  if (!metadataType %in% c('UVA', 'UVB', 'MVA', 'MVB', 'MVC', 'MVD', 'MVE', 'MVF'))
    stop('Must be one of the following: UVA, UVB, MVA, MVB, MVC, MVD, MVE, MVF.')

  # Positive Integer Values
  if (round(nSubjects) != nSubjects ||
      nSubjects<0 ||
      round(nPerSubject) != nPerSubject ||
      nPerSubject<0 ||
      round(nMicrobes) != nMicrobes ||
      nMicrobes<0 ||
      round(nMetadata) != nMetadata ||
      nMetadata<0)
    stop('nSubjects/nPerSubject/nMicrobes/nMetadata must be positive integers.')

  # Proportion Values
  if (spikeMicrobes>1 || spikeMicrobes<=0 || spikeMetadata<=0 || spikeMetadata>1)
    stop('spikeMicrobes/spikeMetadata must be in (0, 1].')

  # Illegal Combinations
  if(RandomEffect==TRUE && nPerSubject==1)
    stop('nPerSubject must be greater 1 when RandomEffect is TRUE.')

  if(metadataType %in% c('UVA', 'UVB') && (nMetadata!=1 || spikeMetadata!=1))
    stop('Both nMetadata and spikeMetadata must be equal to 1 when metadataType is UVA, UVB.')

  if(!metadataType %in% c('UVA', 'UVB') && nMetadata==1)
    stop('nMetadata must be greater than 1 when metadataType is MVA, MVB, MVC, MVD, MVE, MVF.')

  # Define the simulation parameters combinations
  simparams = apply(expand.grid(metadataType,
                                nSubjects,
                                nPerSubject,
                                nMicrobes,
                                spikeMicrobes,
                                nMetadata,
                                spikeMetadata,
                                effectSize,
                                reps), 1, paste, collapse = '_')

  # Define the labels to go with each element of the simulation parameter
  # after splitting on the delimiter
  simparamslabels = c("metadataType","nSubjects", "nPerSubject", "nMicrobes", "spikeMicrobes", "nMetadata", "spikeMetadata", "effectSize", "rep")

  # Set Up Clustering Environment
  no_cores <- nCores
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  # Track Start Time
  cat(c("Job started at:",date()), "\n")
  start.time <- Sys.time()

  ####################
  # Data Generation #
  ###################

  set.seed(rSeed) # For reproducibility

  # Call SparseDOSSA Wrapper
  simlist <- sparseDOSSA_Wrapper(simparams, simparamslabels, noZeroInflate=noZeroInflate)

  # Set Names
  if (noZeroInflate==TRUE && RandomEffect==TRUE) {
    simnames<- paste('noZeroInflate_RandomEffect', simparams, sep='_')}
  if (noZeroInflate==TRUE && RandomEffect==FALSE) {
    simnames<- paste('noZeroInflate_noRandomEffect', simparams, sep='_')}
  if (noZeroInflate==FALSE && RandomEffect==TRUE) {
    simnames<- paste('ZeroInflate_RandomEffect', simparams, sep='_')}
  if (noZeroInflate==FALSE && RandomEffect==FALSE) {
    simnames<- paste('ZeroInflate_noRandomEffect', simparams, sep='_')}
  names(simlist) <- simnames

  # Stop the Cluster
  stopCluster(cl)

  # Track End Time
  cat(c("Job ended at:",date()), "\n")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"),3)
  cat("Computational time:", minutes, "minutes \n")

  # Return
  return(simlist)
}

# Trigger sparseDOSSA
sparseDOSSA_Wrapper<-function(simparams, simparamslabels, noZeroInflate){
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



#####################
# Generate Metadata #
#####################

generateMetadata<-function(metadataType,
                           nSubjects,
                           nPerSubject,
                           nMetadata,
                           spikeMetadata){

  # Calculate Number of Samples
  nSamples = round(nSubjects*nPerSubject)

  # Create Blocking Variable
  if (nPerSubject==1){  # NO RANDOM EFFECTS
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=0))
  }
  if (nPerSubject>1){  # SUBJECT-SPECIFIC RANDOM EFFECTS
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=1))
  }
  BLOCK <- as.vector(matrix(subjectRandomEffects,nrow=nPerSubject,ncol=length(subjectRandomEffects),byrow=TRUE))

  # Specify Mean and Covariance Structure
  mu<-rep(0,nMetadata)
  cov<-diag(1,nMetadata, nMetadata)

  if (metadataType %in% c('MVD', 'MVE', 'MVF')){
    for (i in 1:nMetadata){
      for (j in 1:nMetadata){
        if(i!=j) cov[i,j]=0.5
      }
    }
  }

  # Generate from MVN
  fakeMetadata<-as.matrix(mvrnorm(n=nSamples, mu,cov))

  # Transpose and Add Blocking Structure
  finalMetadata<-apply(fakeMetadata, 2, function(x) x+BLOCK)

  # Dichotomization for Binary Cases
  if (metadataType %in% c('UVB', 'MVB', 'MVE')){
    UserMetadata<-t(apply(finalMetadata, 2, function(x) ifelse(x>median(x), 1, 0)))
  }

  if (metadataType %in% c('MVC', 'MVF')){ # Mixed Scenario - Dichotomize Half of the Features
    t_UserMetadata<-apply(finalMetadata, 2, function(x) ifelse(x>median(x), 1, 0))
    columns_not_to_binarize<-sample(1:nMetadata, nMetadata/2)
    t_UserMetadata[,columns_not_to_binarize]<-finalMetadata[, columns_not_to_binarize]
    UserMetadata<-t(t_UserMetadata)
  }

  if (metadataType %in% c('UVA', 'MVA', 'MVD')){
    UserMetadata<-t(finalMetadata)
  }

  # Collect Relevant Spike-in Information
  spikeCount<- round(nMetadata*spikeMetadata)
  Metadatafrozenidx<-sample(1:nMetadata, spikeCount, replace=FALSE)
  significant_metadata<-paste('Metadata', Metadatafrozenidx, sep='')
  spikeCount<-as.character(spikeCount)

  # Return
  return(list(UserMetadata=UserMetadata,
              Metadatafrozenidx=Metadatafrozenidx,
              significant_metadata=significant_metadata,
              spikeCount=spikeCount))
}
