load(
  "results/quick_and_easy/qe_biom.RData"
)
qe_biom_genus_tmp <- qe_biom_genus %>% subset_samples(Project != "BIDMC-FMT")
feature.count <- qe_biom_genus_tmp %>% otu_table() %>% `@`(`.Data`)
batch <- sample_data(qe_biom_genus_tmp)$Project

ordination <- ordinate(qe_biom_genus_tmp %>%
                         transform_sample_counts(function(x)x/sum(x)),
                       distance = "bray",
                       method = "MDS")
plot_ordination(qe_biom_genus_tmp, ordination, color = "Project") +
  theme_bw()

studies.all.stool <- table(sample_data(all_stool.phyloseq)$studyID)
studies.all.stool <- names(studies.all.stool)[studies.all.stool > 50]
all_stool.phyloseq <- all_stool.phyloseq %>%
  subset_samples(studyID %in% studies.all.stool)
all_stool.phyloseq_ra <- all_stool.phyloseq %>%
  transform_sample_counts(function(x)x/sum(x))
bugs_abd <- all_stool.phyloseq_ra %>%
  filter_taxa(kOverA(5, 2e-05))
all_stool.phyloseq <- all_stool.phyloseq %>% prune_taxa(bugs_abd, .)
feature.count <- all_stool.phyloseq %>% otu_table() %>% `@`(`.Data`)
batch <- sample_data(all_stool.phyloseq)$studyID

feature.count.mod <- MMUPHin::adjust.batch(
  feature.count,
  batch
)
all_stool.phyloseq.batch <- all_stool.phyloseq
otu_table(all_stool.phyloseq.batch) <- otu_table(feature.count.mod, taxa_are_rows = TRUE)
ordination <- ordinate(all_stool.phyloseq.batch %>%
                         transform_sample_counts(function(x)x/sum(x)),
                       distance = "bray",
                       method = "MDS")
plot_ordination(all_stool.phyloseq.batch, ordination, color = "studyID") +
  theme_bw()


feature.count <- feature.count.mod
pseudo.count = 0.5
distance = 2
k.max = 7
cluster.method = fpc::claraCBI
classify.method = "averagedist"
M = 10
nnk = 1
diagnostics = TRUE
verbose = TRUE

  ## Ensure data formatts are as expected
feature.count <- as.matrix(feature.count)
if(is.null(rownames(feature.count)))
  rownames(feauture.count) <- paste0("Feature", 1:nrow(feature.count))
if(any(feature.count < 0, na.rm = TRUE))
  stop("Found negative values in the feature table!")
if(any(is.na(feature.count)))
  stop("Found missing values in the feature table!")
batch <- as.factor(batch)
if(any(is.na(batch))) stop("Found missing values in the batch variable!")

## Data dimensions need to agree with each other
if(ncol(feature.count) != length(batch)) {
  stop("Dimensions of feature table and batch variable do not agree!")
}

lvl.batch <- levels(batch)
n.batch <- length(lvl.batch)
if(n.batch < 3) stop("Only discovery based on three or more studies is supported!")
if(verbose) message("Found ", n.batch, " batches/studies.")

# Number of max clusters to evaluate can't exceed half of smallest sample size
if(k.max > floor(min(table(batch)) / 2) - 1) {
  message("k.max be less than half smallest sample size!")
  k.max <- floor(min(table(batch)) / 2) - 1
  message("Set k.max to ", k.max)
}



if(verbose) message("Calculating all vs. all dissimilarity matrix...")
dist.all <- as.matrix(phyloseq::distance(
  phyloseq::otu_table(apply(feature.count, 2, function(x) x / sum(x)),
                      taxa_are_rows = TRUE),
  method = distance))

# Clustering validation
stats.internal <- replicate(k.max - 1, list(list()))
stats.external <- replicate(k.max - 1, list(list()))
for(k in 2:k.max) {
  if(verbose) message("Now evaluating clustering for k = ", k, "...")

  # Cluster each study first individually
  clusterings <- list()
  for(i in 1:n.batch) {
    i.batch <- lvl.batch[i]
    clusterings[[i]] <- cluster.method(
      data = as.dist(dist.all[batch == i.batch, batch == i.batch]),
      k = k
    )
  }

  # Validation
  classifications <- list()
  for(i in 1:n.batch) {
    i.batch <- lvl.batch[i]
    if(verbose) message("Study is ", i.batch)
    if(verbose) message("Performing internal validation...")
    i.result.internal <- prediction.strength2(
      xdata = as.dist(dist.all[batch == i.batch, batch == i.batch]),
      Gmin = k,
      Gmax = k,
      M = M,
      clustermethod = cluster.method,
      classification = classify.method,
      distances = TRUE,
      count = FALSE
    )
    stats.internal[[k - 1]][[i]] <- i.result.internal$predcorr[[k]]

    if(verbose) message("Performing external validation...")
    i.clustering <- rep(-1, ncol(feature.count))
    i.clustering[batch == i.batch] <- clusterings[[i]]$partition
    classifications[[i]] <- fpc::classifdist(as.dist(dist.all),
                                             clustering = i.clustering,
                                             method = classify.method,
                                             centroids = clusterings[[i]]$result$medoids,
                                             nnk = nnk)

    i.result.external <- matrix(NA, n.batch, k)
    for (j in setdiff(1:n.batch, i)) {
      j.batch <- lvl.batch[[j]]
      ctable <- fpc::xtable(clusterings[[j]]$partition,
                            classifications[[i]][batch == j.batch], k)
      for (kk in 1:k) {
        i.result.external[j, kk] <- sum(ctable[kk, ]^2 - ctable[kk,])
        cpjk <- clusterings[[j]]$partition == kk
        njk <- sum(cpjk)
        if (njk > 1)
          i.result.external[j, kk] <- i.result.external[j, kk]/(njk * (njk - 1))
        else i.result.external[j, kk] <- NA
      }
    }
    i.result.external <- apply(i.result.external[setdiff(1:n.batch, i), ],
                               1, min, na.rm = TRUE)
    stats.external[[k - 1]][[i]] <- i.result.external
  }
}

if(diagnostics) {
  data.internal <- Reduce("rbind",
                          lapply(1:(k.max - 1),
                                 function(i.k) {
                                   data.frame(k = i.k + 1,
                                              Reduce("rbind",
                                                     lapply(1:n.batch,
                                                            function(i)
                                                              data.frame(pred.strength = stats.internal[[i.k]][[i]],
                                                                         batch = lvl.batch[i]))))
                                 }))
  data.external <- Reduce("rbind",
                          lapply(1:(k.max - 1),
                                 function(i.k) {
                                   data.frame(k = i.k + 1,
                                              Reduce("rbind",
                                                     lapply(1:n.batch,
                                                            function(i)
                                                              data.frame(pred.strength = stats.external[[i.k]][[i]],
                                                                         batch = lvl.batch[i]))))
                                 }))
  data.plot <- rbind(data.frame(data.internal, validation = "Internal"),
                     data.frame(data.external, validation = "External"))
  print(ggplot2::ggplot(data.plot, ggplot2::aes(x = as.factor(k), y = pred.strength)) +
          ggplot2::geom_boxplot(ggplot2::aes(color = validation)) +
          ggplot2::facet_grid(.~batch %>% gsub("metaphlan_bugs_list.stool", "", ., fixed = TRUE)) +
          ggplot2::theme_bw())
}
  return(list(internal = stats.internal,
              external = stats.external))
}


prediction.strength2 <- function (xdata, Gmin = 2, Gmax = 10, M = 50, clustermethod = kmeansCBI,
          classification = "centroid", cutoff = 0.8, nnk = 1, distances = inherits(xdata,
                                                                                   "dist"), count = FALSE, ...)
{
  xdata <- as.matrix(xdata)
  n <- nrow(xdata)
  nf <- c(floor(n/2), n - floor(n/2))
  indvec <- clcenters <- clusterings <- jclusterings <- classifications <- list()
  corrpred <- list()
  for (k in Gmin:Gmax) {
    if (count)
      cat(k, " clusters\n")
    corrpred[[k]] <- numeric(0)
    for (l in 1:M) {
      nperm <- sample(n, n)
      if (count)
        cat(" Run ", l, "\n")
      indvec[[l]] <- list()
      indvec[[l]][[1]] <- nperm[1:nf[1]]
      indvec[[l]][[2]] <- nperm[(nf[1] + 1):n]
      for (i in 1:2) {
        if (distances)
          clusterings[[i]] <- clustermethod(as.dist(xdata[indvec[[l]][[i]],
                                                          indvec[[l]][[i]]]), k)
        else clusterings[[i]] <- clustermethod(xdata[indvec[[l]][[i]],
                                                     ], k, ...)
        jclusterings[[i]] <- rep(-1, n)
        jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$partition
        centroids <- NULL
        if (classification == "centroid") {
          if (identical(clustermethod, fpc::kmeansCBI))
            centroids <- clusterings[[i]]$result$centers
          if (identical(clustermethod, fpc::claraCBI))
            centroids <- clusterings[[i]]$result$medoids
        }
        j <- 3 - i
        if (distances)
          classifications[[j]] <- fpc::classifdist(as.dist(xdata),
                                              jclusterings[[i]], method = classification,
                                              centroids = centroids, nnk = nnk)[indvec[[l]][[j]]]
        else classifications[[j]] <- fpc::classifnp(xdata,
                                               jclusterings[[i]], method = classification,
                                               centroids = centroids, nnk = nnk)[indvec[[l]][[j]]]
      }
      ps <- matrix(0, nrow = 2, ncol = k)
      for (i in 1:2) {
        ctable <- fpc::xtable(clusterings[[i]]$partition,
                         classifications[[i]], k)
        for (kk in 1:k) {
          ps[i, kk] <- sum(ctable[kk, ]^2 - ctable[kk,
                                                   ])
          cpik <- clusterings[[i]]$partition == kk
          nik <- sum(cpik)
          if (nik > 1)
            ps[i, kk] <- ps[i, kk]/(nik * (nik - 1))
          else ps[i, kk] <- NA
        }
      }
      corrpred[[k]][l] <- mean(c(min(ps[1, ], na.rm = TRUE), min(ps[2,
                                                      ], na.rm = TRUE)))
    }
  }
  mean.pred <- numeric(0)
  if (Gmin > 1)
    mean.pred <- c(1)
  if (Gmin > 2)
    mean.pred <- c(mean.pred, rep(NA, Gmin - 2))
  for (k in Gmin:Gmax) mean.pred <- c(mean.pred, mean(corrpred[[k]]))
  optimalk <- max(which(mean.pred > cutoff))
  out <- list(predcorr = corrpred, mean.pred = mean.pred,
              optimalk = optimalk, cutoff = cutoff, method = clusterings[[1]]$clustermethod,
              Gmax = Gmax, M = M)
  class(out) <- "predstr"
  out
}

