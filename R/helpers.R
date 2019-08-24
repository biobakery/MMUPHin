## Utilities for normalization/transformation
## These are modified from Maaslin2

#' Normalize feature abundance table
#'
#' @param features feature*sample abundance table.
#' @param normalization normalization method.
#' @param pseudo.count pseudo count to be added to abundances before normalization.
#'
#' @return normalized feature*sample abundance table.
normalizeFeatures <- function(features,
                              normalization = "NONE",
                              pseudo.count = 0) {
  features <- features + pseudo.count
  if (normalization=='TSS')
  {
    features<-apply(features, 2, tss)
  }
  if (normalization=='NONE')
  {
    features<-features
  }
  return(features)
}

#' Transform feature abunadnce table
#'
#' @param features feature*sample abundance table.
#' @param transform transformation method.
#' @param pseudo.count pseudo count to be added to abundances before normalization.
#' @return transformed feature*sample abundance table.
transformFeatures <- function(features,
                              transform = "NONE",
                              pseudo.count = 0) {
  features <- features + pseudo.count
  if (transform =='LOG')   {
    features <- apply(features, 2, LOG)
  }
  if (transform =='AST')   {
    features <- apply(features, 2, AST)
  }
  if (transform =='NONE')   {
    features <- features
  }
  return(features)
}

#' TSS normalization
#'
#' @param x vector of abundance to be normalized.
#'
#' @return normalized vector of abundance.
tss <- function(x) {
  if(all(x == 0)) return(x)
  return(x / sum(x))
}

#' AST transformation
#'
#' @param x vector of abundance to be transformed.
#'
#' @return transformed vector of abundance.
AST<-function(x){
  return(sign(x)*asin(sqrt(abs(x))))
}

#' LOG transformation. This differs from Maaslin2
#'
#' @param x vector of abundance to be transformed.
#'
#' @return transformed vector of abundance.
LOG<-function(x){
  return(log2(x))
}

## Utlities for empirical Bayes estimation for batch.adj
## These are modified from sva


