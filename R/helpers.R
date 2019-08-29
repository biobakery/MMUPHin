## Utilities for normalization/transformation

#' Set pseudo count for an abundance matrix. Pseudo count is currently set to
#' half of minimum non-zero values
#'
#' @param features feature-by-sample matrix of abundances (proportions or
#' counts).
#'
#' @return the pseudo count
set_pseudo <- function(features) {
  type_features <- check_feature_abd(features)
  if(all(features == 0))
    stop("All feature abundances are zero!")

  min(setdiff(features, 0)) / 2
}

#' Match user-specified control parameters with default, and modify if needed
#'
#' @param default list of default control parameters
#' @param control list of user-provided control parameters
#'
#' @return list of control parameters, set to user provided values if specified
#' and default other wise
match_control <- function(default, control) {
  if (missing(control))
    control <- list()

  control_pos <- pmatch(names(control), names(default))
  default[c(na.omit(control_pos))] <- control[!is.na(control_pos)]

  return(default)
}

## These are modified from Maaslin2

#' Normalize feature abundance table
#'
#' @param features feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param normalization normalization method.
#' @param pseudo_count pseudo count to be added to feature_abd.
#'
#' @return normalized abundance table.
normalize_features <- function(features,
                               normalization = "NONE",
                               pseudo_count = 0) {
  type_features <- check_feature_abd(features)
  features <- features + pseudo_count
  if (normalization=='TSS')
    features <- apply(features, 2, TSS)
  if (normalization=='NONE')
    features <- features
  return(features)
}

#' Transform feature abunadnce table
#'
#' @param features feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param transform transformation method.
#' @param pseudo_count pseudo count to be added to feature_abd..
#'
#' @return transformed abundance table.
transform_features <- function(features,
                               transform = "NONE",
                               pseudo_count = 0) {
  type_features <- check_feature_abd(features)
  if(type_features != "proportions")
    stop("Transformation should only be applied to normalized features",
         " (proportions)" )
  features <- features + pseudo_count
  if (transform =='LOG') {
    if(any(features <= 0))
      stop("LOG transformation not applicable to values smaller than or equal",
           " to zero!")
    features <- apply(features, 2, LOG)
  }
  if (transform =='AST')
    features <- apply(features, 2, AST)
  if (transform =='NONE')
    features <- features
  return(features)
}

#' TSS normalization
#'
#' @param x vector of abundance to be normalized.
#'
#' @return normalized vector of abundance.
TSS <- function(x) {
  if(all(x == 0)) return(x)
  return(x / sum(x))
}

#' AST transformation
#'
#' @param x vector of abundance to be transformed.
#'
#' @return transformed vector of abundance.
AST<-function(x){
  return(asin(sqrt(x)))
}

#' LOG transformation. This differs from Maaslin2
#'
#' @param x vector of abundance to be transformed.
#'
#' @return transformed vector of abundance.
LOG<-function(x){
  return(log2(x))
}

