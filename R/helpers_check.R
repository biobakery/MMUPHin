#' Check feature abundance table
#'
#' Given a feature abundance table, make sure that a) it has no missing values,
#' b) all values are non-negative, c) it is either proportions (all no greater
#' than 1) or counts (all integers).
#'
#' @param feature_abd feature-by-sample matrix of abundances (proportions or
#' counts).
#'
#' @return returns an error if any of the check fails. Otherwise either "counts"
#' or "proportions"
check_feature_abd <- function(feature_abd) {
  # Errors if has missing values
  if(any(is.na(feature_abd)))
    stop("Found missing values in the feature table!")
  # Errors if has negative values
  # if statement condition is okay because NA's would've been checked
  if(any(feature_abd < 0))
    stop("Found negative values in the feature table!")
  # Returns "proportions" if all values are less than or equal to one
  # Non-negativeness has been checked by previous if statement
  if(all(feature_abd <= 1)) {
    return("proportions")
  }
  # Returns "counts" if all values are integers
  else if(all(feature_abd == floor(feature_abd))) {
    return("counts")
  }
  # Errors if is neither proportions nor counts (i.e., has values greater than
  # one that are not integers)
  else
    stop("Feature table does not appear to be either proportions or counts!")
}

#' Check that sample numbers and names match between a feature table and a
#' metadata data frame
#'
#' Sample names (column names of the feature table, row names of the metadata
#' data frame) must be matching exactly. Note that this dictates that they
#' cannot be NULL because by design data (a data frame) should have non-empty
#' row names.
#'
#' @param feature_abd feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param data data frame of metadata.
#'
#' @return matched sample names
check_samples <- function(feature_abd, data) {
  # Sample numbers need to agree with each other
  if(ncol(feature_abd) != nrow(data))
    stop("Dimensions of feature table and metadata table do not agree!")

  # By designed use case data should always have row names
  if(is.null(rownames(data)))
    stop("data should not have empty row names!")

  # Sample names need to agree with each other
  if(!identical(colnames(feature_abd), rownames(data)))
    stop("Sample names in feature_abd and data don't agree!")

  return(rownames(data))
}

#' Check that metadata data frame has all the variables and not missing
#'
#' @param data data frame of metadata.
#' @param variables name of variables (batch, covariates, etc.) to check
#'
#' @return data reduced to include only those specified in variables
check_metadata <- function(data, variables) {
  # If variables are NULL return NULL (as in the case of no provided covariates)
  if(is.null(variables)) return(NULL)

  # Variables all be present in data (i.e. in the column names)
  variables_absent <- setdiff(variables, colnames(data))
  if(length(variables_absent) > 0) {
    stop("Following variable(s) not present in data:\n",
         paste(variables_absent, collapse = ","))
  }

  # Variables should not have missing values
  variables_missing <- vapply(variables,
                              function(variable) {
                                any(is.na(data[[variable]]))
                              },
                              TRUE)
  if(any(variables_missing)) {
    stop("Following variable(s) in data have missing values:\n",
         paste(variables[variables_missing], collapse = ","))
  }

  return(data[, variables, drop = FALSE])
}

#' Check batch variable
#'
#' @param x batch variable
#' @param min_n_batch min. number of batches (for MMUPHin functions to run)
#'
#' @return if no errors then the batch variables (factorized if not already)
check_batch <- function(x, min_n_batch = 2) {
  # First ensure batch variable is a factor
  if(!is.factor(x)) {
    warning("Batch variable is not a factor and will be converted to one.")
    x <- as.factor(x)
  }

  # Check min number of batches is satisfied
  if(nlevels(x) < min_n_batch)
    stop("Must have at least ", min_n_batch, " batches!")

  return(x)
}
