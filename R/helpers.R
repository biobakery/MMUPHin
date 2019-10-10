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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
check_metadata <- function(data, variables, no_missing = TRUE) {
  # If variables are NULL return NULL (as in the case of no provided covariates)
  if(is.null(variables)) return(NULL)

  # Variables must all be present in data (i.e. in the column names)
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
  if(any(variables_missing) & no_missing) {
    stop("Following variable(s) in data have missing values:\n",
         paste(variables[variables_missing], collapse = ","))
  }

  return(data[, variables, drop = FALSE])
}

#' Check batch variable
#'
#' @param x batch variable.
#' @param min_n_batch min. number of batches (for MMUPHin functions to run).
#'
#' @return if no errors then the batch variables (factorized if not already)
#' @keywords internal
check_batch <- function(x, min_n_batch = 2) {
  # First ensure batch variable is a factor
  if(!is.factor(x)) {
    warning("Batch variable is not a factor as provided and will be converted ",
            "to one.")
    x <- as.factor(x)
  }

  # Check min number of batches is satisfied
  if(nlevels(x) < min_n_batch)
    stop("Must have at least ", min_n_batch, " batches!")

  return(x)
}

#' Match user-specified control parameters with default, and modify if needed
#'
#' @param default list of default control parameters
#' @param control list of user-provided control parameters
#'
#' @return list of control parameters, set to user provided values if specified
#' and default other wise
#' @keywords internal
match_control <- function(default, control) {
  if (missing(control))
    control <- list()

  control_pos <- pmatch(names(control), names(default))
  default[c(na.omit(control_pos))] <- control[!is.na(control_pos)]

  return(default)
}

#' Set pseudo count for an abundance matrix. Pseudo count is currently set to
#' half of minimum non-zero values
#'
#' @param features feature-by-sample matrix of abundances (proportions or
#' counts).
#'
#' @return the pseudo count
#' @keywords internal
set_pseudo <- function(features) {
  type_features <- check_feature_abd(features)
  if(all(features == 0))
    stop("All feature abundances are zero!")

  min(setdiff(features, 0)) / 2
}

#' Normalize feature abundance table (modified from Maaslin2)
#'
#' @param features feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param normalization normalization method.
#' @param pseudo_count pseudo count to be added to feature_abd.
#'
#' @return normalized abundance table.
#' @keywords internal
normalize_features <- function(features,
                               normalization = "NONE",
                               pseudo_count = 0) {
  if(any(features < 0))
    stop("Feature table must be non-negative for normalization!")
  features <- features + pseudo_count
  if (normalization=='TSS')
    features <- apply(features, 2, TSS)
  if (normalization=='NONE')
    features <- features
  return(features)
}

#' Transform feature abunadnce table (modified from Maaslin2)
#'
#' @param features feature-by-sample matrix of abundances (proportions or
#' counts).
#' @param transform transformation method.
#' @param pseudo_count pseudo count to be added to feature_abd..
#'
#' @return transformed abundance table.
#' @keywords internal
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

#' TSS normalization (modified from Maaslin2)
#'
#' @param x vector of abundance to be normalized.
#'
#' @return normalized vector of abundance.
#' @keywords internal
TSS <- function(x) {
  if(all(x == 0)) return(x)
  return(x / sum(x))
}

#' AST transformation (modified from Maaslin2 and is different)
#'
#' @param x vector of abundance to be transformed.
#'
#' @return transformed vector of abundance.
#' @keywords internal
AST<-function(x){
  return(asin(sqrt(x)))
}

#' LOG transformation (modified from Maaslin2 and is different)
#'
#' @param x vector of abundance to be transformed.
#'
#' @return transformed vector of abundance.
#' @keywords internal
LOG<-function(x){
  return(log2(x))
}

#' Fill in artificial row/column names to a matrix or data frame, if they are
#' missing
#'
#' @param x matrix or data frame
#' @param row_prefix prefix for the artificial row names
#' @param col_prefix prefix for the artificial column names
#'
#' @return x but with the missing dimension names filled in
#' @keywords internal
fill_dimnames <- function(x, row_prefix, col_prefix) {
  if(!(is.matrix(x) | is.data.frame(x)))
    stop("x must either be a matrix or a data frame!")
  if(missing(row_prefix) | missing(col_prefix))
    stop("Row/column prefixes must be specified!")
  if(is.null(rownames(x)))
    rownames(x) <- paste0(row_prefix, seq_len(nrow(x)))
  if(is.null(colnames(x)))
    colnames(x) <- paste0(col_prefix, seq_len(ncol(x)))

  return(x)
}

#' Utility for shorter names
#' Useful when plotting per-feature figures where feature names could be cutoff
#'
#' @param x vector of names
#' @param cutoff number of maximum string length before start cutting off the
#' middle
#'
#' @return vector of new names with .. replacing the middle part if name is
#' longer than cutoff
#' @keywords internal
shorten_name <- function(x, cutoff = 3, replacement = "..") {
  if(anyDuplicated(x))
    stop("x should be unique character strings!")
  x_sub <- x
  length_x <- nchar(x_sub, type = "c")
  ind_change <- length_x > cutoff*2 + nchar(replacement, type = "c")
  x_sub[ind_change] <- paste0(substr(x_sub[ind_change],
                                     start = 1, stop = cutoff),
                              replacement,
                              substr(x_sub[ind_change],
                                     start = length_x[ind_change] - cutoff + 1,
                                     stop = length_x[ind_change]))
  # To avoid the case where shortening makes names indistinguishable
  if(anyDuplicated(x_sub))
    return(x)
  else
    return(x_sub)
}

#' Utility for catching warning/error messages
#'
#' @param expr an expression to run that can generate potential errors/warnings
#'
#' @return a list, capturing both the return value of the expression, as well
#' as generated erros/warnings (\code{NULL} if no errors/warnings)
#' @keywords internal
catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), 
                  error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 

#' Utility for checking options
#'
#' @param x the specified value
#' @param x_name name of the specified value
#' @param options allowed options
#'
#' @return error if \code{x} is not in \code{options}. Otherwise returns
#' \code{x}.
#' @keywords internal
check_options <- function(x, x_name, options) {
  if(!(x %in% options))
    stop(x_name, "can only be one of ", paste(options, collapse = ", "))
  return(x)
}

#' Utility for checking continuous options
#'
#' @param x the specified numeric value
#' @param x_name name of the specified value
#' @param range allowed range
#'
#' @return error if \code{x} is not within \code{range} (boundaries 
#' excluded). Otherwise returns \code{x}.
#' @keywords internal
check_options_continuous <- function(x, x_name, range) {
  if(!(x > range[1] & x < range[2]))
    stop(x_name, "can only be between ", paste(range, collapse = ", "))
  return(x)
}

#' Utility for checking pseudo count
#'
#' @param x the specified pseudo count
#'
#' @return error if pseudo count is smaller than zero. Otherwise returns 
#' \code{x}.
#' @keywords internal
check_pseudo_count <- function(x) {
  if(x < 0)
    stop("pseudo_count must be non-negative")
  return(x)
}