#' Parse Principal Stratification Formula
#'
#' Parses a formula in 2SLS style: Y ~ X1 + X2 | Z | W
#'
#' @param formula A formula with three parts separated by |
#' @param data A data.frame containing the variables
#'
#' @return A list with components X, X_raw, Y, Z, W
#'
#' @keywords internal
parse_psbart_formula <- function(formula, data) {
  if (is.null(data)) {
    stop("data argument is required when using formula interface")
  }

  # Convert to Formula object for multi-part handling
  f <- Formula::Formula(formula)

  # Check formula structure
  n_rhs <- length(f)[2]
  if (n_rhs != 3) {
    stop("Formula must have 3 right-hand-side parts: Y ~ covariates | Z | W\n",
         "  Example: outcome ~ x1 + x2 + x3 | assignment | uptake")
  }

  # Extract outcome (Y)
  Y <- stats::model.response(
    stats::model.frame(f, data = data, lhs = 1, rhs = 0)
  )

  # Extract covariates (X) - first RHS part
  X_formula <- stats::formula(f, lhs = 0, rhs = 1)
  X_frame <- stats::model.frame(X_formula, data = data)
  X <- stats::model.matrix(X_formula, data = X_frame)
  # Remove intercept if present
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  # Extract Z (instrument/assignment) - second RHS part
  Z_formula <- stats::formula(f, lhs = 0, rhs = 2)
  Z_frame <- stats::model.frame(Z_formula, data = data)
  Z <- stats::model.matrix(Z_formula, data = Z_frame)
  if ("(Intercept)" %in% colnames(Z)) {
    Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  }
  if (ncol(Z) != 1) {
    stop("Z (instrument) must be a single variable")
  }
  Z <- as.vector(Z)

  # Extract W (uptake) - third RHS part
  W_formula <- stats::formula(f, lhs = 0, rhs = 3)
  W_frame <- stats::model.frame(W_formula, data = data)
  W <- stats::model.matrix(W_formula, data = W_frame)
  if ("(Intercept)" %in% colnames(W)) {
    W <- W[, colnames(W) != "(Intercept)", drop = FALSE]
  }
  if (ncol(W) != 1) {
    stop("W (uptake) must be a single variable")
  }
  W <- as.vector(W)

  list(X = X, X_raw = as.data.frame(X_frame), Y = Y, Z = Z, W = W)
}


#' @keywords internal
normalize_raw_covariates <- function(X) {
  if (is.data.frame(X)) {
    return(as.data.frame(X, stringsAsFactors = FALSE))
  }
  if (is.matrix(X)) {
    return(as.data.frame(X, stringsAsFactors = FALSE))
  }
  stop("X must be a matrix or data.frame")
}


#' @keywords internal
get_fit_covariates <- function(princebart_fit, type = c("model", "raw")) {
  type <- match.arg(type)

  if (is.null(princebart_fit$data)) {
    return(NULL)
  }

  if (type == "model") {
    return(princebart_fit$data$X_model)
  }

  princebart_fit$data$X_raw
}


#' @keywords internal
validate_fit_covariate_contract <- function(princebart_fit, require_raw = TRUE) {
  if (is.null(princebart_fit$data)) {
    stop("Fit object is missing data storage; please refit with current prince_BART().")
  }
  if (is.null(princebart_fit$data$X_model)) {
    stop("Fit object is missing data$X_model; please refit with current prince_BART().")
  }
  if (require_raw && is.null(princebart_fit$data$X_raw)) {
    stop("Fit object is missing data$X_raw; please refit with current prince_BART().")
  }
  if (is.null(princebart_fit$scaling) ||
      is.null(princebart_fit$scaling$center) ||
      is.null(princebart_fit$scaling$scale)) {
    stop("Fit object is missing scaling metadata; please refit with current prince_BART().")
  }
}


#' Internal Binary BART Wrapper
#'
#' A custom wrapper for dbarts that properly handles binary outcomes,
#' including edge cases where a stratum may have only 0s or only 1s.
#'
#' This addresses a limitation in dbarts::dbarts where binary detection
#' only triggers when there are exactly 2 unique values (0 and 1).
#' In principal stratification, subsets may sometimes contain only one
#' unique outcome value, which would cause dbarts to treat it as continuous
#' and fail when estimating sigma.
#'
#' Based on the fix from princeB.r lines 679-739.
#'
#' @param x Predictor matrix.
#' @param y Response vector (binary 0/1).
#' @param test Test predictor matrix.
#' @param subset Logical vector for subsetting training data.
#' @param offset Offset values for training data.
#' @param control A dbartsControl object.
#' @param k Prior hyperparameter for node prior (default 2).
#'
#' @return A dbartsSampler object configured for binary probit BART.
#'
#' @keywords internal
dbarts_binary <- function(
  x, y, test = NULL, subset = NULL, offset = NULL,
  control = dbarts::dbartsControl(),
  k = 2
) {
  # Get internal dbarts functions
  validateArgumentsInEnvironment <- utils::getFromNamespace(
    "validateArgumentsInEnvironment", "dbarts"
  )
  redirectCall <- utils::getFromNamespace("redirectCall", "dbarts")
  addCallArgument <- utils::getFromNamespace("addCallArgument", "dbarts")
  quoteInNamespace <- utils::getFromNamespace("quoteInNamespace", "dbarts")
  parsePriors <- utils::getFromNamespace("parsePriors", "dbarts")
  setDefaultsFromFormals <-
    utils::getFromNamespace("setDefaultsFromFormals", "dbarts")
  normal <- utils::getFromNamespace("normal", "dbarts")
  fixed <- utils::getFromNamespace("fixed", "dbarts")
  cgm <- utils::getFromNamespace("cgm", "dbarts")

  # Build the call as dbarts expects
  matchedCall <- match.call()

  # Remap our arguments to dbarts::dbarts argument names
  matchedCall[[1L]] <- quote(dbarts::dbarts)
  names(matchedCall)[names(matchedCall) == "x"] <- "formula"
  names(matchedCall)[names(matchedCall) == "y"] <- "data"

  evalEnv <- parent.frame(1L)

  # Set up control
  control@verbose <- FALSE

  # Create dbartsData - this is the standard data processing
  dataCall <- matchedCall
  dataCall[[1L]] <- quoteInNamespace("dbartsData")
  # Remove arguments that dbartsData doesn't accept
  dataCall$control <- NULL
  dataCall$k <- NULL

  data <- eval(dataCall, evalEnv, getNamespace("dbarts"))
  data@n.cuts <- rep_len(attr(control, "n.cuts"), ncol(data@x))
  # Use NA_real_ as in original; dbarts handles this for binary
  data@sigma <- NA_real_
  attr(control, "n.cuts") <- NULL


  # KEY FIX: Enhanced binary detection that handles edge cases

  # Original dbarts only checks:
  #' length(uniqueResponses) == 2 && all(sort(uniqueResponses) == c(0, 1))
  #' We add: length(uniqueResponses) == 1 && 
  #' all(uniqueResponses == 0 | uniqueResponses == 1)
  uniqueResponses <- unique(data@y)
  if (length(uniqueResponses) == 2 && all(sort(uniqueResponses) == c(0, 1))) {
    control@binary <- TRUE
  }
  # Edge case: stratum has only 0s or only 1s
  #'(can happen in principal stratification)
  if (length(uniqueResponses) == 1 &&
        all(uniqueResponses == 0 | uniqueResponses == 1)) {
    control@binary <- TRUE
  }

  # For binary, clear offsets if all zero
  if (control@binary && !is.null(data@offset) && all(data@offset == 0)) {
    data@offset <- NULL
  }
  if (control@binary &&
      !is.null(data@offset.test) &&
      all(data@offset.test == 0)
  ) {
    data@offset.test <- NULL
  }

  # Get prior objects directly
  # normal(k) returns a list with $node.prior and $node.hyperprior
  # cgm() returns a tree prior object directly
  # fixed(1) returns a resid prior object directly
  node_prior_result <- normal(k)
  tree_prior <- cgm()
  resid_prior <- fixed(1)

  # Create the model with binary-appropriate settings
  model <- methods::new(
    "dbartsModel"
    , tree_prior
    , node_prior_result$node.prior
    , node_prior_result$node.hyperprior
    , resid_prior
    , proposal.probs =
      c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)
    , node.scale = 3.0  # 3.0 for binary, 0.5 for continuous
  )

  # Create and return the sampler
  result <- methods::new("dbartsSampler", control, model, data)
  result
}

#' @keywords internal
resolve_and_validate_uptake <- function(W, uptake_type) {
  W <- as.numeric(W)
  if (any(is.na(W))) {
    stop("W must not contain missing values")
  }

  unique_vals <- unique(W)
  is_binary <- all(unique_vals %in% c(0, 1))

  if (uptake_type == "auto") {
    resolved_type <- if (is_binary) "binary" else "ordinal"
  } else {
    resolved_type <- uptake_type
  }

  if (resolved_type == "binary") {
    if (!is_binary) {
      stop("For uptake_type = 'binary', W must be binary (0/1)")
    }
    return(list(W = W, uptake_type = resolved_type))
  }

  # Ordinal/count uptake: non-negative integers
  is_integer_like <- abs(W - round(W)) <= sqrt(.Machine$double.eps)
  if (!all(is_integer_like) || any(W < 0)) {
    stop("For uptake_type = 'ordinal', W must be non-negative integer values")
  }

  list(W = as.numeric(round(W)), uptake_type = resolved_type)
}

#' @keywords internal
validate_and_prepare_X <- function(X) {
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix or data.frame")
  }
  scale(X)
}

#' @keywords internal
validate_binary <- function(x, name) {
  x <- as.numeric(x)
  unique_vals <- unique(x[!is.na(x)])
  if (!all(unique_vals %in% c(0, 1))) {
    stop(name, " must be binary (0/1)")
  }
  x
}

#' @keywords internal
validate_propensity <- function(propensity, n) {
  if (length(propensity) != n) {
    stop("propensity must have length equal to number of observations")
  }
  if (any(propensity <= 0 | propensity >= 1, na.rm = TRUE)) {
    stop("propensity scores must be strictly between 0 and 1")
  }
}
