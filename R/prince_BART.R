#' @keywords internal
combine_chain_array <- function(lst, var_names) {
  out <- abind::abind(lst, along = 4)
  out <- aperm(out, c(1, 4, 3, 2))
  dimnames(out) <- list(
    iteration = NULL,
    chain = NULL,
    variable = var_names,
    unit = NULL
  )
  out
}

#' @keywords internal
combine_chain_trees <- function(chain_results, keep_trees) {
  if (!keep_trees) {
    return(NULL)
  }

  list_tree <- lapply(chain_results, function(x) x$trees)
  do.call(rbind, Map(function(df, id) {
    df$chain <- id
    df
  }, list_tree, seq_along(list_tree)))
}


#' @keywords internal
.run_psbart_binary_chains <- function(
  X,
  Y,
  Z,
  W,
  n_warmup,
  n_samples,
  n_chains,
  keep_trees,
  k,
  n_trees,
  n_initial,
  verbose,
  workers
) {
  res0 <- future.apply::future_lapply(
    seq_len(n_chains),
    function(chain_id) {
      .fit_psbart_binary(
        X = X,
        Y = Y,
        Z = Z,
        W = W,
        n_warmup = n_warmup,
        n_samples = n_samples,
        save_trees = keep_trees,
        k = k,
        n_trees = n_trees,
        n_initial = n_initial,
        verbose = FALSE
      )
    },
    future.seed = TRUE
  )

  list(
    trees = combine_chain_trees(res0, keep_trees)
    , imp = combine_chain_array(
      lapply(res0, function(x) x$imputed), c("nt", "at")
    )
    , probs = combine_chain_array(
      lapply(res0, function(x) x$probs)
      , c("p_a", "p_n", "m_y0c", "m_y1c", "m_y0n", "m_y1a")
    )
    , check = NULL
  )
}

#' @keywords internal
.run_psbart_ordinal_chains <- function(
  X,
  Y,
  Z,
  W,
  n_warmup,
  n_samples,
  n_chains,
  keep_trees,
  k,
  n_trees,
  lambda_z1,
  lambda_z0,
  init_w_poisson_lambda,
  n_thin,
  n_threads,
  monotonicity,
  rho,
  verbose,
  workers
) {
  data_for_fit <- data.frame(
    Y = Y,
    Z = Z,
    W = W,
    X,
    check.names = FALSE
  )

  res0 <- future.apply::future_lapply(
    seq_len(n_chains),
    function(chain_id) {
      .fit_psbart_ordinal(
        data = data_for_fit,
        n_warmup = n_warmup,
        n_samples = n_samples,
        k = k,
        n_trees = n_trees,
        lambda_z1 = lambda_z1,
        lambda_z0 = lambda_z0,
        init_w_poisson_lambda = init_w_poisson_lambda,
        n_thin = n_thin,
        n_threads = n_threads,
        monotonicity = monotonicity,
        rho = rho,
        save_trees = keep_trees,
        verbose = FALSE
      )
    },
    future.seed = TRUE
  )

  list(
    trees = combine_chain_trees(res0, keep_trees),
    imp = combine_chain_array(lapply(res0, function(x) x$imputed), c("w0", "w1")),
    probs = combine_chain_array(lapply(res0, function(x) x$probs), c("m_y0", "m_y1")),
    check = combine_chain_array(lapply(res0, function(x) x$check), c("check_w0", "check_w1"))
  )
}

#' Principal Stratification using BART
#'
#' Fits a Bayesian principal stratification model using Bayesian Additive
#' Regression Trees (BART) for causal inference with endogenous treatments.
#' The method is designed for instrumental variable or encouragement designs
#' with noncompliance and targets causal effects for compliers.
#'
#' @param formula A formula following 2SLS convention:
#'   \code{Y ~ X1 + X2 | Z | W}, where Y is the outcome, X1 + X2 are covariates,
#'   Z is the instrument (treatment assignment), and W is treatment uptake.
#' @param data A data.frame containing the variables in the formula.
#' @param X A matrix or data.frame of covariates. Used when \code{formula = NULL}.
#' @param Y A binary outcome vector (0/1). Used when \code{formula = NULL}.
#' @param Z A binary instrument or treatment assignment vector (0/1).
#' @param W Treatment uptake/received vector.
#' @param propensity Optional pre-computed instrument propensity scores
#'   \eqn{e = P(Z \mid X)}. If NULL (default), propensity scores are estimated
#'   internally using BART.
#' @param instrument_overlap Optional bounds for instrument propensity trimming
#'   to enforce overlap. If provided as a length-2 vector (e.g., \code{c(0.1, 0.9)}),
#'   observations with estimated propensity scores outside this range are excluded.
#'   Default is NULL (no trimming).
#' @param n_chains Number of parallel MCMC chains (default: 4).
#' @param n_warmup Number of warmup iterations per chain (default: 1000).
#' @param n_samples Number of posterior samples per chain (default: 1000).
#' @param keep_trees Logical; save fitted BART tree structures for downstream
#'   prediction and generalization (default: FALSE).
#' @param k Prior hyperparameter controlling node shrinkage in BART (default: 2).
#' @param n_trees Number of trees in each BART ensemble (default: 200).
#' @param workers Number of parallel workers. If NULL (default), uses all
#'   available cores up to \code{n_chains}. Set to 1 for sequential execution.
#' @param uptake_type Character scalar controlling uptake model:
#'   \code{"auto"} (default) chooses binary if \code{W in \{0,1\}} and ordinal otherwise;
#'   \code{"binary"} enforces binary uptake model;
#'   \code{"ordinal"} enforces ordinal/count uptake model.
#' @param rho Numeric in (-1, 1); residual correlation between the latent
#'   propensity scores for \eqn{W(1)} and \eqn{W(0)} in the bivariate ordinal
#'   sampler. Only used when \code{uptake_type = "ordinal"}. Default is
#'   \code{0} (conditionally independent potential treatments given X).
#'   Positive values encode a tendency for units with high baseline treatment
#'   to also have high treated-condition treatment; because this assumption is
#'   not testable from the observed data, varying \code{rho} is a natural
#'   sensitivity analysis.
#' @param verbose Logical; print progress messages (default: FALSE).
#'
#' @details
#' Shared preprocessing is applied in all modes: formula parsing, validation,
#' covariate scaling, propensity estimation for the binary instrument Z,
#' optional overlap trimming, and propensity augmentation of X.
#' After preprocessing, model fitting dispatches to either binary or ordinal
#' PS-BART chain runners based on \code{uptake_type}.
#'
#' For fits created from a formula or data.frame input, the returned object
#' stores both a model-space covariate matrix and, when available, a raw
#' covariate data.frame. The model-space representation is used internally for
#' fitting and prediction; the raw representation is retained for downstream
#' tasks that benefit from original factor/ordered-factor classes, such as
#' \code{segment_heterogeneity()}.
#'
#' The model jointly estimates:
#' \itemize{
#'   \item Principal stratum membership probabilities
#'     (\eqn{P(\text{complier} \mid X)}, \eqn{P(\text{never-taker} \mid X)},
#'      \eqn{P(\text{always-taker} \mid X)}),
#'   \item Potential outcome regressions within strata, including
#'     \eqn{E[Y(0) \mid \text{complier}, X]} and
#'     \eqn{E[Y(1) \mid \text{complier}, X]}.
#' }
#'
#' The primary identified causal estimand is the conditional average treatment
#' effect among compliers, \eqn{\mathrm{CATE}_C(x)}. Mixed (sample-based) averages
#' of these conditional effects, including LATE-like estimands, can be obtained
#' using \code{summary()} and \code{coef()}.
#'
#' Identification relies on standard instrumental variable assumptions:
#' conditional independence of the instrument given X, exclusion restriction,
#' and monotonicity (no defiers).
#'
#' Parallel computation is handled via
#' \code{future.apply::future_lapply}. Before calling this function, set your
#' preferred parallel backend:
#'
#' \preformatted{
#' # For Unix/Mac (forked processes)
#' future::plan(future::multicore)
#'
#' # For Windows or any platform
#' future::plan(future::multisession)
#' }
#'
#' @return A list of class \code{"prince_bart"} containing posterior draws from
#'   all chains, including:
#'   \itemize{
#'     \item \code{imp}: Imputed latent quantities
#'       (iteration x chain x variable x unit).
#'     \item \code{probs}: Posterior draws of probabilities/outcome means.
#'     \item \code{check}: Ordinal diagnostic array (ordinal mode only).
#'     \item \code{trees}: Fitted BART trees (if \code{keep_trees = TRUE}).
#'     \item \code{data}: Stored input data, including \code{X_model}
#'       (processed covariates used by the fitted model), \code{X_raw}
#'       (raw covariates when available), \code{X} (compatibility alias to
#'       \code{X_model}), instrument propensity scores, and outcomes.
#'   }
#'
#' @examples
#' \dontrun{
#' library(future)
#' plan(multisession, workers = 4)
#'
#' fit <- prince_BART(
#'   Y ~ x1_1 + x2_1 | Z | W,
#'   data = ps_simulate_data(seed = 1)$data,
#'   uptake_type = "auto",
#'   n_chains = 4,
#'   n_warmup = 1000,
#'   n_samples = 1500
#' )
#'
#' summary(fit)
#' }
#'
#' @seealso \code{\link{segment_heterogeneity}},
#'   \code{\link{general_BART}}
#'
#' @export
prince_BART <- function(
  formula = NULL,
  data = NULL,
  X = NULL,
  Y = NULL,
  Z = NULL,
  W = NULL,
  propensity = NULL,
  instrument_overlap = NULL,
  n_warmup = 1000L,
  n_samples = 1000L,
  n_chains = 4L,
  keep_trees = FALSE,
  k = 2,
  n_trees = 200L,
  workers = NULL,
  uptake_type = c("auto", "binary", "ordinal"),
  rho = 0,
  verbose = FALSE
) {

  uptake_type <- match.arg(uptake_type)

  # Hardcode n_initial = 0 (MoM offsets not clearly helpful)
  n_initial <- 0L

  # Ordinal defaults (internal)
  ordinal_defaults <- list(
    lambda_z1 = 1,
    lambda_z0 = 1,
    init_w_poisson_lambda = 3,
    n_thin = 1L,
    n_threads = 1L,
    monotonicity = TRUE
  )

  # Handle formula interface for both binary and ordinal
  X_raw <- NULL
  if (!is.null(formula)) {
    parsed <- parse_psbart_formula(formula, data)
    X <- parsed$X
    X_raw <- parsed$X_raw
    Y <- parsed$Y
    Z <- parsed$Z
    W <- parsed$W
  } else if (!is.null(X)) {
    X_raw <- normalize_raw_covariates(X)
  }

  # Input validation
  if (is.null(X) || is.null(Y) || is.null(Z) || is.null(W)) {
    stop("Must provide either a formula + data, or X, Y, Z, W directly")
  }
  if (is.null(X_raw)) {
    stop("Could not construct raw covariates (X_raw); please provide valid X/data inputs.")
  }

  # Common preprocessing before dispatch
  X <- validate_and_prepare_X(X) #df > m > scaled(m)
  Y <- validate_binary(Y, "Y")
  Z <- validate_binary(Z, "Z")

  uptake_info <- resolve_and_validate_uptake(W, uptake_type)
  W <- uptake_info$W
  uptake_type <- uptake_info$uptake_type

  n <- nrow(X)
  if (length(Y) != n || length(Z) != n || length(W) != n) {
    stop("X, Y, Z, W must all have the same number of observations")
  }

  # Store scaling attributes before adding e
  scaled_center <- attr(X, "scaled:center")
  scaled_scale <- attr(X, "scaled:scale")

  # Compute propensity for the instrument Z in both modes
  if (is.null(propensity)) {
    if (verbose) message("Computing propensity scores for Z...")
    e <- dbarts::bart2(X, Z, verbose = FALSE) |>
      stats::fitted() |>
      stats::qnorm()
  } else {
    validate_propensity(propensity, n)
    e <- stats::qnorm(propensity)
  }

  # Enforce instrument overlap by trimming extreme propensity scores
  if (!is.null(instrument_overlap)) {
    if (length(instrument_overlap) != 2) {
      stop("instrument_overlap must be a length-2 vector, e.g., c(0.1, 0.9)")
    }
    e_prob <- stats::pnorm(e)
    keep <- e_prob >= instrument_overlap[1] & e_prob <= instrument_overlap[2]
    n_trimmed <- sum(!keep)
    if (verbose) {
      message(
        "Trimming ", n_trimmed, " observations (",
        round(100 * n_trimmed / n, 1), "%) outside propensity range [",
        instrument_overlap[1], ", ", instrument_overlap[2], "]"
      )
    }
    if (sum(keep) < 10) {
      stop("Too few observations remain after instrument overlap trimming")
    }
    X <- X[keep, , drop = FALSE]
    # make sure trim raw covariates (alternative save keep)
    X_raw <- X_raw[keep, , drop = FALSE]
    Y <- Y[keep]
    Z <- Z[keep]
    W <- W[keep]
    e <- e[keep]
    n <- sum(keep)
  }

  # Append propensity to X in both modes
  X <- cbind(X, e = e)

  # Determine number of workers
  if (is.null(workers)) {
    workers <- min(parallel::detectCores(), n_chains)
  }

  # Set up future plan if not already configured
  old_plan <- future::plan()
  if (inherits(old_plan, "sequential") && workers > 1) {
    if (verbose) message(
      "Setting up multisession plan with ", workers, " workers"
    )
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  if (verbose) message(
    "Running ", n_chains, " chains (", uptake_type, " uptake)..."
  )

  # Dispatch to modality-specific chain runner
  chain_results <- if (uptake_type == "binary") {
    .run_psbart_binary_chains(
      X = X,
      Y = Y,
      Z = Z,
      W = W,
      n_warmup = n_warmup,
      n_samples = n_samples,
      n_chains = n_chains,
      keep_trees = keep_trees,
      k = k,
      n_trees = n_trees,
      n_initial = n_initial,
      verbose = verbose,
      workers = workers
    )
  } else {
    .run_psbart_ordinal_chains(
      X = X,
      Y = Y,
      Z = Z,
      W = W,
      n_warmup = n_warmup,
      n_samples = n_samples,
      n_chains = n_chains,
      keep_trees = keep_trees,
      k = k,
      n_trees = n_trees,
      lambda_z1 = ordinal_defaults$lambda_z1,
      lambda_z0 = ordinal_defaults$lambda_z0,
      init_w_poisson_lambda = ordinal_defaults$init_w_poisson_lambda,
      n_thin = ordinal_defaults$n_thin,
      n_threads = ordinal_defaults$n_threads,
      monotonicity = ordinal_defaults$monotonicity,
      rho = rho,
      verbose = verbose,
      workers = workers
    )
  }

  # Unified output assembly
  res <- list(
    trees = chain_results$trees,
    imp = chain_results$imp,
    probs = chain_results$probs
  )

  if (!is.null(chain_results$check)) {
    res$check <- chain_results$check
  }

  # Store data reference with UNSCALED X (for general_BART compatibility)
  # X is currently scaled with e appended; unscale covariates before storing
  X_unscaled <- X[, -ncol(X), drop = FALSE]
  if (!is.null(scaled_center) && !is.null(scaled_scale)) {
    for (j in seq_len(ncol(X_unscaled))) {
      v <- colnames(X_unscaled)[j]
      if (v %in% names(scaled_center) && v %in% names(scaled_scale)) {
        X_unscaled[, j] <- X_unscaled[, j] * scaled_scale[v] + scaled_center[v]
      }
    }
  }

  # Add back propensity (never scaled)
  X_unscaled <- cbind(X_unscaled, e = X[, "e"])

  res$data <- list(
    X_model = X_unscaled,
    X_raw = X_raw,
    Y = Y,
    Z = Z,
    W = W,
    e = X[, "e"]
  )
  res$scaling <- list(center = scaled_center, scale = scaled_scale)
  res$uptake_type <- uptake_type
  res$call <- match.call()

  class(res) <- if (uptake_type == "binary") {
    c("prince_bart", "prince_bart_binary")
  } else {
    c("prince_bart", "prince_bart_ordinal")
  }

  if (verbose) message("Done.")
  res
}
