# Weight-shift sensitivity analysis for target-population effects.
sensitivity_analysis <- function(tau, subpop, psu, weights, gamma, seed = NULL) {

  # Sample subset of draws for computational efficiency
  n_draws <- ncol(tau)
  n_sample <- min(n_draws, 100)

  if (!is.null(seed)) set.seed(seed)
  sample_idx <- sample(n_draws, n_sample)
  tau_sample <- tau[, sample_idx, drop = FALSE]

  # Find lower bound (inf)
  lower <- compute_shift_pate(
    tau = tau_sample,
    subpop = subpop,
    psu = psu,
    weights = weights,
    gamma = gamma,
    inf = TRUE,
    seed = seed
  )

  # Find upper bound (sup)
  upper <- compute_shift_pate(
    tau = tau_sample,
    subpop = subpop,
    psu = psu,
    weights = weights,
    gamma = gamma,
    inf = FALSE,
    seed = seed
  )

  list(
    gamma = gamma,
    lower = lower,
    upper = upper
  )
}
compute_shift_pate <- function(tau, subpop, psu, weights, gamma, inf, seed = NULL) {
  tau <- as.matrix(tau)
  tau[!subpop, ] <- NA
  n_draws <- ncol(tau)

  n_psu <- length(unique(psu))
  p_c <- tapply(weights, psu, mean)
  n_c <- tapply(!is.na(psu), psu, sum)

  if (!is.null(seed)) set.seed(seed)
  bb_list <- lapply(seq_len(n_draws), function(r) bayesian_bootstrap(n_psu))

  pate_draws <- sapply(seq_len(n_draws), function(k) {
    y <- tau[, k]
    shift_wts <- find_shift_weights(inf = inf, gamma = gamma, y = y, wts = weights)
    psu_tau <- tapply(y * shift_wts, psu, mean, na.rm = TRUE)
    stats::weighted.mean(psu_tau, bb_list[[k]] * p_c * n_c, na.rm = TRUE)
  })

  list(
    estimate = mean(pate_draws, na.rm = TRUE),
    sd = stats::sd(pate_draws, na.rm = TRUE),
    ci = stats::quantile(pate_draws, c(0.025, 0.975), na.rm = TRUE)
  )
}
#' Find optimal weight shift for sensitivity bounds
#' @keywords internal
find_shift_weights <- function(inf, gamma, y, wts = NULL) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    warning("CVXR not available, returning uniform weights")
    return(rep(1, length(y)))
  }

  if (is.null(wts)) wts <- rep(1, length(y))

  new_wts <- rep(0, length(y))
  s <- !is.na(y)

  wts_s <- as.numeric(wts)[s]
  sum_wts <- sum(wts_s)
  wts_s <- wts_s / sum(wts_s)
  y_s <- y[s]
  n <- sum(s)

  r <- CVXR::Variable(n)

  if (inf) {
    objective <- CVXR::Minimize(sum(y_s * wts_s * r))
  } else {
    objective <- CVXR::Maximize(sum(y_s * wts_s * r))
  }

  constraints <- list(
    sum(wts_s * r) == 1,
    r <= gamma,
    r >= 1 / gamma
  )

  problem <- CVXR::Problem(objective, constraints = constraints)
  result <- solve(problem, solver = "ECOS")

  if (result$status != "optimal") {
    warning("Optimization did not converge, returning uniform weights")
    return(rep(1, length(y)))
  }

  r_val <- round(result$getValue(r), 4)
  new_wts[s] <- wts_s * r_val * sum_wts
  new_wts
}

#' Sensitivity Analysis for Unmeasured Confounding
#'
#' Compute sensitivity bounds for the PATE estimate under potential
#' unmeasured confounding using weight-shift optimization.
#'
#' @param object A `general_pate` object from [general_BART()]
#' @param gamma Numeric sensitivity parameter > 1. Represents the maximum
#'   ratio by which observation weights can be shifted. Larger values allow
#'   more severe confounding.
#' @param n_sample Number of posterior draws to use for sensitivity analysis.
#'   Default is 100.
#' @param verbose Logical; print progress messages. Default is FALSE.
#'
#' @return A list with components:
#' \describe{
#'   \item{gamma}{The sensitivity parameter used}
#'   \item{lower}{Lower bound estimate with CI}
#'   \item{upper}{Upper bound estimate with CI}
#' }
#'
#' @details
#' The sensitivity analysis follows the weight-shift framework where
#' observation weights can be multiplied by a factor between 1/gamma and
#' gamma. The lower and upper bounds represent the extremes of the PATE
#' estimate under this weight perturbation.
#'
#' Requires the CVXR package for optimization.
#'
#' @seealso [general_BART()], [general_BART_overlap()]
#'
#' @examples
#' \dontrun{
#' fit <- general_BART(princebart_fit, newdata, ...)
#' sens <- general_BART_transportability(fit, gamma = 1.5)
#' }
#'
#' @export
general_BART_transportability <- function(
    object,
    gamma,
    n_sample = 100,
    verbose = FALSE
) {
  if (!inherits(object, "general_pate")) {
    stop("object must be a 'general_pate' object from general_BART()")
  }

  if (!is.numeric(gamma) || gamma <= 1) {
    stop("gamma must be numeric and > 1")
  }

  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' is required for sensitivity analysis")
  }

  if (verbose) message("Running sensitivity analysis (gamma = ", gamma, ")...")

  # Get tau from object
  tau <- object$tau
  n_new <- dim(tau)[1]
  n_iter <- dim(tau)[2]
  n_ch <- dim(tau)[3]
  tau_2d <- matrix(tau, nrow = n_new, ncol = n_iter * n_ch)

  # Sample subset of draws for computational efficiency
  n_draws <- ncol(tau_2d)
  n_sample <- min(n_draws, n_sample)

  set.seed(object$seed)
  sample_idx <- sample(n_draws, n_sample)
  tau_sample <- tau_2d[, sample_idx, drop = FALSE]

  if (verbose) message("  Computing lower bound...")

  # Find lower bound (inf)
  lower <- compute_shift_pate(
    tau = tau_sample,
    subpop = object$subpop,
    psu = object$psu,
    weights = object$weights,
    gamma = gamma,
    inf = TRUE,
    seed = object$seed
  )

  if (verbose) message("  Computing upper bound...")

  # Find upper bound (sup)
  upper <- compute_shift_pate(
    tau = tau_sample,
    subpop = object$subpop,
    psu = object$psu,
    weights = object$weights,
    gamma = gamma,
    inf = FALSE,
    seed = object$seed
  )

  result <- list(
    gamma = gamma,
    lower = lower,
    upper = upper
  )

  if (verbose) message("Done.")
  result
}
