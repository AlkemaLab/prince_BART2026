#' Generate Coefficients for Principal Stratification Simulations
#'
#' Generate a reusable set of coefficients for [ps_simulate_data()]. Coefficients
#' are drawn once per call, so the same returned list can be used to generate
#' independent datasets under the same model.
#'
#' @param nx1 Nonnegative integer. Number of independent standard normal
#'   covariates, named `x1_1`, `x1_2`, and so on. Default is 1.
#' @param nx2 Nonnegative integer. Number of independent Bernoulli(0.5)
#'   covariates, named `x2_1`, `x2_2`, and so on. Default is 1. Both covariate
#'   counts may be zero.
#' @param w_levels Integer at least 2. Number of treatment uptake levels,
#'   labeled `0`, ..., `w_levels - 1`. This controls uptake, not the outcome.
#' @param weights Either `"random"` (the default), a finite numeric scalar,
#'   or a named list of coefficient blocks. `"random"` draws each coefficient
#'   independently from N(0, 1). A scalar fills every coefficient with that
#'   constant. A list supplies individual blocks using the names and dimensions
#'   below; omitted blocks are drawn independently from N(0, 1), with a warning.
#'   Unknown names, duplicate names, and incorrect dimensions are errors.
#' @param seed `NULL` (default) or a nonnegative integer seed. With `NULL`,
#'   random draws use and advance the caller's random-number state. With a
#'   supplied seed, results are reproducible and the caller's state is restored
#'   on exit, including when an error occurs.
#'
#' @details
#' Write K = `w_levels`. In every stratum-indexed block, the two stratum axes
#'   are W(0), then W(1), in the order `0`, ..., `K - 1`.
#'
#' The returned coefficient blocks are:
#' \describe{
#'   \item{`g_intercept`, `g_u`}{K by K matrices for stratum score intercepts
#'     and the unobserved covariate U.}
#'   \item{`g_x1`, `g_x2`}{Arrays with dimensions (nx1, K, K) and (nx2, K, K)
#'     for transformed continuous and binary covariates in stratum scores.}
#'   \item{`z_intercept`}{Scalar intercept for instrument assignment.}
#'   \item{`z_x1`, `z_x2`}{Vectors of lengths nx1 and nx2 for instrument
#'     assignment. U does not enter this model.}
#'   \item{`y_intercept`, `y_u`, `y_w`}{K by K matrices for outcome intercepts,
#'     U coefficients, and uptake coefficients, respectively.}
#'   \item{`y_x1`, `y_x2`}{Arrays with dimensions (nx1, K, K) and (nx2, K, K)
#'     for transformed continuous and binary covariates in the outcome model.}
#' }
#'
#' Zero covariate counts produce zero-length vectors and arrays with a zero
#' first dimension. Coefficients are returned even for strata excluded by a
#' subsequent simulation's `stratum_weights`.
#'
#' @return A named list containing all coefficient blocks listed above. The
#'   list can be passed unchanged as `weights` to [ps_simulate_data()] with
#'   matching `nx1`, `nx2`, and `w_levels`.
#' @examples
#' weights <- ps_simulate_weights(seed = 42)
#' first <- ps_simulate_data(n = 100, weights = weights, seed = 1)
#' second <- ps_simulate_data(n = 100, weights = weights, seed = 2)
#' identical(first$weights, second$weights)
#'
#' # Specify a scenario with equal allowed-stratum probabilities,
#' # randomized assignment, and a positive uptake effect on log odds.
#' weights <- ps_simulate_weights(weights = 0)
#' weights$y_w[,] <- 0.8
#' sim <- ps_simulate_data(weights = weights, seed = 3)
#' table(sim$alldata$G)
#' @seealso [ps_simulate_data()]
#' @export
ps_simulate_weights <- function(nx1 = 1L, nx2 = 1L, w_levels = 2L,
                                weights = "random", seed = NULL) {
  nx1 <- .ps_sim_integer(nx1, "nx1", 0)
  nx2 <- .ps_sim_integer(nx2, "nx2", 0)
  w_levels <- .ps_sim_integer(w_levels, "w_levels", 2)
  .ps_sim_with_seed(seed, .ps_sim_weights(nx1, nx2, w_levels, weights))
}

#' Simulate Data for Principal Stratification
#'
#' Generate observed data and latent quantities under an instrumental-variable
#' model with binary or ordinal uptake. Continuous covariates can have additive
#' nonlinear effects. The output supports examples, tests, and repeated
#' simulations with fixed coefficients.
#'
#' @param n Positive integer. Number of observations. Default is 100.
#' @inheritParams ps_simulate_weights
#' @param y_type Outcome distribution: `"binary"` (default) or `"continuous"`.
#'   Binary outcomes use a logistic link; continuous outcomes use a normal
#'   distribution. The current [prince_BART()] fitter accepts binary outcomes
#'   only, including when uptake is ordinal.
#' @param x1_transform Elementwise transformation of continuous covariates:
#'   `"linear"` (identity, default), `"sin"` (sine in radians), `"relu"`
#'   (maximum of zero and x), or `"sigmoid"` (inverse logit). The same
#'   transformation is used in all three models. Returned covariates are on
#'   their original scale. No interaction features are added.
#' @param stratum_weights `"default"` (default), `"full"`, or a finite,
#'   nonnegative `w_levels` by `w_levels` matrix with at least one positive
#'   entry. Rows index W(0) and columns index W(1), in order `0`, ...,
#'   `w_levels - 1`; matrix entries are matched by position. Zeros exclude
#'   strata; positive entries multiply their exponentiated scores before
#'   normalization. `"full"` uses a matrix of ones. To match the current
#'   package, `"default"` permits W(1) >= W(0) for binary uptake and
#'   W(1) <= W(0) for ordinal uptake. These are opposite directions of
#'   monotonicity. Custom matrices can describe departures from the fitting
#'   model's monotonicity assumption.
#' @param sigma_y Positive finite scalar. Residual standard deviation for
#'   continuous outcomes; default is 1. Validated but unused for binary
#'   outcomes, whose variance is determined by their success probabilities.
#'
#' @details
#' For each person, draw independent standard normal continuous covariates and
#' U, and independent Bernoulli(0.5) binary covariates. U is unobserved in the
#' analysis dataset. Write H for the transformed continuous covariates and B
#' for the binary covariates. For stratum g = (a, b), let
#' \deqn{\eta_g = \alpha_g + H^T\beta_g + B^T\delta_g + U\gamma_g.}
#' With M equal to the resolved `stratum_weights` matrix, draw G according to
#' \deqn{P(G=(a,b)\mid X,U) =
#'   \frac{M_{ab}\exp(\eta_{ab})}{\sum_{r,s}M_{rs}\exp(\eta_{rs})}.}
#' M sets support and relative weights; its entries are not generally the
#' marginal stratum probabilities.
#'
#' Independently of U and G conditional on X, draw Z with probability
#' \deqn{e(X) = \operatorname{logit}^{-1}
#'   (\alpha_Z + H^T\beta_Z + B^T\delta_Z).}
#' Observed uptake is W = (1 - Z)a + Zb. Define the outcome predictor at uptake w
#' in stratum g by
#' \deqn{m(w,g) = \alpha^Y_g + H^T\beta^Y_g + B^T\delta^Y_g
#'   + U\gamma^Y_g + w\tau_g.}
#' Draw the observed outcome as Bernoulli(inverse-logit(m(W,G))) or
#' normal(m(W,G), sigma_y squared). The outcome predictor has no direct Z
#' term, so the construction satisfies exclusion. U can affect both uptake
#' stratum and outcomes without confounding instrument assignment given X.
#'
#' Random coefficients are drawn once per dataset. For repeated datasets under
#' one model, generate coefficients with [ps_simulate_weights()] and pass the
#' result to each call. Large coefficients may produce extreme probabilities;
#' no rejection sampling is used to ensure every stratum or instrument arm is
#' observed. Small samples may also omit some uptake levels. Use an explicit
#' `uptake_type` when fitting an intended ordinal simulation.
#'
#' `mean_y0` and `mean_y1` are conditional outcome means under Z = 0 and Z = 1,
#' given the realized X, U, and G. Their difference is an instrument contrast.
#' In the default ordinal model Z reduces uptake. These conditional means are
#' not effects averaged over U given X and G, or population-level estimands.
#'
#' Both covariate counts may be zero for intercept-and-U simulations. Such
#' datasets are useful for checking the generator; support for fitting them
#' depends on the fitting function's covariate requirements.
#'
#' @return A named list with four components:
#' \describe{
#'   \item{data}{A data.frame containing raw covariates `x1_1`, ... and `x2_1`,
#'     ..., followed by `Z`, `W`, and `Y`. It contains only observed variables.}
#'   \item{alldata}{The same data plus transformed continuous covariates
#'     `x1_t_1`, ... (when present), `G` (a factor with labels `(a,b)`),
#'     `U`, `W0`, `W1`, `propensity`, `mean_y0`, and `mean_y1`. Factor levels
#'     include every allowed stratum, including those absent from this sample.}
#'   \item{weights}{The complete named list of coefficients actually used,
#'     in the format returned by [ps_simulate_weights()].}
#'   \item{stratum_weights}{The resolved `w_levels` by `w_levels` matrix used
#'     to set allowed strata and their relative weights. Rows index `W0` and
#'     columns index `W1`; these entries are not marginal stratum probabilities.}
#' }
#' @examples
#' sim <- ps_simulate_data(seed = 42)
#' head(sim$data)
#' head(sim$alldata)
#'
#' # Reuse coefficients with nonlinear additive covariate effects.
#' weights <- ps_simulate_weights(nx1 = 2, nx2 = 0, seed = 10)
#' sim <- ps_simulate_data(
#'   nx1 = 2, nx2 = 0, x1_transform = "sin", weights = weights, seed = 11
#' )
#'
#' # Ordinal uptake; the default has W(1) <= W(0).
#' ordinal <- ps_simulate_data(w_levels = 4, seed = 12)
#' all(ordinal$alldata$W1 <= ordinal$alldata$W0)
#'
#' # Continuous outcomes and no observed covariates.
#' continuous <- ps_simulate_data(
#'   nx1 = 0, nx2 = 0, y_type = "continuous", sigma_y = 0.5, seed = 13
#' )
#'
#' # Specify exact allowed strata: rows are W(0), columns are W(1).
#' compliers <- matrix(0, 2, 2)
#' compliers[1, 2] <- 1
#' sim <- ps_simulate_data(stratum_weights = compliers, seed = 14)
#' all(sim$data$W == sim$data$Z)
#'
#' \dontrun{
#' # A small fit for illustrating the interface; assess convergence separately.
#' sim <- ps_simulate_data(n = 300, seed = 42)
#' fit <- prince_BART(
#'   Y ~ x1_1 + x2_1 | Z | W, data = sim$data,
#'   propensity = sim$alldata$propensity, uptake_type = "binary",
#'   n_chains = 2, workers = 1, n_warmup = 100, n_samples = 100
#' )
#' }
#' @seealso [ps_simulate_weights()], [prince_BART()]
#' @export
ps_simulate_data <- function(n = 100L, nx1 = 1L, nx2 = 1L, w_levels = 2L,
                             y_type = c("binary", "continuous"),
                             seed = NULL,
                             x1_transform = c("linear", "sin", "relu", "sigmoid"),
                             weights = "random", stratum_weights = "default",
                             sigma_y = 1) {
  n <- .ps_sim_integer(n, "n", 1)
  nx1 <- .ps_sim_integer(nx1, "nx1", 0)
  nx2 <- .ps_sim_integer(nx2, "nx2", 0)
  w_levels <- .ps_sim_integer(w_levels, "w_levels", 2)
  y_type <- match.arg(y_type)
  x1_transform <- match.arg(x1_transform)
  if (!is.numeric(sigma_y) || length(sigma_y) != 1L ||
      !is.finite(sigma_y) || sigma_y <= 0) {
    stop("sigma_y must be a positive finite number.", call. = FALSE)
  }
  support <- .ps_sim_support(stratum_weights, w_levels)

  .ps_sim_with_seed(seed, {
    weights <- .ps_sim_weights(nx1, nx2, w_levels, weights)
    x1 <- matrix(stats::rnorm(n * nx1), nrow = n, ncol = nx1)
    x2 <- matrix(stats::rbinom(n * nx2, 1L, 0.5), nrow = n, ncol = nx2)
    u <- stats::rnorm(n)
    h <- switch(x1_transform,
      linear = x1, sin = sin(x1), relu = pmax(x1, 0),
      sigmoid = stats::plogis(x1)
    )

    # Array flattening follows R's column order: W(0) varies fastest.
    score <- .ps_sim_predictor(h, x2, u, weights, "g", w_levels)
    allowed <- which(as.vector(support) > 0)
    log_score <- sweep(score[, allowed, drop = FALSE], 2L,
                       log(support[allowed]), "+")
    mass <- exp(log_score - apply(log_score, 1L, max))
    probability <- mass / rowSums(mass)
    g <- allowed[vapply(seq_len(n), function(i) {
      sample.int(length(allowed), size = 1L, prob = probability[i, ])
    }, integer(1L))]
    w0 <- (g - 1L) %% w_levels
    w1 <- (g - 1L) %/% w_levels

    z_score <- as.vector(weights$z_intercept + h %*% weights$z_x1 +
                         x2 %*% weights$z_x2)
    if (any(!is.finite(z_score))) {
      stop("Instrument predictors are nonfinite; reduce coefficient magnitudes.",
           call. = FALSE)
    }
    propensity <- stats::plogis(z_score)
    z <- stats::rbinom(n, 1L, propensity)
    w <- ifelse(z == 0L, w0, w1)

    baseline <- .ps_sim_predictor(h, x2, u, weights, "y", w_levels)
    baseline <- baseline[cbind(seq_len(n), g)]
    mean_y0 <- baseline + weights$y_w[g] * w0
    mean_y1 <- baseline + weights$y_w[g] * w1
    if (any(!is.finite(c(mean_y0, mean_y1)))) {
      stop("Outcome predictors are nonfinite; reduce coefficient magnitudes.",
           call. = FALSE)
    }
    if (y_type == "binary") {
      mean_y0 <- stats::plogis(mean_y0)
      mean_y1 <- stats::plogis(mean_y1)
    }
    observed_mean <- ifelse(z == 0L, mean_y0, mean_y1)
    y <- if (y_type == "binary") {
      stats::rbinom(n, 1L, observed_mean)
    } else {
      stats::rnorm(n, observed_mean, sigma_y)
    }
    data <- as.data.frame(cbind(x1, x2))
    names(data) <- c(if (nx1 > 0L) paste0("x1_", seq_len(nx1)),
                     if (nx2 > 0L) paste0("x2_", seq_len(nx2)))
    data$Z <- z
    data$W <- w
    data$Y <- y
    labels <- paste0("(", (allowed - 1L) %% w_levels, ",",
                     (allowed - 1L) %/% w_levels, ")")
    alldata <- data
    for (j in seq_len(nx1)) alldata[[paste0("x1_t_", j)]] <- h[, j]
    alldata$G <- factor(g, levels = allowed, labels = labels)
    alldata$U <- u
    alldata$W0 <- w0
    alldata$W1 <- w1
    alldata$propensity <- propensity
    alldata$mean_y0 <- mean_y0
    alldata$mean_y1 <- mean_y1
    list(data = data, alldata = alldata, weights = weights,
         stratum_weights = support)
  })
}

.ps_sim_integer <- function(x, name, lower) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > .Machine$integer.max || x != floor(x)) {
    stop(name, " must be an integer >= ", lower,
         " and <= .Machine$integer.max.", call. = FALSE)
  }
  as.integer(x)
}

.ps_sim_with_seed <- function(seed, expr) {
  if (!is.null(seed)) {
    seed <- .ps_sim_integer(seed, "seed", 0)
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(list = ".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  force(expr)
}

.ps_sim_weights <- function(nx1, nx2, k, weights) {
  shapes <- list(
    g_intercept = c(k, k), g_x1 = c(nx1, k, k),
    g_x2 = c(nx2, k, k), g_u = c(k, k),
    z_intercept = 1L, z_x1 = nx1, z_x2 = nx2,
    y_intercept = c(k, k), y_x1 = c(nx1, k, k),
    y_x2 = c(nx2, k, k), y_u = c(k, k), y_w = c(k, k)
  )
  supplied <- is.list(weights) && !is.data.frame(weights)
  random <- identical(weights, "random")
  constant <- is.numeric(weights) && length(weights) == 1L &&
    is.null(dim(weights)) && is.finite(weights)
  if (!supplied && !random && !constant) {
    stop("weights must be 'random', a finite numeric scalar, or a named list.",
         call. = FALSE)
  }
  if (supplied) {
    nm <- names(weights)
    if (length(weights) && (is.null(nm) || anyNA(nm) || any(nm == "") ||
                            anyDuplicated(nm) || any(!nm %in% names(shapes)))) {
      stop("weights must have unique names from the documented coefficient blocks.",
           call. = FALSE)
    }
    for (name in nm) {
      value <- weights[[name]]
      shape <- shapes[[name]]
      valid_shape <- if (length(shape) == 1L) {
        is.null(dim(value)) && length(value) == shape
      } else {
        identical(dim(value), as.integer(shape))
      }
      if (!is.numeric(value) || any(!is.finite(value)) || !valid_shape) {
        stop("weights$", name, " must be finite numeric with ",
             if (length(shape) == 1L) "length " else "dimensions ",
             paste(shape, collapse = " x "), ".", call. = FALSE)
      }
    }
    missing <- setdiff(names(shapes), nm)
    if (length(missing)) {
      warning("Missing coefficient blocks: ", paste(missing, collapse = ", "),
              ". Drawing them from N(0, 1).", call. = FALSE)
    }
  }
  out <- lapply(names(shapes), function(name) {
    if (supplied && name %in% names(weights)) return(weights[[name]])
    shape <- shapes[[name]]
    value <- if (constant) rep(as.numeric(weights), prod(shape)) else
      stats::rnorm(prod(shape))
    if (length(shape) == 1L) return(value)
    levels <- as.character(seq_len(k) - 1L)
    axes <- list(W0 = levels, W1 = levels)
    if (length(shape) == 3L) axes <- c(list(covariate = NULL), axes)
    array(value, dim = shape, dimnames = axes)
  })
  stats::setNames(out, names(shapes))
}

.ps_sim_support <- function(stratum_weights, k) {
  if (is.character(stratum_weights) && length(stratum_weights) == 1L &&
      !is.na(stratum_weights)) {
    mode <- match.arg(stratum_weights, c("default", "full"))
    support <- matrix(1, k, k)
    if (mode == "default") {
      support <- if (k == 2L) 1 * (row(support) <= col(support)) else
        1 * (row(support) >= col(support))
    }
    return(support)
  }
  if (!is.matrix(stratum_weights) || !is.numeric(stratum_weights) ||
      !identical(dim(stratum_weights), c(k, k)) ||
      any(!is.finite(stratum_weights)) || any(stratum_weights < 0) ||
      !any(stratum_weights > 0)) {
    stop("stratum_weights must be 'default', 'full', or a finite nonnegative ",
         k, " x ", k, " matrix with at least one positive entry.", call. = FALSE)
  }
  stratum_weights
}

.ps_sim_predictor <- function(h, x2, u, weights, prefix, k) {
  block <- function(suffix) weights[[paste0(prefix, "_", suffix)]]
  score <- h %*% matrix(block("x1"), nrow = ncol(h), ncol = k * k) +
    x2 %*% matrix(block("x2"), nrow = ncol(x2), ncol = k * k) +
    outer(u, as.vector(block("u")))
  score <- sweep(score, 2L, as.vector(block("intercept")), "+")
  if (any(!is.finite(score))) {
    stop("Simulation predictors are nonfinite; reduce coefficient magnitudes.",
         call. = FALSE)
  }
  score
}
