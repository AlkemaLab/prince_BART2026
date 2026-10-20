#' Set Up a Reusable Principal Stratification Simulation
#'
#' Fix coefficients and settings for [ps_sim_data()]. Random coefficients use
#' the caller's RNG; `seed` controls subsequent data generation only.
#'
#' @param n Positive integer sample size; default 100.
#' @param nx1,nx2 Numbers of continuous and binary predictors; defaults 1 each.
#' @param w_levels Number of uptake levels: 2 for binary, greater than 2 for count.
#' @param x1_transform Function applied to each continuous predictor value;
#'   defaults to [base::identity()]. Must return one finite number per value.
#' @param y_type Outcome distribution: `"binary"` (default) or `"continuous"`.
#' @param sigma_y Positive residual standard deviation for continuous outcomes.
#' @param x1_mean,x1_var Normal means and nonnegative variances for continuous
#'   predictors. Each accepts a scalar or a vector of length `nx1`.
#' @param x2_p Bernoulli probabilities: a scalar or a vector of length `nx2`.
#'   Predictors are drawn independently.
#' @param coefficients `"random"` for N(0,1) draws, a finite scalar filling all
#'   blocks, or a named list. Missing blocks are drawn randomly with a warning.
#' @param seed Data-generation seed, or `NULL` to use the caller's RNG.
#'   Seeded data generation restores the caller's RNG state.
#' @param stratum_weights `"default"` or a nonnegative `w_levels` by `w_levels`
#'   matrix with positive total weight. Rows index W(0), columns W(1).
#'   Zeros exclude strata; positive entries set relative weights.
#'
#' @details
#' Coefficient blocks (K = `w_levels`) are `g_intercept`, `g_u`, `y_intercept`,
#' `y_u`, `y_w` (K by K); `g_x1`, `y_x1` (`nx1` by K by K); `g_x2`, `y_x2`
#' (`nx2` by K by K); `z_intercept` (scalar); `z_x1` (length `nx1`); and `z_x2`
#' (length `nx2`). Stratum axes are W(0), then W(1).
#' Default support permits W(1) >= W(0) for binary uptake and W(1) <= W(0)
#' for count uptake, matching the package's fitting conventions.
#'
#' @return A list of validated settings, resolved coefficient blocks, and the
#'   resolved stratum weight matrix. Edit settings and reuse with [ps_sim_data()].
#' @examples
#' set.seed(10) # Controls coefficient randomization.
#' setup <- ps_sim_setup(n = 100, x1_transform = sin, seed = 1)
#' first <- ps_sim_data(setup)
#' setup$seed <- 2
#' second <- ps_sim_data(setup) # Same coefficients, new data.
#' @export
ps_sim_setup <- function(n = 100L, nx1 = 1L, nx2 = 1L, w_levels = 2L,
                         x1_transform = identity, y_type = "binary",
                         sigma_y = 1, x1_mean = 0, x1_var = 1, x2_p = 0.5,
                         coefficients = "random", seed = NULL,
                         stratum_weights = "default") {
  setup <- list(n = n, nx1 = nx1, nx2 = nx2, w_levels = w_levels,
    x1_transform = x1_transform, y_type = y_type, sigma_y = sigma_y,
    x1_mean = x1_mean, x1_var = x1_var, x2_p = x2_p,
    coefficients = coefficients, seed = seed, stratum_weights = stratum_weights)
  .ps_sim_validate_setup(setup, initialize = TRUE)
}

#' Generate Data From a Simulation Setup
#'
#' Generate observed and latent data using fixed coefficients from [ps_sim_setup()].
#' @param setup A setup returned by [ps_sim_setup()], or `NULL` to create one.
#'   Supplied setups are checked again for consistency.
#' @param ... Arguments to [ps_sim_setup()] when `setup = NULL`. When supplying
#'   a setup, edit its fields directly instead of passing additional arguments.
#'
#' @details
#' Stratum probabilities are proportional to stratum weights times exponentiated
#' linear scores. Instrument assignment uses a logistic model of transformed
#' covariates. Outcomes depend on covariates, an unobserved N(0,1) variable,
#' stratum, and received uptake, using a logistic link for binary outcomes or
#' an identity link for continuous outcomes.
#'
#' @return A list containing `data` (observed covariates, Z, W, Y), `alldata`
#'   (also transformed predictors, G, U, W0, W1, propensity, mean_y0, mean_y1),
#'   and `setup` (the input or generated setup). `mean_y0` and `mean_y1` are
#'   conditional outcome means under assignment Z = 0 and Z = 1.
#' @examples
#' sim <- ps_sim_data(n = 100, coefficients = 0, seed = 42)
#' head(sim$data)
#' shifted <- sim$setup
#' shifted$x1_mean <- 2
#' target <- ps_sim_data(shifted)
#' @export
ps_sim_data <- function(setup = NULL, ...) {
  if (is.null(setup)) {
    setup <- ps_sim_setup(...)
  } else {
    if (length(list(...))) {
      stop("Supply either setup or setup arguments, not both.", call. = FALSE)
    }
    setup <- .ps_sim_validate_setup(setup)
  }
  n <- setup$n
  nx1 <- setup$nx1
  nx2 <- setup$nx2
  w_levels <- setup$w_levels
  y_type <- setup$y_type
  sigma_y <- setup$sigma_y
  coefficients <- setup$coefficients
  stratum_weights <- setup$stratum_weights

  .ps_sim_with_seed(setup$seed, {
    x1 <- matrix(stats::rnorm(n * nx1,
      mean = rep(setup$x1_mean, each = n),
      sd = rep(sqrt(setup$x1_var), each = n)), nrow = n, ncol = nx1)
    x2 <- matrix(stats::rbinom(n * nx2, 1L,
      prob = rep(setup$x2_p, each = n)), nrow = n, ncol = nx2)
    u <- stats::rnorm(n)
    h <- matrix(vapply(as.vector(x1), function(x) {
      value <- setup$x1_transform(x)
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        stop("x1_transform must return one finite number per value.", call. = FALSE)
      }
      as.numeric(value)
    }, numeric(1L)), nrow = n, ncol = nx1)

    # Array flattening follows R's column order: W(0) varies fastest.
    score <- .ps_sim_predictor(h, x2, u, coefficients, "g", w_levels)
    allowed <- which(as.vector(stratum_weights) > 0)
    log_score <- sweep(score[, allowed, drop = FALSE], 2L,
                       log(stratum_weights[allowed]), "+")
    mass <- exp(log_score - apply(log_score, 1L, max))
    probability <- mass / rowSums(mass)
    g <- allowed[vapply(seq_len(n), function(i) {
      sample.int(length(allowed), size = 1L, prob = probability[i, ])
    }, integer(1L))]
    w0 <- (g - 1L) %% w_levels
    w1 <- (g - 1L) %/% w_levels

    z_score <- as.vector(coefficients$z_intercept + h %*% coefficients$z_x1 +
                         x2 %*% coefficients$z_x2)
    if (any(!is.finite(z_score))) {
      stop("Instrument predictors are nonfinite; reduce coefficient magnitudes.",
           call. = FALSE)
    }
    propensity <- stats::plogis(z_score)
    z <- stats::rbinom(n, 1L, propensity)
    w <- ifelse(z == 0L, w0, w1)

    baseline <- .ps_sim_predictor(h, x2, u, coefficients, "y", w_levels)
    baseline <- baseline[cbind(seq_len(n), g)]
    mean_y0 <- baseline + coefficients$y_w[g] * w0
    mean_y1 <- baseline + coefficients$y_w[g] * w1
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
    list(data = data, alldata = alldata, setup = setup)
  })
}

.ps_sim_validate_setup <- function(setup, initialize = FALSE) {
  required <- names(formals(ps_sim_setup))
  nm <- names(setup)
  if (!is.list(setup) || is.data.frame(setup) || is.null(nm) ||
      anyNA(nm) || anyDuplicated(nm) || !setequal(nm, required)) {
    stop("setup must contain exactly the named settings returned by ps_sim_setup().",
         call. = FALSE)
  }
  setup$n <- .ps_sim_integer(setup$n, "n", 1)
  setup$nx1 <- .ps_sim_integer(setup$nx1, "nx1", 0)
  setup$nx2 <- .ps_sim_integer(setup$nx2, "nx2", 0)
  setup$w_levels <- .ps_sim_integer(setup$w_levels, "w_levels", 2)
  if (!is.function(setup$x1_transform)) {
    stop("x1_transform must be a function.", call. = FALSE)
  }
  setup$y_type <- match.arg(setup$y_type, c("binary", "continuous"))
  if (!is.numeric(setup$sigma_y) || length(setup$sigma_y) != 1L ||
      !is.finite(setup$sigma_y) || setup$sigma_y <= 0) {
    stop("sigma_y must be a positive finite number.", call. = FALSE)
  }
  setup$x1_mean <- .ps_sim_parameter(setup$x1_mean, "x1_mean", setup$nx1)
  setup$x1_var <- .ps_sim_parameter(setup$x1_var, "x1_var", setup$nx1, lower = 0)
  setup$x2_p <- .ps_sim_parameter(setup$x2_p, "x2_p", setup$nx2, lower = 0, upper = 1)
  if (!is.null(setup$seed)) setup$seed <- .ps_sim_integer(setup$seed, "seed", 0)
  setup$stratum_weights <- .ps_sim_support(setup$stratum_weights, setup$w_levels)
  setup$coefficients <- .ps_sim_coefficients(setup$nx1, setup$nx2, setup$w_levels,
    setup$coefficients, initialize = initialize)
  setup
}

.ps_sim_parameter <- function(value, name, count, lower = -Inf, upper = Inf) {
  if (!is.numeric(value) || !is.null(dim(value)) ||
      !length(value) %in% c(1L, count) || any(!is.finite(value)) ||
      any(value < lower | value > upper)) {
    stop(name, " must be finite numeric of length 1 or ", count,
         " with values in [", lower, ", ", upper, "].", call. = FALSE)
  }
  rep(as.numeric(value), length.out = count)
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

.ps_sim_coefficients <- function(nx1, nx2, k, coefficients, initialize) {
  shapes <- list(
    g_intercept = c(k, k), g_x1 = c(nx1, k, k),
    g_x2 = c(nx2, k, k), g_u = c(k, k),
    z_intercept = 1L, z_x1 = nx1, z_x2 = nx2,
    y_intercept = c(k, k), y_x1 = c(nx1, k, k),
    y_x2 = c(nx2, k, k), y_u = c(k, k), y_w = c(k, k)
  )
  supplied <- is.list(coefficients) && !is.data.frame(coefficients)
  random <- identical(coefficients, "random")
  constant <- is.numeric(coefficients) && length(coefficients) == 1L &&
    is.null(dim(coefficients)) && is.finite(coefficients)
  if (!initialize && !supplied) {
    stop("setup$coefficients must be a complete named list.", call. = FALSE)
  }
  if (!supplied && !random && !constant) {
    stop("coefficients must be 'random', a finite numeric scalar, or a named list.",
         call. = FALSE)
  }
  if (supplied) {
    nm <- names(coefficients)
    if (length(coefficients) && (is.null(nm) || anyNA(nm) || any(nm == "") ||
                            anyDuplicated(nm) || any(!nm %in% names(shapes)))) {
      stop("coefficients must have unique names from the documented coefficient blocks.",
           call. = FALSE)
    }
    for (name in nm) {
      value <- coefficients[[name]]
      shape <- shapes[[name]]
      valid_shape <- if (length(shape) == 1L) {
        is.null(dim(value)) && length(value) == shape
      } else {
        identical(dim(value), as.integer(shape))
      }
      if (!is.numeric(value) || any(!is.finite(value)) || !valid_shape) {
        stop("coefficients$", name, " must be finite numeric with ",
             if (length(shape) == 1L) "length " else "dimensions ",
             paste(shape, collapse = " x "), ".", call. = FALSE)
      }
    }
    missing <- setdiff(names(shapes), nm)
    if (length(missing)) {
      if (!initialize) {
        stop("Missing coefficient blocks in setup: ", paste(missing, collapse = ", "),
             call. = FALSE)
      }
      warning("Missing coefficient blocks: ", paste(missing, collapse = ", "),
              ". Drawing them from N(0, 1).", call. = FALSE)
    }
  }
  out <- lapply(names(shapes), function(name) {
    if (supplied && name %in% names(coefficients)) return(coefficients[[name]])
    shape <- shapes[[name]]
    value <- if (constant) rep(as.numeric(coefficients), prod(shape)) else
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
  if (identical(stratum_weights, "default")) {
    support <- matrix(1, k, k)
    support <- if (k == 2L) 1 * (row(support) <= col(support)) else
      1 * (row(support) >= col(support))
    return(support)
  }
  if (!is.matrix(stratum_weights) || !is.numeric(stratum_weights) ||
      !identical(dim(stratum_weights), c(k, k)) ||
      any(!is.finite(stratum_weights)) || any(stratum_weights < 0) ||
      !any(stratum_weights > 0)) {
    stop("stratum_weights must be 'default' or a finite nonnegative ",
         k, " x ", k, " matrix with at least one positive entry.", call. = FALSE)
  }
  stratum_weights
}

.ps_sim_predictor <- function(h, x2, u, coefficients, prefix, k) {
  block <- function(suffix) coefficients[[paste0(prefix, "_", suffix)]]
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
