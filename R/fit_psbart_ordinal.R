#' Fit Principal Stratification BART for Ordinal Uptake (Single Chain)
#'
#' Internal function that fits one MCMC chain of the ordinal principal
#' stratification BART model. Users should typically call
#' \code{\link{prince_BART}}.
#'
#' @param data A data.frame/data.table containing columns \code{Y}, \code{Z},
#'   \code{W}, and covariates.
#' @param n_warmup Number of warmup iterations.
#' @param n_samples Number of posterior samples to retain.
#' @param k BART node prior hyperparameter.
#' @param n_trees Number of trees per BART ensemble.
#' @param lambda_z1 Box-Cox parameter for latent \eqn{z_1}.
#' @param lambda_z0 Box-Cox parameter for latent \eqn{z_0}.
#' @param init_w_poisson_lambda Poisson rate used to jitter initial
#'   \code{w0}/\code{w1} assignments.
#' @param n_thin Thinning interval for retained posterior samples.
#' @param n_threads Number of threads used by dbarts.
#' @param monotonicity Logical; if TRUE enforce \eqn{w_1 \le w_0} when
#'   constructing candidate strata.
#' @param rho Correlation used in bivariate latent sampling.
#' @param save_trees Logical; if TRUE include sampled trees.
#' @param save_trees_interval Save trees every this many retained samples.
#' @param verbose Logical; print progress updates.
#' @param ... Unused.
#'
#' @return A list of iteration-by-unit matrices \code{w0}, \code{w1},
#'   \code{m_y0}, \code{m_y1}, \code{check_w0}, and \code{check_w1}, plus a
#'   \code{trees} table when requested.
#'
#' @keywords internal
.fit_psbart_ordinal <- function(
  data = NULL,
  n_warmup = 20L,
  n_samples = 20L,
  k = 2,
  n_trees = 200L,
  lambda_z1 = 1,
  lambda_z0 = 1,
  init_w_poisson_lambda = 3,
  n_thin = 1L,
  n_threads = 4L,
  monotonicity = TRUE,
  rho = 0,
  save_trees = TRUE,
  save_trees_interval = 1L,
  verbose = FALSE,
  ...
) {
  validate_ordinal_inputs(
    data = data,
    n_warmup = n_warmup,
    n_samples = n_samples,
    n_trees = n_trees,
    n_thin = n_thin,
    n_threads = n_threads,
    save_trees_interval = save_trees_interval
  )

  prepared <- prepare_ordinal_data(data)
  X <- prepared$X
  Y <- prepared$Y
  Z <- prepared$Z
  W <- prepared$W

  n <- nrow(X)
  n_total <- n_warmup + (n_samples * n_thin)

  storage <- initialize_ordinal_storage(n_samples = n_samples, n = n)
  trees <- list()

  w_states <- initialize_w_states(
    W = W,
    Z = Z,
    n = n,
    init_w_poisson_lambda = init_w_poisson_lambda
  )
  w0 <- w_states$w0
  w1 <- w_states$w1
  w_max <- w_states$w_max

  control <- dbartsControl(
    updateState = TRUE,
    keepTrees = save_trees,
    verbose = FALSE,
    n.burn = 0L,
    n.samples = 1L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = n_trees,
    n.threads = n_threads
  )

  latent <- initialize_latent_states(
    w0 = w0,
    w1 = w1,
    lambda_z0 = lambda_z0,
    lambda_z1 = lambda_z1
  )
  z0 <- latent$z0
  z1 <- latent$z1
  sig_z0 <- latent$sig_z0
  sig_z1 <- latent$sig_z1

  x_plus <- cbind(X, w0 = w0, w1 = w1)
  w_grid <- build_w_grid(w_max = w_max, monotonicity = monotonicity)

  test_data <- prepare_test_matrices(
    X = X,
    Z = Z,
    w0 = w0,
    w1 = w1,
    w_grid = w_grid
  )
  X_w_grid <- test_data$X_w_grid
  id <- test_data$id

  samplers <- initialize_ordinal_samplers(
    X = X,
    Y = Y,
    Z = Z,
    z0 = z0,
    z1 = z1,
    sig_z0 = sig_z0,
    sig_z1 = sig_z1,
    x_plus = x_plus,
    X_w_grid = X_w_grid,
    control = control,
    k = k
  )

  set_initial_sampler_data(
    samplers = samplers,
    X = X,
    Y = Y,
    Z = Z,
    z0 = z0,
    z1 = z1,
    x_plus = x_plus,
    X_w_grid = X_w_grid
  )

  for (s in list(samplers$z0, samplers$z1, samplers$y0, samplers$y1)){
    s$sampleTreesFromPrior()
  }

  offset <- 0L
  for (i in seq_len(n_total)) {
    latent_draw <- sample_latent_z_states(
      w0 = w0,
      w1 = w1,
      w_max = w_max,
      mu_z0 = samplers$state$mu_z0,
      mu_z1 = samplers$state$mu_z1,
      sig_z0 = samplers$state$sig_z0,
      sig_z1 = samplers$state$sig_z1,
      lambda_z0 = lambda_z0,
      lambda_z1 = lambda_z1,
      rho = rho,
      prev_z0 = z0,
      prev_z1 = z1
    )
    z0 <- latent_draw$z0
    z1 <- latent_draw$z1

    samplers$z1$setResponse(z1)
    samplers$z0$setResponse(z0)

    sample_z1 <- samplers$z1$run(updateState = TRUE)
    sample_z0 <- samplers$z0$run(updateState = TRUE)

    samplers$state$mu_z1 <- sample_z1$test
    samplers$state$sig_z1 <- sample_z1$sigma
    samplers$state$mu_z0 <- sample_z0$test
    samplers$state$sig_z0 <- sample_z0$sigma

    check_w1 <- round_floor(
      ibc(rnorm(n, samplers$state$mu_z1, samplers$state$sig_z1), lambda_z1),
      y_max = w_max
    )
    check_w0 <- round_floor(
      ibc(rnorm(n, samplers$state$mu_z0, samplers$state$sig_z0), lambda_z0),
      y_max = w_max
    )

    sample_y1 <- samplers$y1$run()
    sample_y0 <- samplers$y0$run()
    mean_y1 <- pnorm(sample_y1$test)
    mean_y0 <- pnorm(sample_y0$test)

    strata_probs_dt <- build_strata_posterior_dt(
      X_w_grid = X_w_grid,
      id = id,
      Y = Y,
      Z = Z,
      mean_y0 = mean_y0,
      mean_y1 = mean_y1,
      mu_z0 = samplers$state$mu_z0,
      mu_z1 = samplers$state$mu_z1,
      sig_z0 = samplers$state$sig_z0,
      sig_z1 = samplers$state$sig_z1
    )

    sampled_strata <- sample_strata_membership(strata_probs_dt)
    w0 <- sampled_strata$w0
    w1 <- sampled_strata$w1

    x_plus <- cbind(X, w0 = w0, w1 = w1)
    update_ordinal_outcome_samplers(
      samplers = samplers,
      x_plus = x_plus,
      Y = Y,
      Z = Z,
      X_w_grid = X_w_grid
    )

    if (i > n_warmup && i %% n_thin == 0) {
      offset <- offset + 1L
      storage$samples_w0[offset, ] <- w0
      storage$samples_w1[offset, ] <- w1
      storage$mean_y0[offset, ] <- sampled_strata$my0
      storage$mean_y1[offset, ] <- sampled_strata$my1
      storage$check_w0[offset, ] <- check_w0
      storage$check_w1[offset, ] <- check_w1

      if (save_trees && offset %% save_trees_interval == 0) {
        trees[[offset]] <- extract_ordinal_trees(
          sampler_z0 = samplers$z0,
          sampler_z1 = samplers$z1,
          sampler_y0 = samplers$y0,
          sampler_y1 = samplers$y1,
          iteration = offset
        )
      }
    }

    if (verbose && i %% 100 == 0) {
      message("Iteration ", i, " of ", n_total)
    }
  }

  if (save_trees) {
    trees <- do.call(rbind, trees)
    if (!is.null(trees)) trees$sample <- NULL
  } else {
    trees <- NULL
  }

  list(
    w0 = storage$samples_w0,
    w1 = storage$samples_w1,
    m_y0 = storage$mean_y0,
    m_y1 = storage$mean_y1,
    check_w0 = storage$check_w0,
    check_w1 = storage$check_w1,
    trees = trees
  )
}


# -----------------------------------------------------------------------------
# Validation and setup helpers
# -----------------------------------------------------------------------------

validate_ordinal_inputs <- function(
  data, n_warmup, n_samples, n_trees, n_thin, n_threads, save_trees_interval
) {
  if (is.null(data)) {
    stop("data must not be NULL")
  }
  if (!all(c("Y", "Z", "W") %in% names(data))) {
    stop("data must contain columns Y, Z, and W")
  }
  numeric_scalars <- c(n_warmup, n_samples, n_trees, n_thin, n_threads
   , save_trees_interval
  )
  if (any(!is.finite(numeric_scalars))) {
    stop("iteration/tree/thread arguments must be finite")
  }
  if (n_warmup < 0 || n_samples <= 0 || n_trees <= 0 ||
      n_thin <= 0 || n_threads <= 0 || save_trees_interval <= 0
  ) {
    stop("n_warmup >= 0 and 
      n_samples/n_trees/n_thin/n_threads/save_trees_interval must be > 0"
    )
  }
}

prepare_ordinal_data <- function(data) {
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  covariate_cols <- setdiff(names(data), c("Y", "Z", "W"))
  X <- scale(as.matrix(data[, covariate_cols, with = FALSE]))

  list(
    X = X,
    Y = as.numeric(data$Y),
    Z = as.numeric(data$Z),
    W = as.numeric(data$W)
  )
}

initialize_ordinal_storage <- function(n_samples, n) {
  list(
    samples_w0 = matrix(NA_real_, n_samples, n),
    samples_w1 = matrix(NA_real_, n_samples, n),
    mean_y0 = matrix(NA_real_, n_samples, n),
    mean_y1 = matrix(NA_real_, n_samples, n),
    check_w0 = matrix(NA_real_, n_samples, n),
    check_w1 = matrix(NA_real_, n_samples, n)
  )
}

initialize_w_states <- function(W, Z, n, init_w_poisson_lambda) {
  w1 <- w0 <- W
  d <- ifelse(runif(n) < 0.03, 0, rpois(n, init_w_poisson_lambda))
  w1[Z == 0] <- pmax(0, (W - d)[Z == 0])
  w0[Z == 1] <- pmax(0, (W + d)[Z == 1])

  list(w0 = w0, w1 = w1, w_max = max(w0, w1))
}

build_w_grid <- function(w_max, monotonicity) {
  if (monotonicity) {
    data.table::rbindlist(lapply(0:w_max, function(x) data.table::CJ(w0 = x, w1 = 0:x)))
  } else {
    expand.grid(w0 = 0:w_max, w1 = 0:w_max) |> data.table::as.data.table()
  }
}

prepare_test_matrices <- function(X, Z, w0, w1, w_grid) {
  z1_idx <- which(Z == 1)
  z0_idx <- which(Z == 0)

  covariate_names <- colnames(X)
  X0 <- as.data.frame(X[z0_idx, , drop = FALSE], check.names = FALSE)
  X1 <- as.data.frame(X[z1_idx, , drop = FALSE], check.names = FALSE)

  xw0 <- cbind(data.frame(id = z0_idx, w0 = w0[z0_idx]), X0)
  xw1 <- cbind(data.frame(id = z1_idx, w1 = w1[z1_idx]), X1)

  grid_df <- as.data.frame(w_grid)
  xw0 <- merge(xw0, grid_df, by = "w0", all.x = TRUE, sort = FALSE)
  xw1 <- merge(xw1, grid_df, by = "w1", all.x = TRUE, sort = FALSE)

  cols <- c("id", covariate_names, "w0", "w1")
  xw <- rbind(xw1[, cols, drop = FALSE], xw0[, cols, drop = FALSE])
  xw <- xw[order(xw$id), , drop = FALSE]

  id <- xw$id
  xw$id <- NULL

  list(X_w_grid = as.matrix(xw), id = id)
}


# -----------------------------------------------------------------------------
# Sampler initialization and update helpers
# -----------------------------------------------------------------------------

initialize_latent_states <- function(w0, w1, lambda_z0, lambda_z1) {
  z1 <- bc(w1 + abs(rnorm(length(w1))), lambda_z1)
  z0 <- bc(w0 + abs(rnorm(length(w0))), lambda_z0)

  list(
    z0 = z0,
    z1 = z1,
    sig_z0 = sd(bc(w0 + 1, lambda_z0), na.rm = TRUE),
    sig_z1 = sd(bc(w1 + 1, lambda_z1), na.rm = TRUE)
  )
}

initialize_ordinal_samplers <- function(
  X, Y, Z, z0, z1, sig_z0, sig_z1, x_plus, X_w_grid, control, k
) {
  normal <- utils::getFromNamespace("normal", "dbarts")
  sampler_z1 <- dbarts(X, z1, test = X, control = control
    , node.prior = normal(k), sigma = sig_z1
  )
  sample_z1 <- sampler_z1$run()

  sampler_z0 <- dbarts(X, z0, test = X, control = control
    , node.prior = normal(k), sigma = sig_z0
  )
  sample_z0 <- sampler_z0$run()

  sampler_y1 <- dbarts(x_plus, Y, subset = Z == 1
    , test = X_w_grid, control = control, node.prior = normal(k)
  )
  sampler_y0 <- dbarts(x_plus, Y, subset = Z == 0
    , test = X_w_grid, control = control, node.prior = normal(k)
  )

  list(
    z0 = sampler_z0,
    z1 = sampler_z1,
    y0 = sampler_y0,
    y1 = sampler_y1,
    state = list(
      mu_z0 = sample_z0$test,
      mu_z1 = sample_z1$test,
      sig_z0 = sample_z0$sigma,
      sig_z1 = sample_z1$sigma
    )
  )
}

set_initial_sampler_data <- function(
  samplers, X, Y, Z, z0, z1, x_plus, X_w_grid
) {
  dt_z1 <- dbartsData(X, z1, test = X)
  dt_z0 <- dbartsData(X, z0, test = X)
  dt_y1 <- dbartsData(x_plus, Y, subset = Z == 1, test = X_w_grid)
  dt_y0 <- dbartsData(x_plus, Y, subset = Z == 0, test = X_w_grid)

  invisible(samplers$z1$setData(dt_z1))
  invisible(samplers$z0$setData(dt_z0))
  invisible(samplers$y1$setData(dt_y1))
  invisible(samplers$y0$setData(dt_y0))
}

sample_latent_z_states <- function(
  w0, w1, w_max, mu_z0, mu_z1, sig_z0, sig_z1,
  lambda_z0, lambda_z1, rho, prev_z0, prev_z1
) {
  z1_lw <- bc(a_j(w1, w_max), lambda_z1)
  z1_up <- bc(a_j(w1 + 1, w_max), lambda_z1)
  z0_lw <- bc(a_j(w0, w_max), lambda_z0)
  z0_up <- bc(a_j(w0 + 1, w_max), lambda_z0)

  if (rho == 0) {
    z1 <- length(w1) |>
      truncnorm::rtruncnorm(a = z1_lw, b = z1_up, mean = mu_z1, sd = sig_z1)
    z0 <- length(w0) |>
      truncnorm::rtruncnorm(a = z0_lw, b = z0_up, mean = mu_z0, sd = sig_z0)
  } else {
    z1z0 <- sample_bivariate_gibbs(
      mu1 = mu_z1,
      sd1 = sig_z1,
      lw1 = z1_lw,
      up1 = z1_up,
      mu2 = mu_z0,
      sd2 = sig_z0,
      lw2 = z0_lw,
      up2 = z0_up,
      rho = rho,
      x1_prev = prev_z1,
      x2_prev = prev_z0
    )
    z1 <- z1z0[, 1]
    z0 <- z1z0[, 2]
  }

  list(z0 = z0, z1 = z1)
}

build_strata_posterior_dt <- function(
  X_w_grid, id, Y, Z, mean_y0, mean_y1, mu_z0, mu_z1, sig_z0, sig_z1
) {
  my0 <- my1 <- my <- pg <- NULL
  numerator_y1 <- numerator_y0 <- sum_num_y1 <- sum_num_y0 <- NULL
  ppw <- safe_ppw <- row_sum <- w0 <- w1 <- NULL
  strata_probs_dt <- data.table::data.table(X_w_grid)
  strata_probs_dt[, id := id]
  strata_probs_dt[, Z := Z[id]]
  strata_probs_dt[, Y := Y[id]]

  strata_probs_dt[, my0 := mean_y0[id]]
  strata_probs_dt[, my1 := mean_y1[id]]
  strata_probs_dt[, my := ifelse(Z == 1, my1, my0)]

  strata_probs_dt[
    , pg := p_wpair(w0, w1, mu_z1[id], mu_z0[id], sig_z1, sig_z0, rho = 0)
  ]
  strata_probs_dt[, pg := pg / sum(pg), by = id]

  strata_probs_dt[, numerator_y1 := pg * my]
  strata_probs_dt[, numerator_y0 := pg * (1 - my)]
  strata_probs_dt[, sum_num_y1 := sum(numerator_y1), by = id]
  strata_probs_dt[, sum_num_y0 := sum(numerator_y0), by = id]

  strata_probs_dt[
    , ppw := ifelse(Y == 1
      , numerator_y1 / sum_num_y1
      , numerator_y0 / sum_num_y0
    )
  ]
  strata_probs_dt[, safe_ppw := ppw]
  strata_probs_dt[is.na(safe_ppw), safe_ppw := 0]
  strata_probs_dt[, row_sum := sum(safe_ppw), by = id]
  strata_probs_dt[row_sum == 0, safe_ppw := 1 / .N, by = id]

  strata_probs_dt
}

sample_strata_membership <- function(strata_probs_dt) {
  id <- safe_ppw <- NULL
  sampled <- strata_probs_dt[, .SD[sample.int(.N, 1, prob = safe_ppw)], by = id]
  list(
    w0 = sampled$w0,
    w1 = sampled$w1,
    my0 = sampled$my0,
    my1 = sampled$my1
  )
}

update_ordinal_outcome_samplers <- function(samplers, x_plus, Y, Z, X_w_grid) {
  dt_y1 <- dbartsData(x_plus, Y, subset = Z == 1, test = X_w_grid)
  dt_y0 <- dbartsData(x_plus, Y, subset = Z == 0, test = X_w_grid)

  invisible(samplers$y1$setData(dt_y1))
  invisible(samplers$y0$setData(dt_y0))
}

extract_ordinal_trees <- function(sampler_z0, sampler_z1, sampler_y0, sampler_y1, iteration) {
  list_trees <- list(
    cbind(m = "z0", sampler_z0$getTrees()),
    cbind(m = "z1", sampler_z1$getTrees()),
    cbind(m = "y0", sampler_y0$getTrees()),
    cbind(m = "y1", sampler_y1$getTrees())
  )
  result <- do.call(rbind, list_trees)
  result$iteration <- iteration
  result
}


# -----------------------------------------------------------------------------
# Statistical utility functions
# -----------------------------------------------------------------------------

sample_bivariate_gibbs <- function(
  mu1, sd1, lw1, up1,
  mu2, sd2, lw2, up2,
  rho = 0.5,
  n_iter = 1,
  x1_prev = 0,
  x2_prev = 0
) {
  n <- length(mu1)
  x1 <- pmax(lw1, pmin(up1, x1_prev))
  x2 <- pmax(lw2, pmin(up2, x2_prev))

  for (iter in seq_len(n_iter)) {
    cond_mu1 <- mu1 + rho * (sd1 / sd2) * (x2 - mu2)
    cond_sd1 <- sqrt(1 - rho^2) * sd1
    x1 <- truncnorm::rtruncnorm(n, a = lw1, b = up1, mean = cond_mu1, sd = cond_sd1)

    cond_mu2 <- mu2 + rho * (sd2 / sd1) * (x1 - mu1)
    cond_sd2 <- sqrt(1 - rho^2) * sd2
    x2 <- truncnorm::rtruncnorm(n, a = lw2, b = up2, mean = cond_mu2, sd = cond_sd2)
  }

  cbind(x1, x2)
}

round_floor <- function(z, y_max = Inf) {
  pmin(floor(z) * as.numeric(z > 0), y_max)
}

bc <- function(t, lambda) {
  if (lambda == 1) {
    t
  } else if (lambda == 0) {
    sign(t) * log(abs(t))
  } else {
    (sign(t) * abs(t)^lambda - 1) / lambda
  }
}

ibc <- function(s, lambda) {
  if (lambda == 1) {
    s
  } else if (lambda == 0) {
    exp(s)
  } else {
    sign(lambda * s + 1) * abs(lambda * s + 1)^(1 / lambda)
  }
}

a_j <- function(j, y_max = Inf) {
  val <- j
  val[j <= 0] <- -Inf
  val[j >= y_max + 1] <- Inf
  val
}

# nothing is calling or using this function. removed.
# mdiff <- function(x) c(x[1], diff(x))

p_wpair <- function(w0, w1, mu_z1, mu_z0, sig_z1, sig_z0, rho = 0) {
  n <- max(length(w0), length(w1))
  result <- numeric(n)

  both_zero <- w0 == 0 & w1 == 0
  w0_zero <- w0 == 0 & !both_zero
  w1_zero <- w1 == 0 & !both_zero
  standard_case <- !both_zero & !w0_zero & !w1_zero

  if (any(both_zero)) {
    mu_z1_vals <- mu_z1[both_zero]
    mu_z0_vals <- mu_z0[both_zero]
    result[both_zero] <- VGAM:::pbinorm(
      0,
      0,
      mu_z1_vals,
      mu_z0_vals,
      sig_z1^2,
      sig_z0^2,
      rho * sig_z1 * sig_z0
    )
  }

  if (any(w0_zero)) {
    w1_vals <- w1[w0_zero]
    mu_z1_vals <- mu_z1[w0_zero]
    mu_z0_vals <- mu_z0[w0_zero]
    result[w0_zero] <- VGAM:::pbinorm(
      0,
      w1_vals + 1,
      mu_z1_vals,
      mu_z0_vals,
      sig_z1^2,
      sig_z0^2,
      rho * sig_z1 * sig_z0
    ) - VGAM:::pbinorm(
      0,
      w1_vals,
      mu_z1_vals,
      mu_z0_vals,
      sig_z1^2,
      sig_z0^2,
      rho * sig_z1 * sig_z0
    )
  }

  if (any(w1_zero)) {
    w0_vals <- w0[w1_zero]
    mu_z1_vals <- mu_z1[w1_zero]
    mu_z0_vals <- mu_z0[w1_zero]
    result[w1_zero] <- VGAM:::pbinorm(
      w0_vals + 1,
      0,
      mu_z1_vals,
      mu_z0_vals,
      sig_z1^2,
      sig_z0^2,
      rho * sig_z1 * sig_z0
    ) - VGAM:::pbinorm(
      w0_vals,
      0,
      mu_z1_vals,
      mu_z0_vals,
      sig_z1^2,
      sig_z0^2,
      rho * sig_z1 * sig_z0
    )
  }

  if (any(standard_case)) {
    w0_std <- w0[standard_case]
    w1_std <- w1[standard_case]
    mu_z1_vals <- mu_z1[standard_case]
    mu_z0_vals <- mu_z0[standard_case]

    upper_prob <- VGAM:::pbinorm(
      w0_std + 1,
      w1_std + 1,
      mu_z1_vals,
      mu_z0_vals,
      sig_z1^2,
      sig_z0^2,
      rho * sig_z1 * sig_z0
    )
    p_w0_w1 <- VGAM:::pbinorm(
      w0_std,
      w1_std,
      mu_z1_vals,
      mu_z0_vals,
      sig_z1^2,
      sig_z0^2,
      rho * sig_z1 * sig_z0
    )
    p_w0p1_w1 <- VGAM:::pbinorm(
      w0_std + 1,
      w1_std,
      mu_z1_vals,
      mu_z0_vals,
      sig_z1^2,
      sig_z0^2,
      rho * sig_z1 * sig_z0
    )
    p_w0_w1p1 <- VGAM:::pbinorm(
      w0_std,
      w1_std + 1,
      mu_z1_vals,
      mu_z0_vals,
      sig_z1^2,
      sig_z0^2,
      rho * sig_z1 * sig_z0
    )
    result[standard_case] <- pmax(0
      , upper_prob - p_w0p1_w1 - p_w0_w1p1 + p_w0_w1
    )
  }

  result
}
