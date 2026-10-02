test_that("simulation exposes observed data and separate latent quantities", {
  sim <- ps_simulate_data(n = 20, seed = 1)
  expect_named(sim, c("data", "alldata", "weights"))
  expect_named(sim$data, c("x1_1", "x2_1", "Z", "W", "Y"))
  expect_equal(nrow(sim$data), 20L)
  expect_identical(sim$alldata[names(sim$data)], sim$data)
  expect_true(all(sim$data$Y %in% 0:1))
  expect_true(all(sim$data$Z %in% 0:1))
  expect_true(all(sim$alldata$W1 >= sim$alldata$W0))
  expect_equal(sim$data$W, with(sim$alldata, (1 - Z) * W0 + Z * W1))
  expect_identical(as.character(sim$alldata$G),
                   with(sim$alldata, paste0("(", W0, ",", W1, ")")))
})

test_that("coefficients are reusable and seeds preserve caller state", {
  set.seed(804)
  before <- .Random.seed
  weights <- ps_simulate_weights(seed = 10)
  first <- ps_simulate_data(weights = weights, seed = 11)
  expect_identical(.Random.seed, before)
  expect_identical(first, ps_simulate_data(weights = weights, seed = 11))
  second <- ps_simulate_data(weights = weights, seed = 12)
  expect_identical(first$weights, weights)
  expect_identical(second$weights, weights)
  expect_false(identical(first$data, second$data))
  expect_error(ps_simulate_data(weights = "invalid", seed = 13), "weights")
  expect_identical(.Random.seed, before)
  invisible(ps_simulate_data(seed = NULL))
  expect_false(identical(.Random.seed, before))
})

test_that("seeded calls preserve the absence of an initial RNG state", {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)
  if (had_seed) rm(list = ".Random.seed", envir = globalenv())
  invisible(ps_simulate_weights(seed = 20))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  invisible(ps_simulate_data(seed = 20))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("zero covariates and single observations retain their dimensions", {
  for (counts in list(c(0, 0), c(0, 2), c(2, 0))) {
    for (transform in c("linear", "sin", "relu", "sigmoid")) {
      sim <- ps_simulate_data(n = 1, nx1 = counts[1], nx2 = counts[2],
                              x1_transform = transform, seed = 1)
      expect_equal(nrow(sim$data), 1L)
      expect_equal(ncol(sim$data), sum(counts) + 3L)
      expect_equal(dim(sim$weights$g_x1), c(counts[1], 2L, 2L))
      expect_equal(dim(sim$weights$g_x2), c(counts[2], 2L, 2L))
      if (counts[1] == 0 && counts[2] == 2) {
        expect_named(sim$data, c("x2_1", "x2_2", "Z", "W", "Y"))
      }
    }
  }
})

test_that("support zeros exclude strata and custom weights set relative mass", {
  support <- matrix(0, 3, 3)
  support[3, 1] <- 1
  sim <- ps_simulate_data(w_levels = 3, stratum_weights = support, seed = 1)
  expect_true(all(sim$alldata$W0 == 2L))
  expect_true(all(sim$alldata$W1 == 0L))
  expect_identical(levels(sim$alldata$G), "(2,0)")
  ordinal <- ps_simulate_data(w_levels = 4, seed = 2)
  expect_true(all(ordinal$alldata$W1 <= ordinal$alldata$W0))

  # Zero scores leave probabilities proportional to the supplied matrix.
  support <- matrix(c(1, 2, 3, 4), 2, 2)
  sim <- ps_simulate_data(n = 10000, weights = 0,
                          stratum_weights = support, seed = 3)
  expect_equal(as.numeric(table(sim$alldata$G)) / 10000,
               c(0.1, 0.2, 0.3, 0.4), tolerance = 0.02)
  full <- ps_simulate_data(stratum_weights = "full", weights = 0, seed = 4)
  expect_length(levels(full$alldata$G), 4L)
})

test_that("known outcome and instrument models agree with returned truth", {
  weights <- ps_simulate_weights(nx1 = 2, nx2 = 1, weights = 0)
  weights$z_intercept <- -0.2
  weights$z_x1 <- c(0.4, -0.6)
  weights$z_x2 <- 0.3
  # Force compliers and distinguish their coefficients from other strata.
  support <- matrix(0, 2, 2)
  support[1, 2] <- 1
  weights$y_intercept[1, 2] <- -0.5
  weights$y_x1[, 1, 2] <- c(0.7, -0.4)
  weights$y_x2[, 1, 2] <- 0.2
  weights$y_u[1, 2] <- 0.3
  weights$y_w[1, 2] <- 0.8
  for (transform in c("linear", "sin", "relu", "sigmoid")) {
    sim <- ps_simulate_data(nx1 = 2, weights = weights, seed = 5,
                            stratum_weights = support, x1_transform = transform)
    d <- sim$alldata
    f <- switch(transform, linear = identity, sin = sin,
                relu = function(x) pmax(x, 0), sigmoid = plogis)
    expected_e <- plogis(-0.2 + 0.4 * f(d$x1_1) - 0.6 * f(d$x1_2) +
                          0.3 * d$x2_1)
    baseline <- -0.5 + 0.7 * f(d$x1_1) - 0.4 * f(d$x1_2) +
      0.2 * d$x2_1 + 0.3 * d$U
    expect_equal(d$propensity, expected_e)
    expect_equal(d$mean_y0, plogis(baseline))
    expect_equal(d$mean_y1, plogis(baseline + 0.8))
  }
  # With no uptake response, exclusion gives equal conditional outcome means.
  support <- diag(2)
  sim <- ps_simulate_data(nx1 = 2, weights = weights, seed = 6,
                          stratum_weights = support)
  expect_equal(sim$alldata$mean_y0, sim$alldata$mean_y1)
})

test_that("continuous residual scale does not alter the data-generating means", {
  weights <- ps_simulate_weights(weights = 0)
  weights$y_intercept[,] <- 2
  weights$y_w[,] <- 0.5
  one <- ps_simulate_data(weights = weights, y_type = "continuous",
                          sigma_y = 1, seed = 7)
  two <- ps_simulate_data(weights = weights, y_type = "continuous",
                          sigma_y = 2, seed = 7)
  expect_equal(one$alldata$mean_y0, 2 + 0.5 * one$alldata$W0)
  expect_equal(one$alldata$mean_y1, 2 + 0.5 * one$alldata$W1)
  mu <- 2 + 0.5 * one$data$W
  expect_equal(two$data$Y - mu, 2 * (one$data$Y - mu))
  expect_identical(ps_simulate_data(sigma_y = 1, seed = 7),
                   ps_simulate_data(sigma_y = 2, seed = 7))
})

test_that("invalid inputs fail clearly and missing blocks are disclosed", {
  expect_error(ps_simulate_data(n = 0), "n must")
  expect_error(ps_simulate_data(nx1 = -1), "nx1 must")
  expect_error(ps_simulate_data(nx2 = 0.5), "nx2 must")
  expect_error(ps_simulate_data(w_levels = 1), "w_levels must")
  expect_error(ps_simulate_data(seed = "random"), "seed must")
  expect_error(ps_simulate_data(seed = NA_real_), "seed must")
  expect_error(ps_simulate_data(sigma_y = 0), "sigma_y must")
  expect_error(ps_simulate_data(stratum_weights = matrix(0, 2, 2)), "positive")
  expect_error(ps_simulate_data(stratum_weights = matrix(-1, 2, 2)), "nonnegative")
  expect_error(ps_simulate_data(stratum_weights = diag(3)), "2 x 2")
  expect_error(ps_simulate_weights(weights = list(typo = 1)), "unique names")
  expect_error(ps_simulate_weights(weights = list(z_x1 = c(1, 2))), "z_x1")
  expect_error(ps_simulate_weights(weights = list(g_u = matrix(NA, 2, 2))), "g_u")
  expect_error(ps_simulate_weights(weights = list(z_x1 = 1, z_x1 = 2)), "unique")
  expect_warning(w <- ps_simulate_weights(weights = list(z_intercept = 2), seed = 8),
                  "Missing coefficient blocks")
  expect_identical(w$z_intercept, 2)
  expect_length(w, 12L)
  expect_silent(ps_simulate_data(weights = w, seed = 9))
})
