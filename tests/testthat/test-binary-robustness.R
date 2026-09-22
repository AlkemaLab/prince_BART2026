# Run from the package root:
# testthat::test_local(".", filter = "binary-robustness")
#
# Regression expectations describe correct behavior. Some should fail against
# the current implementation; they are intentionally not skipped or inverted.

test_that("binary fits preserve draw dimensions and observed stratum constraints", {
  d <- binary_robustness_data()
  fit <- fit_binary_robustness(d, n_chains = 2L)
  expect_valid_binary_draws(fit, nrow(d), chains = 2L)

  # Monotonicity identifies these strata directly from observed (Z, W).
  forced_nt <- d$Z == 1 & d$W == 0
  forced_at <- d$Z == 0 & d$W == 1
  expect_true(all(fit$imp[, , "nt", forced_nt, drop = FALSE] == 1))
  expect_true(all(fit$imp[, , "at", forced_at, drop = FALSE] == 1))
  expect_true(all(fit$imp[, , "nt", d$W == 1, drop = FALSE] == 0))
  expect_true(all(fit$imp[, , "at", d$W == 0, drop = FALSE] == 0))
})

test_that("one-chain mixed effects aggregate units without creating extra chains", {
  fit <- fit_binary_robustness()
  expect_valid_binary_draws(fit, nrow(fit$data$X_raw))

  # This is the internal reduction used by coef(type = "mixed") and summary().
  # Compare each draw to the defining weighted mean, using explicit indices
  # so the reference calculation does not share the array-simplification bug.
  mixed <- get_mix_tau(fit$probs)
  expect_equal(dim(mixed[[1]]), c(4L, 1L, 3L))
  expect_equal(dim(mixed[[2]]), c(4L, 1L, 5L))

  expected <- vapply(seq_len(4L), function(draw) {
    p <- fit$probs[draw, 1L, , ]
    p_c <- 1 - p["p_a", ] - p["p_n", ]
    stats::weighted.mean(p["m_y1c", ] - p["m_y0c", ], p_c)
  }, numeric(1L))
  actual <- mixed[[2]][, , "Mixed ATE for compliers", drop = FALSE]
  expect_equal(as.numeric(actual), expected, tolerance = 1e-10)
})

test_that("binary inputs reject invalid values, lengths, and propensity bounds", {
  d <- binary_robustness_data()
  for (variable in c("Y", "Z", "W")) {
    override <- stats::setNames(
      list(replace(d[[variable]], 1L, 2)), variable
    )
    expect_error(
      do.call(fit_binary_robustness, c(list(data = d), override)),
      paste0(variable, ".*binary"),
      info = paste("Invalid binary variable:", variable)
    )
  }

  expect_error(
    fit_binary_robustness(d, W = d$W[-1L]),
    "same number of observations"
  )
  expect_error(
    fit_binary_robustness(d, propensity = rep(0.5, nrow(d) - 1L)),
    "propensity.*length"
  )
  for (boundary in c(0, 1)) {
    expect_error(
      fit_binary_robustness(
        d, propensity = replace(rep(0.5, nrow(d)), 1L, boundary)
      ),
      "propensity.*strictly between"
    )
  }
})

test_that("missing binary inputs produce informative validation errors", {
  d <- binary_robustness_data()
  inputs <- list(Y = d$Y, Z = d$Z, W = d$W,
                 propensity = rep(0.5, nrow(d)))
  for (variable in names(inputs)) {
    override <- stats::setNames(
      list(replace(inputs[[variable]], 1L, NA_real_)), variable
    )
    # A downstream generic sampler error is not sufficient: name the input
    # and explain that missing/nonfinite values are unsupported.
    expect_error(
      do.call(fit_binary_robustness, c(list(data = d), override)),
      paste0("(?i)", variable, ".*(missing|NA|finite)"),
      info = paste("Missing value in:", variable)
    )
  }
})

test_that("constant binary outcomes still produce valid probability draws", {
  d <- binary_robustness_data()
  for (outcome in c(0, 1)) {
    fit <- fit_binary_robustness(d, Y = rep(outcome, nrow(d)))
    expect_valid_binary_draws(fit, nrow(d))
  }
})

test_that("perfect compliance does not fail on an empty initial stratum", {
  d <- binary_robustness_data()
  d$W <- d$Z
  # Valid observed binary data. The initial noncomplier subset is empty;
  # fitting must handle it rather than passing an empty response to dbarts.
  fit <- fit_binary_robustness(d)
  expect_valid_binary_draws(fit, nrow(d))
  # Do not assert that finite-sample W = Z proves every latent unit complies.
})

test_that("binary class probabilities use only the observed outcome likelihood", {
  # Hand-calculated Bayes probabilities with unequal class priors.
  interior <- compute_posterior_class_prob(
    Y = c(1, 0), pco = c(0.25, 0.25), pother = c(0.75, 0.75),
    myco = c(0.8, 0.8), myother = c(0.2, 0.2)
  )
  expect_equal(as.numeric(interior), c(4 / 7, 1 / 13))

  # At these boundaries the observed event is certain in either stratum.
  # Its posterior equals its prior; an unused 0/0 branch must not yield NaN.
  boundary <- compute_posterior_class_prob(
    Y = c(0, 1), pco = c(0.3, 0.7), pother = c(0.7, 0.3),
    myco = c(0, 1), myother = c(0, 1)
  )
  expect_true(all(is.finite(boundary)))
  expect_equal(as.numeric(boundary), c(0.3, 0.7))
})
