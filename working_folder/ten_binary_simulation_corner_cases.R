# Run from the package root:
# Rscript working_folder/ten_binary_simulation_corner_cases.R
# Or source this file in R after setting the working directory to the package root.
# Each test_that() block contains exactly one expect_*() assertion.

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the princeBART package root.")
}
pkgload::load_all(".", quiet = TRUE)

binary_data <- function() {
  n <- 96L
  data.frame(
    x1 = rep(seq(-1, 1, length.out = 16L), 6L),
    x2 = sin(seq_len(n)),
    Y = rep(c(0, 1), n / 2L),
    Z = rep(c(0, 1), each = n / 2L),
    W = c(rep(0, 32L), rep(1, 16L), rep(0, 16L), rep(1, 32L))
  )
}

fit_binary <- function(data = binary_data(), ...) {
  set.seed(20260922L)
  args <- utils::modifyList(list(
    X = as.matrix(data[c("x1", "x2")]),
    Y = data$Y, Z = data$Z, W = data$W,
    propensity = rep(0.5, nrow(data)),
    uptake_type = "binary", n_warmup = 3L, n_samples = 4L,
    n_chains = 1L, n_trees = 5L, workers = 1L,
    keep_trees = FALSE, verbose = FALSE
  ), list(...), keep.null = TRUE)
  do.call(princeBART::prince_BART, args)
}

run_tests <- function() {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  reporter <- testthat::with_reporter(testthat::SummaryReporter$new(), {
    testthat::test_that("1. observed Z=1, W=0 forces never-taker draws", {
      d <- binary_data()
      fit <- fit_binary(d)
      forced_nt <- d$Z == 1 & d$W == 0
      testthat::expect_true(all(fit$chains[[1]]$nt[, forced_nt, drop = FALSE] == 1))
    })

    testthat::test_that("2. one-chain mixed effect is the weighted contrast", {
      fit <- fit_binary()
      chain <- fit$chains[[1]]
      expected <- vapply(seq_len(4L), function(draw) {
        p_c <- 1 - chain$p_a[draw, ] - chain$p_n[draw, ]
        stats::weighted.mean(chain$m_y1c[draw, ] - chain$m_y0c[draw, ], p_c)
      }, numeric(1L))
      actual <- princeBART:::get_mix_tau(fit$chains)[[2]][,
        , "Mixed ATE for compliers", drop = FALSE]
      testthat::expect_equal(as.numeric(actual), expected, tolerance = 1e-10)
    })

    testthat::test_that("3. missing binary outcome gets an informative error", {
      d <- binary_data()
      d$Y[1] <- NA_real_
      testthat::expect_error(fit_binary(d), "(?i)Y.*(missing|NA|finite)")
    })

    testthat::test_that("4. perfect compliance fits with an empty initial stratum", {
      d <- binary_data()
      d$W <- d$Z
      testthat::expect_s3_class(fit_binary(d), "prince_bart_binary")
    })

    testthat::test_that("5. constant outcomes give finite probability draws", {
      d <- binary_data()
      fits <- lapply(c(0, 1), function(outcome) {
        fit_binary(d, Y = rep(outcome, nrow(d)))
      })
      probability_names <- c("p_a", "p_n", "m_y0c", "m_y1c", "m_y0n", "m_y1a")
      finite <- vapply(fits, function(fit) {
        all(vapply(fit$chains[[1]][probability_names],
          function(x) all(is.finite(x)), logical(1L)))
      }, logical(1L))
      testthat::expect_true(all(finite))
    })

    testthat::test_that("6. boundary Bayes probabilities retain class priors", {
      boundary <- princeBART:::compute_posterior_class_prob(
        Y = c(0, 1), pco = c(0.3, 0.7), pother = c(0.7, 0.3),
        myco = c(0, 1), myother = c(0, 1)
      )
      testthat::expect_equal(as.numeric(boundary), c(0.3, 0.7))
    })

    testthat::test_that("7. observed uptake follows assigned potential uptake", {
      sim <- princeBART::ps_simulate_data(n = 20, seed = 1)
      testthat::expect_equal(sim$data$W,
        with(sim$alldata, (1 - Z) * W0 + Z * W1))
    })

    testthat::test_that("8. one unit and no covariates retain data dimensions", {
      sim <- princeBART::ps_simulate_data(n = 1, nx1 = 0, nx2 = 0, seed = 1)
      testthat::expect_equal(dim(sim$data), c(1L, 3L))
    })

    testthat::test_that("9. custom binary stratum weights set relative mass", {
      support <- matrix(c(1, 2, 3, 4), 2, 2)
      sim <- princeBART::ps_simulate_data(n = 10000, weights = 0,
        stratum_weights = support, seed = 3)
      testthat::expect_equal(as.numeric(table(sim$alldata$G)) / 10000,
        c(0.1, 0.2, 0.3, 0.4), tolerance = 0.02)
    })

    testthat::test_that("10. binary potential-outcome mean follows model coefficients", {
      weights <- princeBART::ps_simulate_weights(weights = 0)
      weights$y_intercept[1, 2] <- -0.5
      weights$y_w[1, 2] <- 0.8
      support <- matrix(0, 2, 2)
      support[1, 2] <- 1
      sim <- princeBART::ps_simulate_data(weights = weights,
        stratum_weights = support, seed = 5)
      testthat::expect_equal(sim$alldata$mean_y1,
        rep(stats::plogis(-0.5 + 0.8), nrow(sim$data)))
    })
  })
  if (reporter$failures$size() > 0L) {
    stop("One or more corner-case expectations failed; see the report above.",
      call. = FALSE)
  }
}

run_tests()
