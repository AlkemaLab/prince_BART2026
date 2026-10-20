test_that("zero covariates and single observations retain their dimensions", {
  for (counts in list(c(0, 0), c(0, 2), c(2, 0))) {
    for (transform in list(identity, sin, function(x) max(0, x), plogis)) {
      sim <- ps_sim_data(n = 1, nx1 = counts[1], nx2 = counts[2],
                              x1_transform = transform, seed = 1)
      expect_equal(nrow(sim$data), 1L)
      expect_equal(ncol(sim$data), sum(counts) + 3L)
      expect_equal(dim(sim$setup$coefficients$g_x1), c(counts[1], 2L, 2L))
      expect_equal(dim(sim$setup$coefficients$g_x2), c(counts[2], 2L, 2L))
      if (counts[1] == 0 && counts[2] == 2) {
        expect_named(sim$data, c("x2_1", "x2_2", "Z", "W", "Y"))
      }
      expect_equal(sum(grepl("^x1_t_", names(sim$alldata))), counts[1])
    }
  }
})

test_that("support zeros exclude strata and custom weights set relative mass", {
  support <- matrix(0, 3, 3)
  support[3, 1] <- 1
  sim <- ps_sim_data(w_levels = 3, stratum_weights = support, seed = 1)
  expect_identical(sim$setup$stratum_weights, support)
  expect_true(all(sim$alldata$W0 == 2L))
  expect_true(all(sim$alldata$W1 == 0L))
  expect_identical(levels(sim$alldata$G), "(2,0)")
  ordinal <- ps_sim_data(w_levels = 4, seed = 2)
  expect_true(all(ordinal$alldata$W1 <= ordinal$alldata$W0))

  # Zero scores leave probabilities proportional to the supplied matrix.
  support <- matrix(c(1, 2, 3, 4), 2, 2)
  sim <- ps_sim_data(n = 10000, coefficients = 0,
                          stratum_weights = support, seed = 3)
  expect_equal(as.numeric(table(sim$alldata$G)) / 10000,
               c(0.1, 0.2, 0.3, 0.4), tolerance = 0.02)
  full <- ps_sim_data(stratum_weights = matrix(1, 2, 2), coefficients = 0, seed = 4)
  expect_length(levels(full$alldata$G), 4L)
})


test_that("invalid inputs fail clearly and missing blocks are disclosed", {
  expect_error(ps_sim_data(n = 0), "n must")
  expect_error(ps_sim_data(nx1 = -1), "nx1 must")
  expect_error(ps_sim_data(nx2 = 0.5), "nx2 must")
  expect_error(ps_sim_data(w_levels = 1), "w_levels must")
  expect_error(ps_sim_data(seed = "random"), "seed must")
  expect_error(ps_sim_data(seed = NA_real_), "seed must")
  expect_error(ps_sim_data(sigma_y = 0), "sigma_y must")
  expect_error(ps_sim_data(stratum_weights = matrix(0, 2, 2)), "positive")
  expect_error(ps_sim_data(stratum_weights = matrix(-1, 2, 2)), "nonnegative")
  expect_error(ps_sim_data(stratum_weights = diag(3)), "2 x 2")
  expect_error(ps_sim_setup(coefficients = list(typo = 1)), "unique names")
  expect_error(ps_sim_setup(coefficients = list(z_x1 = c(1, 2))), "z_x1")
  expect_error(ps_sim_setup(coefficients = list(g_u = matrix(NA, 2, 2))), "g_u")
  expect_error(ps_sim_setup(coefficients = list(z_x1 = 1, z_x1 = 2)), "unique")
  expect_warning(w <- ps_sim_setup(coefficients = list(z_intercept = 2), seed = 8),
                  "Missing coefficient blocks")
  expect_identical(w$coefficients$z_intercept, 2)
  expect_length(w$coefficients, 12L)
  w$seed <- 9L
  expect_silent(ps_sim_data(w))
})

test_that("data seeds do not randomize or replace setup coefficients", {
  set.seed(20)
  setup <- ps_sim_setup(seed = 1)
  set.seed(20)
  other <- ps_sim_setup(seed = 2)
  expect_identical(setup$coefficients, other$coefficients)

  before <- .Random.seed
  first <- ps_sim_data(setup)
  second <- ps_sim_data(other)
  expect_identical(.Random.seed, before)
  expect_named(first, c("data", "alldata", "setup"))
  expect_identical(first$setup, setup)
  expect_identical(second$setup$coefficients, first$setup$coefficients)
  expect_false(identical(first$data, second$data))
  expect_equal(first$data$W, with(first$alldata, (1 - Z) * W0 + Z * W1))
})

test_that("covariate distributions accept per-predictor settings and scalar functions", {
  sim <- ps_sim_data(n = 12, nx1 = 2, nx2 = 2,
    x1_mean = c(-2, 3), x1_var = 0, x2_p = c(0, 1),
    x1_transform = function(x) max(0, x), coefficients = 0, seed = 1)
  expect_equal(sim$data$x1_1, rep(-2, 12))
  expect_equal(sim$data$x1_2, rep(3, 12))
  expect_equal(sim$data$x2_1, rep(0, 12))
  expect_equal(sim$data$x2_2, rep(1, 12))
  expect_equal(sim$alldata$x1_t_1, rep(0, 12))
  expect_equal(sim$alldata$x1_t_2, rep(3, 12))

  shifted <- sim$setup
  shifted$x1_mean <- 5
  target <- ps_sim_data(shifted)
  expect_equal(target$data$x1_1, rep(5, 12))
  expect_identical(target$setup$coefficients, sim$setup$coefficients)

  normal <- ps_sim_data(n = 5000, x1_mean = 3, x1_var = 4,
    coefficients = 0, seed = 1)
  expect_equal(mean(normal$data$x1_1), 3, tolerance = 0.05)
  expect_equal(stats::var(normal$data$x1_1), 4, tolerance = 0.1)
})

test_that("edited setups must keep support and coefficient dimensions consistent", {
  setup <- ps_sim_setup(coefficients = 0)
  bad <- setup
  bad$w_levels <- 3
  expect_error(ps_sim_data(bad), "3 x 3")
  bad$stratum_weights <- matrix(1, 3, 3)
  expect_error(ps_sim_data(bad), "coefficients\\$g_intercept")

  bad <- setup
  bad$nx1 <- 2
  expect_error(ps_sim_data(bad), "coefficients\\$g_x1")
  bad <- setup
  bad$coefficients$y_w <- NULL
  expect_error(ps_sim_data(bad), "Missing coefficient blocks")
  bad$coefficients <- "random"
  expect_error(ps_sim_data(bad), "complete named list")
  expect_error(ps_sim_data(setup, n = 5), "either setup")
  expect_error(ps_sim_data(list(n = 5)), "named settings")
})

test_that("invalid covariate settings and transformation results fail clearly", {
  expect_error(ps_sim_setup(x1_var = -1), "x1_var")
  expect_error(ps_sim_setup(nx1 = 2, x1_mean = 1:3), "x1_mean")
  expect_error(ps_sim_setup(x1_mean = NA_real_), "x1_mean")
  expect_error(ps_sim_setup(x2_p = 1.1), "x2_p")
  expect_error(ps_sim_setup(x1_transform = "sin"), "function")
  expect_error(ps_sim_data(x1_transform = function(x) c(x, x)), "one finite number")
  expect_error(ps_sim_data(x1_transform = function(x) NA_real_), "one finite number")
})

test_that("fixed coefficients generate the specified instrument and outcome means", {
  setup <- ps_sim_setup(coefficients = 0, x1_transform = sin, seed = 4)
  setup$coefficients$z_intercept <- -0.2
  setup$coefficients$z_x1 <- 0.4
  setup$coefficients$y_intercept[,] <- -0.5
  setup$coefficients$y_x1[,,] <- 0.7
  setup$coefficients$y_u[,] <- 0.3
  setup$coefficients$y_w[,] <- 0.8
  sim <- ps_sim_data(setup)
  d <- sim$alldata
  baseline <- -0.5 + 0.7 * sin(d$x1_1) + 0.3 * d$U
  expect_equal(d$propensity, plogis(-0.2 + 0.4 * sin(d$x1_1)))
  expect_equal(d$mean_y0, plogis(baseline + 0.8 * d$W0))
  expect_equal(d$mean_y1, plogis(baseline + 0.8 * d$W1))

  setup$y_type <- "continuous"
  setup$sigma_y <- 0.5
  continuous <- ps_sim_data(setup)
  expect_equal(continuous$alldata$mean_y0, baseline + 0.8 * d$W0)
  expect_equal(continuous$alldata$mean_y1, baseline + 0.8 * d$W1)
})
