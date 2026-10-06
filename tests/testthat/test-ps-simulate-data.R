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
      expect_equal(sum(grepl("^x1_t_", names(sim$alldata))), counts[1])
    }
  }
})

test_that("support zeros exclude strata and custom weights set relative mass", {
  support <- matrix(0, 3, 3)
  support[3, 1] <- 1
  sim <- ps_simulate_data(w_levels = 3, stratum_weights = support, seed = 1)
  expect_identical(sim$stratum_weights, support)
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
