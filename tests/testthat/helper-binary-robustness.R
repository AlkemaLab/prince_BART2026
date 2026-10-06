# Small deterministic fixture: both outcomes in every (Z, W) cell, with
# treatment uptake more common for Z = 1. No external dataset is required.
binary_robustness_data <- function() {
  n <- 96L
  data.frame(
    x1 = rep(seq(-1, 1, length.out = 16L), 6L),
    x2 = sin(seq_len(n)),
    Y = rep(c(0, 1), n / 2L),
    Z = rep(c(0, 1), each = n / 2L),
    W = c(rep(0, 32L), rep(1, 16L), rep(0, 16L), rep(1, 32L))
  )
}

# These short chains test execution and invariants, not convergence or recovery
# of a true causal effect. Supplying propensity avoids a separate BART fit.
fit_binary_robustness <- function(data = binary_robustness_data(), ...) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  old_seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(list = ".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(20260922L)

  args <- utils::modifyList(list(
    X = as.matrix(data[c("x1", "x2")]),
    Y = data$Y,
    Z = data$Z,
    W = data$W,
    propensity = rep(0.5, nrow(data)),
    uptake_type = "binary",
    n_warmup = 3L,
    n_samples = 4L,
    n_chains = 1L,
    n_trees = 5L,
    workers = 1L,
    keep_trees = FALSE,
    verbose = FALSE
  ), list(...), keep.null = TRUE)

  do.call(prince_BART, args)
}

expect_valid_binary_draws <- function(fit, n, chains = 1L) {
  expect_s3_class(fit, "prince_bart")
  expect_s3_class(fit, "prince_bart_binary")
  expect_identical(fit$uptake_type, "binary")
  expect_null(fit$imp)
  expect_null(fit$probs)
  expect_null(fit$trees)
  expect_length(fit$chains, chains)
  for (chain in fit$chains) {
    expect_named(chain,
      c("nt", "at", "p_a", "p_n", "m_y0c", "m_y1c",
        "m_y0n", "m_y1a", "trees"))
    for (name in setdiff(names(chain), "trees")) {
      expect_equal(dim(chain[[name]]), c(4L, n))
    }
    for (name in c("p_a", "p_n", "m_y0c", "m_y1c", "m_y0n", "m_y1a")) {
      expect_true(all(is.finite(chain[[name]])))
      expect_true(all(chain[[name]] >= 0 & chain[[name]] <= 1))
    }
    expect_true(all(chain$p_a + chain$p_n <= 1 + 1e-12))
    expect_true(all(chain$nt %in% c(0, 1)))
    expect_true(all(chain$at %in% c(0, 1)))
    expect_true(all(chain$nt + chain$at <= 1))
    expect_null(chain$trees)
  }
}
