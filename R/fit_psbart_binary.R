#' Fit Principal Stratification BART (Single Chain)
#'
#' Internal function that fits a single chain of the principal stratification
#' BART model. Users should typically use \code{\link{prince_BART}} instead,
#' which handles data preparation and runs multiple chains.
#'
#' @param X A scaled numeric matrix of covariates with propensity \code{e}
#'   already appended as the last column.
#' @param Y A binary outcome vector (0/1).
#' @param Z A binary treatment assignment/instrument vector (0/1).
#' @param W A binary treatment uptake/received vector (0/1).
#' @param n_warmup Number of warmup/burn-in iterations (default: 1000).
#' @param n_samples Number of posterior samples to collect (default: 1000).
#' @param save_trees Logical; save tree structures for prediction (default: FALSE).
#' @param k Prior hyperparameter for node mean prior (default: 2).
#' @param n_trees Number of trees in the BART ensemble (default: 200).
#' @param verbose Logical; print progress (default: FALSE).
#' @param n_initial Number of initial iterations using MoM offsets (internal, default: 0).
#'
#' @return A list containing:
#'   \item{imputed}{Array of imputed compliance class memberships with dimensions
#'     (n_samples, n, 2) and named variables: "nt" (never-takers), "at" (always-takers)}
#'   \item{probs}{Array of posterior probabilities and outcome means with dimensions
#'     (n_samples, n, 6) and named variables: "p_a" (P(always-taker)), "p_n" (P(never-taker)),
#'     "m_y0c" (\eqn{E[Y(0) \mid complier]}), "m_y1c" (\eqn{E[Y(1) \mid complier]}),
#'     "m_y0n" (\eqn{E[Y(0) \mid never-taker]}), "m_y1a" (\eqn{E[Y(1) \mid always-taker]})}
#'   \item{trees}{Tree structures if \code{save_trees = TRUE}}
#'
#' @keywords internal
.fit_psbart_binary <- function(
  X,
  Y,
  Z,
  W,
  n_warmup = 1000L,
  n_samples = 1000L,
  save_trees = FALSE,
  k = 2,
  n_trees = 200L,
  verbose = FALSE,
  n_initial = 0L
) {

  # X should already be a scaled matrix with e appended
  n <- nrow(X)
  n_total <- n_warmup + n_samples

  # Storage arrays
  s_nt <- s_at <- matrix(NA_real_, n_samples, n)
  trees <- list()
  p_at <- p_nt <- m_y0co <- m_y1co <- m_y0nt <- m_y1at <-
    matrix(NA_real_, n_samples, n)

  # Method-of-moment initial estimates
  if (n_initial > 0) {
    mom <- compute_mom_estimates(X, Y, Z, W)
    intat <- mom$intat; intnt <- mom$intnt; intco <- mom$intco
    inty1at <- mom$inty1at; inty0nt <- mom$inty0nt
    inty1co <- mom$inty1co; inty0co <- mom$inty0co
  } else {
    intat <- intnt <- intco <- inty1at <- inty0nt <- inty1co <- inty0co <- 0.5
  }

  # Initial compliance class assignments
  nt <- (Z == 1) * (W == 0)
  at <- (Z == 0) * (W == 1)
  co <- 1 - nt - at

  control <- dbarts::dbartsControl(
    updateState = TRUE,
    keepTrees = save_trees,
    verbose = FALSE,
    n.burn = 0L,
    n.samples = 1L,
    n.thin = 20L,
    n.chains = 1L,
    n.trees = n_trees,
    n.threads = 1L
  )

  # Initialize samplers
  samplers <- initialize_samplers(X, Y, Z, W, co, at, nt, control, k,
    intco, intat, intnt, inty0nt, inty1at, inty0co, inty1co
  )

  # Sample from prior
  for (s in samplers) s$sampleTreesFromPrior()

  # Main MCMC loop
  for (i in seq_len(n_total)) {
    samples <- run_all_samplers(samplers)

    # Compute probabilities
    pco <- stats::pnorm(samples$co$test)
    patnoco <- stats::pnorm(samples$atnoco$test)
    pat <- patnoco * (1 - pco)
    pnt <- (1 - patnoco) * (1 - pco)

    my0nt <- stats::pnorm(samples$y0nt$test)
    my1at <- stats::pnorm(samples$y1at$test)
    my0co <- stats::pnorm(samples$y0co$test)
    my1co <- stats::pnorm(samples$y1co$test)

    # Update class probabilities using Bayes rule
    pcoy0 <- compute_posterior_class_prob(Y, pco, pnt, my0co, my0nt)
    pcoy1 <- compute_posterior_class_prob(Y, pco, pat, my1co, my1at)

    # Impute compliance classes
    nty <- stats::rbinom(n, 1, (1 - pcoy0))
    aty <- stats::rbinom(n, 1, (1 - pcoy1))
    nt <- (Z == 1) * (W == 0) + (Z == 0) * (W == 0) * nty
    at <- (Z == 0) * (W == 1) + (Z == 1) * (W == 1) * aty
    co <- 1 - nt - at

    # Update sampler data
    update_samplers(samplers, X, Y, Z, co, nt, at, i, n_initial, verbose)

    # Store post-warmup samples
    if (i > n_warmup) {
      offset <- i - n_warmup
      s_nt[offset, ] <- nt
      s_at[offset, ] <- at
      p_at[offset, ] <- pat
      p_nt[offset, ] <- pnt
      m_y0co[offset, ] <- my0co
      m_y1co[offset, ] <- my1co
      m_y0nt[offset, ] <- my0nt
      m_y1at[offset, ] <- my1at

      if (save_trees) {
        trees[[offset]] <- extract_trees(samplers, offset)
      }
    }

    if (verbose && i %% 100 == 0) message("Iteration ", i, " of ", n_total)
  }

  # Prepare output
  if (save_trees) {
    trees <- do.call(rbind, trees)
    trees$sample <- NULL
  } else {
    trees <- NULL
  }

  imputed <- array(c(s_nt, s_at), dim = c(n_samples, n, 2))
  dimnames(imputed) <- list(
    iteration = NULL,
    unit = NULL,
    variable = c("nt", "at")
  )

  probs <- array(
    c(p_at, p_nt, m_y0co, m_y1co, m_y0nt, m_y1at),
    dim = c(n_samples, n, 6)
  )
  dimnames(probs) <- list(
    iteration = NULL,
    unit = NULL,
    variable = c("p_a", "p_n", "m_y0c", "m_y1c", "m_y0n", "m_y1a")
  )

  result <- list(
    imputed = imputed,
    trees = trees,
    probs = probs
  )

  if (verbose) message("Done")
  result
}


# -----------------------------------------------------------------------------
# Internal helper functions
# -----------------------------------------------------------------------------

clip <- function(x, lower = 0.001) {
  upper <- 1 - lower
  x[x < lower] <- lower
  x[x > upper] <- upper
  x
}

compute_mom_estimates <- function(X, Y, Z, W) {
  df <- data.frame(W, X)
  intat <- clip(coef(stats::lm(W ~ ., df[Z == 0, ]))[1])
  intnt <- clip(coef(stats::lm(I(W == 0) ~ ., df[Z == 1, ]))[1])
  intco <- clip(1 - intat - intnt)

  get_int <- function(a, b) {
    s <- W == a & Z == b
    df_y <- data.frame(Y, X)
    if (sum(s) > 4) {
      clip(coef(stats::lm(Y ~ ., df_y[s, ]))[1])
    } else {
      clip(coef(stats::lm(Y ~ ., df_y[W == a, ]))[1])
    }
  }

  inty1at <- get_int(1, 0)
  inty0nt <- get_int(0, 1)
  inty1 <- get_int(1, 1)
  inty0 <- get_int(0, 0)

  intco_noat <- clip(intco / (intco + intnt))
  intco_nont <- clip(intco / (intco + intat))
  inty1co <- clip((inty1 - inty1at * (1 - intco_nont)) / intco_nont)
  inty0co <- clip((inty0 - inty0nt * (1 - intco_noat)) / intco_noat)

  list(
    intat = intat, intnt = intnt, intco = intco,
    inty1at = inty1at, inty0nt = inty0nt, inty1co = inty1co, inty0co = inty0co
  )
}

initialize_samplers <- function(X, Y, Z, W, co, at, nt, control, k,
  intco, intat, intnt, inty0nt, inty1at, inty0co, inty1co
) {

  sampler_co <- dbarts_binary(
    X, co, test = X,
    offset = stats::qnorm(intco),
    control = control,
    k = k
  )

  sampler_atnoco <- dbarts_binary(
    X, at, test = X,
    subset = co == 0,
    offset = stats::qnorm(intat / (intat + intnt)),
    control = control,
    k = k
  )

  nt1 <- if (sum(nt) > 0) nt else (W == 0)
  sampler_y0nt <- dbarts_binary(
    X, Y, test = X,
    subset = nt1 == 1,
    offset = stats::qnorm(inty0nt),
    control = control,
    k = k
  )

  at1 <- if (sum(at) > 0) at else (W == 1)
  sampler_y1at <- dbarts_binary(
    X, Y, test = X,
    subset = at1 == 1,
    offset = stats::qnorm(inty1at),
    control = control,
    k = k
  )

  sampler_y1co <- dbarts_binary(
    X, Y, test = X,
    subset = (co == 1) & (Z == 1),
    offset = stats::qnorm(inty1co),
    control = control,
    k = k
  )

  sampler_y0co <- dbarts_binary(
    X, Y, test = X,
    subset = (co == 1) & (Z == 0),
    offset = stats::qnorm(inty0co),
    control = control,
    k = k
  )

  list(
    co = sampler_co,
    atnoco = sampler_atnoco,
    y0nt = sampler_y0nt,
    y1at = sampler_y1at,
    y0co = sampler_y0co,
    y1co = sampler_y1co
  )
}

run_all_samplers <- function(samplers) {
  list(
    co = samplers$co$run(),
    atnoco = samplers$atnoco$run(),
    y0nt = samplers$y0nt$run(),
    y1at = samplers$y1at$run(),
    y0co = samplers$y0co$run(),
    y1co = samplers$y1co$run()
  )
}

compute_posterior_class_prob <- function(Y, pco, pother, myco, myother) {
  Y * (pco * myco / (pco * myco + pother * myother)) +
    (1 - Y) * (pco * (1 - myco) / (pco * (1 - myco) + pother * (1 - myother)))
}

update_samplers <- function(
  samplers, X, Y, Z, co, nt, at, i, n_initial, verbose
) {
  use_offset <- i >= n_initial

  if (sum(co) > 0) {
    if (use_offset) {
      dt_co <- dbarts::dbartsData(X, co, test = X, offset = 0)
      dt_atnoco <- dbarts::dbartsData(X, at, test = X
        , subset = co == 0, offset = 0
      )
    } else {
      dt_co <- dbarts::dbartsData(X, co, test = X)
      dt_atnoco <- dbarts::dbartsData(X, at, test = X, subset = co == 0)
    }
    samplers$co$setData(dt_co)
    samplers$atnoco$setData(dt_atnoco)
  } else if (verbose) {
    message("Skip co update at iteration ", i)
  }

  if (sum(co * Z) > 0) {
    dt <- if (use_offset) {
      dbarts::dbartsData(X, Y, test = X, subset = co == 1 & Z == 1, offset = 0)
    } else {
      dbarts::dbartsData(X, Y, test = X, subset = co == 1 & Z == 1)
    }
    samplers$y1co$setData(dt)
  }

  if (sum(co * (Z == 0)) > 0) {
    dt <- if (use_offset) {
      dbarts::dbartsData(X, Y, test = X, subset = co == 1 & Z == 0, offset = 0)
    } else {
      dbarts::dbartsData(X, Y, test = X, subset = co == 1 & Z == 0)
    }
    samplers$y0co$setData(dt)
  }

  if (sum(nt) > 0) {
    dt <- if (use_offset) {
      dbarts::dbartsData(X, Y, test = X, subset = nt == 1, offset = 0)
    } else {
      dbarts::dbartsData(X, Y, test = X, subset = nt == 1)
    }
    samplers$y0nt$setData(dt)
  }

  if (sum(at) > 0) {
    dt <- if (i >= 40) {
      dbarts::dbartsData(X, Y, test = X, subset = at == 1, offset = 0)
    } else {
      dbarts::dbartsData(X, Y, test = X, subset = at == 1)
    }
    samplers$y1at$setData(dt)
  }
}

extract_trees <- function(samplers, offset) {
  list_trees <- list(
    cbind(m = "co", samplers$co$getTrees()),
    cbind(m = "atnoco", samplers$atnoco$getTrees()),
    cbind(m = "y1co", samplers$y1co$getTrees()),
    cbind(m = "y0co", samplers$y0co$getTrees()),
    cbind(m = "y1at", samplers$y1at$getTrees()),
    cbind(m = "y0nt", samplers$y0nt$getTrees())
  )
  result <- do.call(rbind, list_trees)
  result$iteration <- offset
  result
}
