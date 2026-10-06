draws_from_chain_matrices <- function(matrices) {
  n_iter <- nrow(matrices[[1]])
  variables <- colnames(matrices[[1]])
  out <- array(NA_real_, c(n_iter, length(matrices), length(variables)),
    dimnames = list(iteration = NULL, chain = NULL, variable = variables))
  for (i in seq_along(matrices)) {
    out[, i, ] <- matrices[[i]]
  }
  posterior::as_draws_array(out)
}

#' Treatment Effect Estimands for Compliers
#'
#' Internal functions to compute various treatment effect estimands for the complier
#' stratum from a fitted prince_bart model.
#'
#' @param prince_bart_fit A fitted object from \code{prince_BART}.
#' @param induce_residual_corr Logical; for sample estimands, whether to induce
#'   residual correlation between potential outcomes (default: FALSE).
#'
#' @return A \code{posterior} summary object with posterior mean, median,
#'   standard deviation, and quantiles.
#'
#' @name estimands
#' @keywords internal
NULL


#' @describeIn estimands Mixed Average Treatment Effect for Compliers (MATE_C)
#' @keywords internal
mate_c <- function(prince_bart_fit) {
  res <- get_mix_tau(prince_bart_fit$chains)
  lapply(res, summary)
}


#' @describeIn estimands Mixed Average Treatment Effect on the Treated
#'   for Compliers (MATT_C)
#' @keywords internal
matt_c <- function(prince_bart_fit) {
  Z <- prince_bart_fit$data$Z
  treated <- Z == 1
  res <- get_mix_tau(prince_bart_fit$chains, treated)
  lapply(res, summary)
}


#' @describeIn estimands Sample Average Treatment Effect for Compliers (SATE_C)
#' @keywords internal
sate_c <- function(prince_bart_fit, induce_residual_corr = FALSE) {
  impo <- imput_potentialoutcomes_c(prince_bart_fit)
  res <- get_sample_tau(prince_bart_fit$chains, impo,
    include_corr = induce_residual_corr)
  lapply(res, summary)
}


#' @describeIn estimands Sample Average Treatment Effect on the Treated
#'   for Compliers (SATT_C)
#' @keywords internal
satt_c <- function(prince_bart_fit, induce_residual_corr = FALSE) {
  Z <- prince_bart_fit$data$Z
  treated <- Z == 1
  impo <- imput_potentialoutcomes_c(prince_bart_fit)
  res <- get_sample_tau(prince_bart_fit$chains, impo, treated,
    include_corr = induce_residual_corr)
  lapply(res, summary)
}


#' Compute Mixed Treatment Effects
#'
#' Internal function to compute mixture-based estimands from posterior
#' per-chain posterior matrices.
#'
#' @param chains List of per-chain draw matrices from a prince_bart fit.
#' @param treated Optional logical vector indicating treated units.
#'
#' @return List with strata probabilities and mean effects as draws arrays.
#'
#' @keywords internal
get_mix_tau <- function(chains, treated = NULL) {
  n <- ncol(chains[[1]]$p_a)
  if (is.null(treated)) {
    treated <- rep(TRUE, n)
    outname <- c("Y(0) | compliers", "Y(1) | compliers",
      "Y(0) | never-takers", "Y(1) | always-takers",
      "Mixed ATE for compliers")
    sname <- c("compliers", "never-takers", "always-takers")
  } else {
    outname <- c("Y(0) | compliers, Z=1", "Y(1) | compliers, Z=1",
      "Y(0) | never-takers, Z=1", "Y(1) | always-takers, Z=1",
      "Mixed ATT for compliers")
    sname <- c("compliers, Z=1", "never-takers, Z=1",
      "always-takers, Z=1")
  }

  quantities <- lapply(chains, function(chain) {
    p_a <- chain$p_a[, treated, drop = FALSE]
    p_n <- chain$p_n[, treated, drop = FALSE]
    p_c <- 1 - p_a - p_n
    mean_p <- cbind(rowMeans(p_c), rowMeans(p_n), rowMeans(p_a))
    colnames(mean_p) <- sname

    weighted_mean <- function(outcome, weight) {
      numerator <- rowMeans(outcome[, treated, drop = FALSE] * weight)
      numerator / rowMeans(weight)
    }
    y0c <- weighted_mean(chain$m_y0c, p_c)
    y1c <- weighted_mean(chain$m_y1c, p_c)
    y0n <- weighted_mean(chain$m_y0n, p_n)
    y1a <- weighted_mean(chain$m_y1a, p_a)
    effect <- cbind(y0c, y1c, y0n, y1a, y1c - y0c)
    colnames(effect) <- outname
    list(strata = mean_p, effect = effect)
  })

  list(
    draws_from_chain_matrices(lapply(quantities, `[[`, "strata")),
    draws_from_chain_matrices(lapply(quantities, `[[`, "effect"))
  )
}


#' Compute Sample-Based Treatment Effects
#'
#' Internal function to compute sample-based estimands from imputed
#' compliance classes and outcomes.
#'
#' @param chains List of per-chain draw matrices.
#' @param imp_o List of per-chain potential-outcome matrices.
#' @param treated Optional logical vector indicating treated units.
#' @param include_corr Logical; use correlated imputations.
#'
#' @return A draws array with posterior samples.
#'
#' @keywords internal
get_sample_tau <- function(chains, imp_o, treated = NULL, include_corr = TRUE) {
  n <- ncol(chains[[1]]$nt)
  if (is.null(treated)) {
    treated <- rep(TRUE, n)
    outname <- c("Y(0) | compliers", "Y(1) | compliers",
      "Y(0) | never-takers", "Y(1) | always-takers",
      "Sample ATE for compliers")
    sname <- c("compliers", "never-takers", "always-takers")
  } else {
    outname <- c("Y(0) | compliers, Z=1", "Y(1) | compliers, Z=1",
      "Y(0) | never-takers, Z=1", "Y(1) | always-takers, Z=1",
      "Sample ATT for compliers")
    sname <- c("compliers, Z=1", "never-takers, Z=1",
      "always-takers, Z=1")
  }

  quantities <- lapply(seq_along(chains), function(i) {
    chain <- chains[[i]]
    co <- 1 - chain$nt - chain$at
    strata <- cbind(rowMeans(co), rowMeans(chain$nt), rowMeans(chain$at))
    colnames(strata) <- sname

    y0 <- imp_o[[i]][[if (include_corr) "cy0" else "y0"]]
    y1 <- imp_o[[i]][[if (include_corr) "cy1" else "y1"]]
    y0n <- y0
    y1a <- y1
    y0n[co == 1] <- NA_real_
    y1a[co == 1] <- NA_real_
    y0[co == 0] <- NA_real_
    y1[co == 0] <- NA_real_

    mean_rows <- function(x) rowMeans(x[, treated, drop = FALSE], na.rm = TRUE)
    effect <- cbind(
      mean_rows(y0), mean_rows(y1), mean_rows(y0n), mean_rows(y1a),
      mean_rows(y1 - y0))
    colnames(effect) <- outname
    list(strata = strata, effect = effect)
  })

  list(
    draws_from_chain_matrices(lapply(quantities, `[[`, "strata")),
    draws_from_chain_matrices(lapply(quantities, `[[`, "effect"))
  )
}


# =============================================================================
# Ordinal Uptake Estimands
# =============================================================================

#' Ordinal Complier-Like Contrast Estimands
#'
#' Compute posterior summaries for contrasts in ordinal uptake settings.
#' Focuses on the monotone compliance group (W(0) - W(1) = 1) and
#' stratifies by baseline uptake W(0) up to level 5.
#'
#' @param prince_bart_fit A fitted ordinal prince_bart object.
#' @param adaptive_levels Logical; if TRUE (default), choose reported levels
#'   adaptively by cumulative affected-unit mass.
#' @param cumulative_mass Numeric in (0, 1]; target cumulative mass used when
#'   \\code{adaptive_levels = TRUE}. Default is \\code{0.80}.
#' @param level_threshold Optional integer threshold K. If supplied, levels
#'   \\code{1..K} are shown individually and higher levels are pooled as
#'   \\code{"Level > K"}; this overrides adaptive selection.
#'
#' @return A list with the following elements:
#'   \item{overall}{Data frame of posterior summary for
#'     \\eqn{E[Y(w=1)-Y(w=0) \\mid W(0)-W(1)=1]}.}
#'   \item{by_w0}{Data frame of posterior summaries by reported level groups.}
#'   \item{grouping}{List with grouping metadata used for level reporting.}
#'
#' @keywords internal
estimands_ordinal_mixed <- function(
  prince_bart_fit,
  adaptive_levels = TRUE,
  cumulative_mass = 0.80,
  level_threshold = NULL
) {
  if (!is.logical(adaptive_levels) || length(adaptive_levels) != 1L ||
      is.na(adaptive_levels)) {
    stop("adaptive_levels must be a single TRUE/FALSE value")
  }
  if (!is.numeric(cumulative_mass) || length(cumulative_mass) != 1L ||
      is.na(cumulative_mass) || cumulative_mass <= 0 || cumulative_mass > 1) {
    stop("cumulative_mass must be a single numeric value in (0, 1]")
  }
  if (!is.null(level_threshold)) {
    if (!is.numeric(level_threshold) || length(level_threshold) != 1L ||
        is.na(level_threshold) || level_threshold < 1 ||
        as.integer(level_threshold) != level_threshold) {
      stop("level_threshold must be NULL or a single integer >= 1")
    }
    level_threshold <- as.integer(level_threshold)
  }

  chains <- prince_bart_fit$chains
  affected <- lapply(chains, function(chain) (chain$w0 - chain$w1) == 1)
  contrast <- lapply(chains, function(chain) chain$m_y1 - chain$m_y0)
  row_effect <- function(values, included) {
    denom <- rowSums(included)
    result <- rowSums(values * included) / denom
    result[denom == 0] <- NA_real_
    result
  }
  overall <- lapply(seq_along(chains), function(i) {
    result <- matrix(row_effect(contrast[[i]], affected[[i]]), ncol = 1)
    colnames(result) <- "Mixed ATE among affected units"
    result
  })
  overall_summary <- summary(draws_from_chain_matrices(overall))

  affected_w0 <- unlist(lapply(seq_along(chains), function(i) {
    chains[[i]]$w0[affected[[i]]]
  }), use.names = FALSE)
  affected_w0 <- affected_w0[!is.na(affected_w0) & affected_w0 >= 1]
  rule <- if (!is.null(level_threshold)) "manual" else
    if (adaptive_levels) "adaptive" else "all_levels"

  if (length(affected_w0) == 0L) {
    return(list(
      overall = overall_summary,
      by_w0 = data.frame(variable = character(0), mean = numeric(0)),
      grouping = list(rule = rule, cumulative_mass = cumulative_mass,
        threshold = if (is.null(level_threshold)) NA_integer_ else level_threshold,
        pooled = FALSE,
        level_masses = data.frame(level = integer(0), mass = numeric(0),
          cumulative_mass = numeric(0)))
    ))
  }

  level_table <- table(affected_w0)
  levels_observed <- as.integer(names(level_table))
  level_mass <- as.numeric(level_table) / sum(level_table)
  cum_mass <- cumsum(level_mass)
  k <- if (!is.null(level_threshold)) {
    level_threshold
  } else if (adaptive_levels) {
    first_reach <- which(cum_mass >= cumulative_mass)[1]
    max(levels_observed[first_reach],
      levels_observed[min(2L, length(levels_observed))])
  } else {
    max(levels_observed)
  }

  levels_to_report <- levels_observed[levels_observed <= k]
  pooled <- any(levels_observed > k)
  labels <- paste0("Level ", levels_to_report)
  if (pooled) labels <- c(labels, paste0("Level > ", k))

  by_level <- lapply(seq_along(chains), function(i) {
    chain <- chains[[i]]
    values <- vapply(levels_to_report, function(j) {
      row_effect(contrast[[i]], affected[[i]] & chain$w0 == j)
    }, numeric(nrow(chain$w0)))
    values <- matrix(values, nrow = nrow(chain$w0))
    if (pooled) {
      values <- cbind(values,
        row_effect(contrast[[i]], affected[[i]] & chain$w0 > k))
    }
    colnames(values) <- labels
    values
  })

  list(
    overall = overall_summary,
    by_w0 = summary(draws_from_chain_matrices(by_level)),
    grouping = list(
      rule = rule, cumulative_mass = cumulative_mass, threshold = k,
      pooled = pooled,
      level_masses = data.frame(level = levels_observed, mass = level_mass,
        cumulative_mass = cum_mass)
    )
  )
}


#' Impute Potential Outcomes for Compliers
#'
#' Imputes potential outcomes Y(0) and Y(1) for compliers based on
#' posterior samples from a fitted prince_bart model. Used internally
#' by sample-based estimand functions.
#'
#' @param prince_bart_fit A fitted object from \code{prince_BART}.
#'
#' @return A list of chains, each containing iteration-by-unit matrices:
#'   \item{y0}{Imputed Y(0) for compliers}
#'   \item{y1}{Imputed Y(1) for compliers}
#'   \item{cy0}{Imputed Y(0) with induced correlation}
#'   \item{cy1}{Imputed Y(1) with induced correlation}
#'
#' @keywords internal
imput_potentialoutcomes_c <- function(prince_bart_fit) {
  Y <- prince_bart_fit$data$Y
  Z <- prince_bart_fit$data$Z
  lapply(prince_bart_fit$chains, function(chain) {
    my0 <- chain$m_y0c
    my1 <- chain$m_y1c
    n_iter <- nrow(my0)
    y_observed <- matrix(rep(Y, each = n_iter), nrow = n_iter)

    draw_binary <- function(prob) {
      matrix(stats::rbinom(length(prob), 1, as.vector(prob)),
        nrow = n_iter, ncol = length(Y))
    }
    y0 <- draw_binary(my0)
    y1 <- draw_binary(my1)
    y0[, Z == 0] <- y_observed[, Z == 0, drop = FALSE]
    y1[, Z == 1] <- y_observed[, Z == 1, drop = FALSE]

    k0 <- (my0 / (my1 + .0001)) * (my0 < my1) +
      ((1 - my0) / (1 - my1 + .0001)) * (my0 > my1)
    k1 <- (my1 / (my0 + .0001)) * (my1 < my0) +
      ((1 - my1) / (1 - my0 + .0001)) * (my1 > my0)
    mu0cc <- my0 + (y_observed - my1) * k0
    mu1cc <- my1 + (y_observed - my0) * k1
    cy0 <- draw_binary(mu0cc)
    cy1 <- draw_binary(mu1cc)
    cy0[, Z == 0] <- y_observed[, Z == 0, drop = FALSE]
    cy1[, Z == 1] <- y_observed[, Z == 1, drop = FALSE]

    list(y0 = y0, y1 = y1, cy0 = cy0, cy1 = cy1)
  })
}
