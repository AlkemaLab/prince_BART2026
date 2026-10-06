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
  prob <- prince_bart_fit$probs
  res <- get_mix_tau(prob)
  lapply(res, summary)
}


#' @describeIn estimands Mixed Average Treatment Effect on the Treated
#'   for Compliers (MATT_C)
#' @keywords internal
matt_c <- function(prince_bart_fit) {
  Z <- prince_bart_fit$data$Z
  treated <- Z == 1
  prob <- prince_bart_fit$probs
  res <- get_mix_tau(prob, treated)
  lapply(res, summary)
}


#' @describeIn estimands Sample Average Treatment Effect for Compliers (SATE_C)
#' @keywords internal
sate_c <- function(prince_bart_fit, induce_residual_corr = FALSE) {
  impo <- imput_potentialoutcomes_c(prince_bart_fit)
  impg <- prince_bart_fit$imp
  res <- get_sample_tau(impg, impo, include_corr = induce_residual_corr)
  lapply(res, summary)
}


#' @describeIn estimands Sample Average Treatment Effect on the Treated
#'   for Compliers (SATT_C)
#' @keywords internal
satt_c <- function(prince_bart_fit, induce_residual_corr = FALSE) {
  Z <- prince_bart_fit$data$Z
  treated <- Z == 1
  impo <- imput_potentialoutcomes_c(prince_bart_fit)
  impg <- prince_bart_fit$imp
  res <- get_sample_tau(impg, impo, treated, include_corr = induce_residual_corr)
  lapply(res, summary)
}


#' Compute Mixed Treatment Effects
#'
#' Internal function to compute mixture-based estimands from posterior
#' probability arrays.
#'
#' @param p_arr 4D probability array from prince_bart fit.
#' @param treated Optional logical vector indicating treated units.
#'
#' @return List with strata probabilities and mean effects as draws arrays.
#'
#' @keywords internal
get_mix_tau <- function(p_arr, treated = NULL) {
  if (is.null(treated)) {
    treated <- rep(TRUE, dim(p_arr)[4])
    outname <- c("Y(0) | compliers", "Y(1) | compliers", "Y(0) | never-takers",
                 "Y(1) | always-takers", "Mixed ATE for compliers")
    sname <- c("compliers", "never-takers", "always-takers")
  } else {
    outname <- c("Y(0) | compliers, Z=1", "Y(1) | compliers, Z=1", 
                 "Y(0) | never-takers, Z=1", "Y(1) | always-takers, Z=1", 
                 "Mixed ATT for compliers")
    sname <- c("compliers, Z=1", "never-takers, Z=1", "always-takers, Z=1")
  }

  # Extract strata probabilities for treated units
  # p_arr dimensions: (iter, chain, var, units)
  # Variable indices: 1=p_a, 2=p_n, 3:6=outcomes
  p_a_arr <- p_arr[, , 1, treated, drop = FALSE]  # (iter, chain, 1, n_treated)
  p_n_arr <- p_arr[, , 2, treated, drop = FALSE]  # (iter, chain, 1, n_treated)

  # Compute p_c = 1 - p_n - p_a
  p_c_arr <- 1 - p_n_arr - p_a_arr
  dimnames(p_c_arr)$variable <- "p_c" #otherwise inherits incorrect name

  # Reshape: aperm to move variable dim to end, then select that dim
  # From (iter, chain, 1, units) -> (iter, chain, units, 1)
  p_a <- aperm(p_a_arr, c(1, 2, 4, 3))[, , , 1, drop = FALSE]
  p_n <- aperm(p_n_arr, c(1, 2, 4, 3))[, , , 1, drop = FALSE]
  p_c <- aperm(p_c_arr, c(1, 2, 4, 3))[, , , 1, drop = FALSE]

  # Now p_a, p_n, p_c are (iter, chain, units, 1). Use abind to stack along dim 4
  p_g <- abind::abind(p_c, p_n, p_a, along = 4)  # Result: (iter, chain, units, 3)
  
  strata_prob <- apply(p_g, c(1:2, 4), mean)  # Keep iter, chain, strata; average units

  m_y_arr <- p_arr[, , 3:6, treated, drop = FALSE]  # (iter, chain, 4_vars, n_treated)
  # Reshape: aperm to move var dim to end
  # From (iter, chain, 4, units) -> (iter, chain, units, 4)
  m_y <- aperm(m_y_arr, c(1, 2, 4, 3))

  str_y <- c("m_y0c", "m_y1c", "m_y0n", "m_y1a")
  str_p <- c("p_c", "p_c", "p_n", "p_a")

  mean_pout <- lapply(1:4, function(g) {
    # m_y[, , , g] is (iter, chain, units)
    # p_g[, , , str[g]] is (iter, chain, units)
    numer <- apply(m_y[, , , str_y[g]] * p_g[, , , str_p[g]], 1:2, mean)
    numer / strata_prob[, , str_p[g]]
  })


  mean_effect <- mean_pout
  mean_effect$tau <- mean_pout[[2]] - mean_pout[[1]]
  mean_effect <- abind::abind(mean_effect, along = 3)

  dimnames(mean_effect) <- list(
    iteration = NULL,
    chain = NULL,
    variable = outname
  )
  mean_effect <- posterior::as_draws_array(mean_effect)

  dimnames(strata_prob) <- list(
    iteration = NULL,
    chain = NULL,
    variable = sname
  )
  strata_prob <- posterior::as_draws_array(strata_prob)

  list(strata_prob, mean_effect)
}


#' Compute Sample-Based Treatment Effects
#'
#' Internal function to compute sample-based estimands from imputed
#' compliance classes and outcomes.
#'
#' @param imp_g Imputed compliance class array.
#' @param imp_o Imputed potential outcomes array.
#' @param treated Optional logical vector indicating treated units.
#' @param include_corr Logical; use correlated imputations.
#'
#' @return A draws array with posterior samples.
#'
#' @keywords internal
get_sample_tau <- function(imp_g, imp_o, treated = NULL, include_corr = TRUE) {
  if (is.null(treated)) {
    treated <- rep(TRUE, dim(imp_g)[4])
    outname <- c("Y(0) | compliers", "Y(1) | compliers", "Y(0) | never-takers", 
                 "Y(1) | always-takers", "Sample ATE for compliers")
    sname <- c("compliers", "never-takers", "always-takers")
  } else {
    outname <- c("Y(0) | compliers, Z=1", "Y(1) | compliers, Z=1", 
                 "Y(0) | never-takers, Z=1", "Y(1) | always-takers, Z=1", 
                 "Sample ATT for compliers")
    sname <- c("compliers, Z=1", "never-takers, Z=1", "always-takers, Z=1")
  }

  # Complier indicator
  nt <- imp_g[, , "nt", ]
  at <- imp_g[, , "at", ]
  co <- 1 - nt - at
  p_g <- abind::abind(co, nt, at, along = 4)
  strata_prob <- apply(p_g, c(1:2, 4), mean)

  # Select outcome imputations
  if (include_corr) {
    Y0 <- imp_o[, , "cy0", ]
    Y1 <- imp_o[, , "cy1", ]
  } else {
    Y0 <- imp_o[, , "y0", ]
    Y1 <- imp_o[, , "y1", ]
  }

  Y0n <- Y0
  Y1a <- Y1
  Y0n[co == 1] <- NA
  Y1a[co == 1] <- NA

  # Mask non-compliers
  Y0[co == 0] <- NA
  Y1[co == 0] <- NA

  mY0c <- apply(Y0[, , treated], 1:2, mean, na.rm = TRUE)
  mY1c <- apply(Y1[, , treated], 1:2, mean, na.rm = TRUE)
  mY0n <- apply(Y0n[, , treated], 1:2, mean, na.rm = TRUE)
  mY1a <- apply(Y1a[, , treated], 1:2, mean, na.rm = TRUE)

  tau <- Y1 - Y0
  tau_c <- apply(tau[, , treated], 1:2, mean, na.rm = TRUE)

  mean_effect <- abind::abind(list(mY0c, mY1c, mY0n, mY1a, tau_c), along = 3)

  dimnames(mean_effect) <- list(
    iteration = NULL,
    chain = NULL,
    variable = outname
  )
  mean_effect <- posterior::as_draws_array(mean_effect)

  dimnames(strata_prob) <- list(
    iteration = NULL,
    chain = NULL,
    variable = sname
  )
  strata_prob <- posterior::as_draws_array(strata_prob)

  list(strata_prob, mean_effect)
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
  if (!is.logical(adaptive_levels) ||
        length(adaptive_levels) != 1L || is.na(adaptive_levels)) {
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

  imp <- prince_bart_fit$imp
  probs <- prince_bart_fit$probs

  w0 <- imp[, , "w0", , drop = FALSE][, , 1, , drop = FALSE]
  w1 <- imp[, , "w1", , drop = FALSE][, , 1, , drop = FALSE]
  y0 <- probs[, , "m_y0", , drop = FALSE][, , 1, , drop = FALSE]
  y1 <- probs[, , "m_y1", , drop = FALSE][, , 1, , drop = FALSE]

  # Affected units under monotonicity where instrument shifts uptake by 1.
  co <- (w0 - w1) == 1

  # Mixed contrast: mean Y(1) - Y(0) among affected units.
  d_p <- y1 - y0
  d_p[!co] <- NA

  # Average over units per iteration-chain pair (keep iter x chain shape)
  md_overall <- apply(d_p, c(1, 2), function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  })
  md_overall <- array(md_overall
    , dim = c(dim(md_overall)[1], dim(md_overall)[2], 1)
  )
  md_overall <- posterior::as_draws_array(md_overall)
  dimnames(md_overall)[[3]] <- "Mixed ATE among affected units"

  overall_summary <- summary(md_overall)

  # Determine affected levels and their masses p_j, with j ordered ascending.
  affected_w0 <- as.numeric(w0[co])
  affected_w0 <- affected_w0[!is.na(affected_w0) & affected_w0 >= 1]

  if (length(affected_w0) == 0L) {
    by_w0_summary <- data.frame(
      variable = character(0),
      mean = numeric(0),
      median = numeric(0),
      sd = numeric(0),
      mad = numeric(0),
      q5 = numeric(0),
      q95 = numeric(0)
    )

    return(list(
      overall = overall_summary,
      by_w0 = by_w0_summary,
      grouping = list(
        rule = if (!is.null(level_threshold)) "manual"
        else if (adaptive_levels) "adaptive" else "all_levels"
        , cumulative_mass = cumulative_mass
        , threshold = if (!is.null(level_threshold)) level_threshold 
        else NA_integer_
        , pooled = FALSE
        , level_masses = data.frame(level = integer(0)
          , mass = numeric(0), cumulative_mass = numeric(0)
        )
      )
    ))
  }

  level_table <- table(affected_w0)
  levels_observed <- as.integer(names(level_table))
  ord <- order(levels_observed)
  levels_observed <- levels_observed[ord]
  level_mass <- as.numeric(level_table)[ord] / sum(as.numeric(level_table))
  cum_mass <- cumsum(level_mass)

  # Choose threshold K following precedence: manual > adaptive > all levels.
  if (!is.null(level_threshold)) {
    k <- level_threshold
    rule <- "manual"
  } else if (adaptive_levels) {
    first_reach_idx <- which(cum_mass >= cumulative_mass)[1]
    if (is.na(first_reach_idx)) {
      first_reach_idx <- length(levels_observed)
    }
    k <- levels_observed[first_reach_idx]
    # Safeguard: keep at least first two observed levels when available.
    if (length(levels_observed) >= 2L) {
      k <- max(k, levels_observed[2])
    }
    rule <- "adaptive"
  } else {
    k <- max(levels_observed)
    rule <- "all_levels"
  }

  levels_to_report <- levels_observed[levels_observed <= k]
  levels_above_k <- levels_observed[levels_observed > k]
  pooled <- length(levels_above_k) > 0L

  md_by_w0_list <- lapply(levels_to_report, function(j) {
    co_j <- ((w0 - w1) == 1) & (w0 == j)
    d_p_j <- y1 - y0
    d_p_j[!co_j] <- NA
    md_j <- apply(d_p_j, c(1, 2), function(x) {
      if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    })
    md_j
  })

  labels <- paste0("Level ", levels_to_report)

  if (pooled) {
    co_pool <- ((w0 - w1) == 1) & (w0 > k)
    d_p_pool <- y1 - y0
    d_p_pool[!co_pool] <- NA
    md_pool <- apply(d_p_pool, c(1, 2), function(x) {
      if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    })
    md_by_w0_list[[length(md_by_w0_list) + 1L]] <- md_pool
    labels <- c(labels, paste0("Level > ", k))
  }

  # Combine into draws array with variable names
  md_by_w0 <- abind::abind(md_by_w0_list, along = 3)
  dimnames(md_by_w0) <- list(
    iteration = NULL,
    chain = NULL,
    variable = labels
  )
  md_by_w0 <- posterior::as_draws_array(md_by_w0)
  by_w0_summary <- summary(md_by_w0)

  level_masses <- data.frame(
    level = levels_observed,
    mass = level_mass,
    cumulative_mass = cum_mass
  )

  list(
    overall = overall_summary,
    by_w0 = by_w0_summary,
    grouping = list(
      rule = rule,
      cumulative_mass = cumulative_mass,
      threshold = k,
      pooled = pooled,
      level_masses = level_masses
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
#' @return A 4D array with dimensions (iteration, chain, variable, unit)
#'   containing imputed potential outcomes:
#'   \item{y0}{Imputed Y(0) for compliers}
#'   \item{y1}{Imputed Y(1) for compliers}
#'   \item{cy0}{Imputed Y(0) with induced correlation}
#'   \item{cy1}{Imputed Y(1) with induced correlation}
#'
#' @keywords internal
imput_potentialoutcomes_c <- function(prince_bart_fit) {
  data  <- prince_bart_fit$data
  probs <- prince_bart_fit$probs

  Y <- data$Y
  Z <- data$Z
  N <- length(Y)

  my0co <- probs[, , "m_y0c", ]
  my1co <- probs[, , "m_y1c", ]

  # Standard imputations
  y0c <- apply(my0co, 1:2, function(p) {
    (Z == 1) * stats::rbinom(N, 1, p) + (Z == 0) * Y
  })
  y1c <- apply(my1co, 1:2, function(p) {
    (Z == 1) * Y + (Z == 0) * stats::rbinom(N, 1, p)
  })

  # Imputations with induced correlation
  ycc <- apply(probs[, , c("m_y0c", "m_y1c"), ], 1:2, function(m) {
    my0co <- m["m_y0c", ]
    my1co <- m["m_y1c", ]

    s0 <- cbind(my0co, 1 - my0co) < cbind(my1co, 1 - my1co)
    m0 <- cbind(my0co, 1 - my0co) / (cbind(my1co, 1 - my1co) + .0001)
    k0 <- rowSums(m0 * s0)

    s1 <- cbind(my1co, 1 - my1co) < cbind(my0co, 1 - my0co)
    m1 <- cbind(my1co, 1 - my1co) / (cbind(my0co, 1 - my0co) + .0001)
    k1 <- rowSums(m1 * s1)

    mu0cc <- my0co + (Y - my1co) * k0
    mu1cc <- my1co + (Y - my0co) * k1

    y0cc <- (Z == 1) * stats::rbinom(N, 1, mu0cc) + (Z == 0) * Y
    y1cc <- (Z == 1) * Y + (Z == 0) * stats::rbinom(N, 1, mu1cc)

    data.frame(y0cc, y1cc)
  }, simplify = FALSE)

  y0cc <- apply(ycc, 1:2, function(x) x[[1]]$y0cc)
  y1cc <- apply(ycc, 1:2, function(x) x[[1]]$y1cc)

  imp <- abind::abind(list(y0c, y1c, y0cc, y1cc), along = 4)
  imp <- aperm(imp, c(2, 3, 4, 1))

  dimnames(imp) <- list(
    iteration = NULL,
    chain = NULL,
    variable = c("y0", "y1", "cy0", "cy1"),
    unit = NULL
  )

  imp
}



