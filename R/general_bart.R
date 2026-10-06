#' Generalize Treatment Effects to External Populations (Transportability)
#'
#' Estimates Population Average Treatment Effects (PATE) in an external
#' population using a fitted \code{prince_bart} model and external data
#' , possibly from complex sample surveys.
#' Assumes that covariates X capture all sources of effect heterogeneity,
#' allowing the conditional complier effect CATE_C(x) to generalize to the
#' target population.
#'
#' @param princebart_fit A fitted \code{prince_bart} object with \code{keep_trees = TRUE}.
#' @param newdata A data.frame containing the external population (e.g., survey data).
#'   Covariates present in the source data but missing from \code{newdata} will be
#'   automatically detected and multiply imputed using auxiliary BART models.
#' @param subpop Logical vector of length \code{nrow(newdata)} indicating which
#'   units belong to the target subpopulation. Default is all units.
#' @param psu Vector of primary sampling unit identifiers for \code{newdata}.
#'   Required for complex survey inference.
#' @param weights Vector of survey weights for \code{newdata}. Default is equal weights.
#' @param fast_propensity Logical; if TRUE (default), compute instrument propensity
#'   e = P(Z|X) once using only covariates common to both source and external data,
#'   before imputation. This is much faster. If FALSE, compute e for each 
#'   MI-completed dataset (slower but more statistically more faithful 
#'  when the propensity depends strongly on imputed variables).
#' @param n_cores Number of cores for parallel computation. Default is 1.
#' @param seed Random seed for reproducibility.
#' @param verbose Logical; print progress messages. Default is FALSE
#'
#' @return An object of class \code{general_pate} containing:
#'   \itemize{
#'     \item \code{pate}: Point estimate (posterior mean) of the PATE
#'     \item \code{ci}: 95\% credible interval for PATE
#'     \item \code{sd}: Posterior standard deviation
#'     \item \code{draws}: Vector of posterior draws of the PATE
#'     \item \code{y0}: Array of predicted Y(0) values (units x iterations x chains)
#'     \item \code{y1}: Array of predicted Y(1) values (units x iterations x chains)
#'     \item \code{subpop}: Subpopulation indicator used
#'     \item \code{psu}: PSU identifiers used
#'     \item \code{weights}: Survey weights used
#'   }
#'
#' @details
#' This function implements a multi-step procedure:
#' \enumerate{
#'   \item Multiple imputation of missing covariates in external data using
#'     auxiliary BART models fit on source data.
#'   \item Feature expansion: compute instrument propensity e = P(Z|X) in
#'     external data using BART fit on source data.
#'   \item Predict potential outcomes Y(0) and Y(1) using saved princeBART trees.
#'   \item Estimate PATE using Bayesian bootstrap for complex survey data.
#' }
#'
#' The resulting PATE is a population-level estimand defined over the specified target population.
#'
#' The key identifying assumption is that \eqn{\mathrm{CATE}_C(x)} is transportable,
#' meaning that conditional on X, treatment effects for compliers in the source
#' study equal conditional effects in the target population.
#'
#' For sensitivity analyses (overlap trimming, confounding bounds), use
#' \code{\link{general_BART_overlap}} and \code{\link{general_BART_transportability}}
#' on the returned object.
#'
#' @examples
#' \dontrun{
#' # Fit princeBART on source study
#' fit <- prince_BART(Y ~ X1 + X2 + X3 | Z | W, data = source_data,
#'                  keep_trees = TRUE, n_samples = 1000)
#'
#' # Generalize to external survey population
#' pate <- general_BART(
#'   princebart_fit = fit,
#'   newdata = survey_data,
#'   subpop = survey_data$eligible == 1,
#'   psu = survey_data$cluster_id,
#'   weights = survey_data$survey_weight,
#'   n_cores = 4
#' )
#'
#' # Run sensitivity analyses on the result
#' overlap <- general_BART_overlap(pate, threshold = 0.05)
#' sens <- general_BART_transportability(pate, gamma = 2)
#'
#' print(pate)
#' }
#'
#' @seealso \code{\link{general_BART_overlap}}, \code{\link{general_BART_transportability}}
#'
#' @export
general_BART <- function(
  princebart_fit,
  newdata,
  subpop = NULL,
  psu = NULL,
  weights = NULL,
  fast_propensity = TRUE,
  n_cores = 1L,
  seed = NULL,
  verbose = FALSE
) {

  # Input validation
  if (!inherits(princebart_fit, "prince_bart")) {
    stop("princebart_fit must be a 'prince_bart' object")
  }
  trees_df <- collect_fit_trees(princebart_fit)
  if (is.null(trees_df)) {
    stop("princebart_fit must have saved trees (use keep_trees = TRUE)")
  }

  newdata <- as.data.frame(newdata)
  n_new <- nrow(newdata)

  # Get chain/sample structure from princebart fit
  trees_df <- as.data.frame(trees_df)
  n_chains <- length(unique(trees_df$chain))
  n_samples <- length(unique(trees_df$iteration))

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Default subpopulation: all units
  if (is.null(subpop)) {
    subpop <- rep(TRUE, n_new)
  }

  # Default weights: equal
  if (is.null(weights)) {
    weights <- rep(1, n_new)
  }

  # Default PSU: each unit is its own PSU
  if (is.null(psu)) {
    psu <- seq_len(n_new)
  }

  # Extract source data from princebart fit
  validate_fit_covariate_contract(princebart_fit, require_raw = TRUE)
  source_X <- get_fit_covariates(princebart_fit, type = "model")
  source_Z <- princebart_fit$data$Z
  fit_uptake_type <- if (inherits(princebart_fit, "prince_bart_ordinal")) {
    "ordinal"
  } else if (inherits(princebart_fit, "prince_bart_binary")) {
    "binary"
  } else {
    NA_character_
  }

  source_group_prob <- NULL
  if (identical(fit_uptake_type, "ordinal")) {
    source_group_prob <- Reduce(`+`, lapply(princebart_fit$chains,
      function(chain) colMeans((chain$w0 - chain$w1) == 1,
        na.rm = TRUE))) / length(princebart_fit$chains)
  }

  # Get scaling parameters from canonical fit metadata.
  x_cols <- setdiff(colnames(source_X), "e")
  scaled_center <- princebart_fit$scaling$center
  scaled_scale  <- princebart_fit$scaling$scale

  lambda_z0 <- 1
  lambda_z1 <- 1
  if (!is.null(princebart_fit$call)) {
    if (!is.null(princebart_fit$call$lambda_z0)) {
      lambda_z0 <- suppressWarnings(as.numeric(princebart_fit$call$lambda_z0))
      if (!is.finite(lambda_z0)) lambda_z0 <- 1
    }
    if (!is.null(princebart_fit$call$lambda_z1)) {
      lambda_z1 <- suppressWarnings(as.numeric(princebart_fit$call$lambda_z1))
      if (!is.finite(lambda_z1)) lambda_z1 <- 1
    }
  }
  w_max <- suppressWarnings(max(princebart_fit$data$W, na.rm = TRUE))
  if (!is.finite(w_max)) w_max <- Inf

  # Expand newdata to model-space columns when possible
  # (e.g., ordered factors -> polynomial contrasts like education.L/Q).
  missing_model_cols <- setdiff(colnames(source_X), c(colnames(newdata), "e"))
  if (length(missing_model_cols) > 0) {
    mm <- NULL

    raw_cols <- colnames(princebart_fit$data$X_raw)
    if (length(raw_cols) > 0 && all(raw_cols %in% colnames(newdata))) {
      x_raw_new <- normalize_raw_covariates(newdata[, raw_cols, drop = FALSE])
      mm_try <- try(stats::model.matrix(~ ., data = x_raw_new), silent = TRUE)
      if (!inherits(mm_try, "try-error")) {
        if ("(Intercept)" %in% colnames(mm_try)) {
          mm_try <- mm_try[, colnames(mm_try) != "(Intercept)", drop = FALSE]
        }
        mm <- mm_try
      }
    }

    if (!is.null(mm) && ncol(mm) > 0) {
      mm_df <- as.data.frame(mm)
      addable_cols <- intersect(missing_model_cols, colnames(mm_df))
      if (length(addable_cols) > 0) {
        for (nm in addable_cols) {
          newdata[[nm]] <- mm_df[[nm]]
        }
      }
    }
  }

  # Auto-detect missing covariates (source vars not in newdata)
  mi_vars <- setdiff(x_cols, colnames(newdata))
  common_vars <- intersect(x_cols, colnames(newdata))

  # Notify user about auto-detected missing variables
  if (length(mi_vars) > 0) {
    message("Variables not present in external data (will be imputed): ",
            paste(mi_vars, collapse = ", "))
  }

  # Pre-compute propensity if using fast mode (before imputation)
  e_precomputed <- NULL
  if (fast_propensity) {
    if (verbose) message("Computing propensity on common covariates...")
    propensity_bart_common <- dbarts::bart2(
      source_X[, common_vars, drop = FALSE],
      source_Z,
      keepTrees = TRUE,
      verbose = FALSE
    )
    e_precomputed <- colMeans(
      stats::pnorm(
        stats::predict(propensity_bart_common, newdata = newdata[, common_vars, drop = FALSE])
      )
    )
    e_precomputed <- stats::qnorm(pmax(pmin(e_precomputed, 0.999), 0.001))
  }

  if (verbose) message("Step 1: Imputing missing covariates...")

  # Step 1: Multiple imputation of missing X
  mi_draws <- NULL
  if (length(mi_vars) > 0) {
    mi_draws <- impute_missing_x(
      source_X = source_X[, x_cols, drop = FALSE],
      newdata = newdata,
      mi_vars = mi_vars,
      n_samples = n_samples,
      n_chains = n_chains,
      n_cores = n_cores,
      verbose = verbose
    )
  }

  if (verbose) message("Step 2: Expanding instrument propensity...")

  # Step 2: Fit propensity BART on source data (all covariates)
  # Only needed if not using fast propensity mode
  propensity_bart <- NULL
  if (!fast_propensity) {
    propensity_bart <- dbarts::bart2(
      source_X[, x_cols, drop = FALSE],
      source_Z,
      keepTrees = TRUE,
      verbose = FALSE
    )
  }

  if (verbose) message("Step 3: Predicting potential outcomes...")

  # Step 3: Predict Y(0) and Y(1) in external data
  preds <- predict_external_outcomes(
    trees = trees_df,
    newdata = newdata,
    mi_draws = mi_draws,
    mi_vars = mi_vars,
    propensity_bart = propensity_bart,
    e_precomputed = e_precomputed,
    scaled_center = scaled_center,
    scaled_scale = scaled_scale,
    x_cols = x_cols,
    lambda_z0 = lambda_z0,
    lambda_z1 = lambda_z1,
    w_max = w_max,
    n_samples = n_samples,
    n_chains = n_chains,
    n_cores = n_cores,
    verbose = verbose
  )

  y0 <- preds$y0
  y1 <- preds$y1

  if (verbose) message("Step 4: Computing PATE...")

  # Step 4: Compute PATE using survey-weighted Bayesian bootstrap
  # y0, y1, tau are 3D: [units, iterations, chains]
  tau <- y1 - y0

  # Reshape tau from 3D [units, iterations, chains] to 2D [units, samples]
  # Each (iteration, chain) pair becomes a separate posterior draw
  n_iter <- dim(tau)[2]
  n_ch <- dim(tau)[3]
  tau_2d <- matrix(tau, nrow = n_new, ncol = n_iter * n_ch)

  # Survey PATE with Bayesian bootstrap
  pate_result <- survey_pate(
    tau = tau_2d,
    subpop = subpop,
    psu = psu,
    weights = weights,
    seed = seed
  )

  # Construct result object with all info needed for sensitivity analyses
  result <- list(
    pate = pate_result$estimate,
    ci = pate_result$ci,
    sd = pate_result$sd,
    draws = pate_result$draws,
    y0 = y0,
    y1 = y1,
    tau = tau,
    subpop = subpop,
    psu = psu,
    weights = weights,
    newdata = newdata,
    source_X = source_X,
    source_Z = source_Z,
    fit_uptake_type = fit_uptake_type,
    source_group_prob = source_group_prob,
    trees = trees_df,
    scaled_center = scaled_center,
    scaled_scale = scaled_scale,
    n_obs = sum(subpop),
    seed = seed,
    call = match.call()
  )
  class(result) <- "general_pate"

  if (verbose) message("Done.")
  result
}


# -----------------------------------------------------------------------------
# Helper: Impute missing covariates using auxiliary BART models
# -----------------------------------------------------------------------------
impute_missing_x <- function(
    source_X,
    newdata,
    mi_vars,
    n_samples,
    n_chains,
    n_cores = 1,
    verbose = FALSE
) {
  # Fit auxiliary BART for each missing variable and impute
  # Match the chain/sample structure of the princebart fit

  # Covariates available in both source and newdata (excluding mi_vars)
  common_vars <- intersect(colnames(source_X), colnames(newdata))
  common_vars <- setdiff(common_vars, mi_vars)

  if (length(common_vars) == 0) {
    stop("No common covariates between source and newdata for imputation")
  }

  n_units <- nrow(newdata)

  mi_list <- lapply(mi_vars, function(v) {
    if (verbose) message("  Imputing: ", v)

    # Fit BART on source data: v ~ common_vars
    x_train <- source_X[, common_vars, drop = FALSE]
    y_train <- source_X[, v]

    # Check if binary or continuous
    is_binary <- all(y_train %in% c(0, 1))

    # Generate n_samples x n_chains draws to match princebart structure
    mi_fit <- dbarts::bart2(
      x_train, y_train,
      newdata[, common_vars, drop = FALSE],
      n.chains = n_chains,
      n.samples = n_samples,
      verbose = FALSE
    )

    # yhat.test structure from dbarts::bart2 with n.chains > 1:
    # - 3D array: chains x samples x units
    # - sigma: chains x samples
    # If n.chains = 1: matrix (samples x units), sigma is vector
    yhat <- mi_fit$yhat.test
    sigma <- mi_fit$sigma
    
    # Normalize to 3D: chains x samples x units
    if (length(dim(yhat)) == 2) {
      # Single chain: samples x units -> 1 x samples x units
      yhat <- array(yhat, dim = c(1, dim(yhat)))
      if (!is.null(sigma)) sigma <- matrix(sigma, nrow = 1)
    }
    
    n_ch <- dim(yhat)[1]
    n_samp <- dim(yhat)[2]
    # n_u <- dim(yhat)[3] # removed unused variable
    
    if (is_binary) {
      # Binary: convert to probabilities and sample
      probs <- stats::pnorm(yhat)
      result <- array(stats::rbinom(length(probs), 1, probs), dim = dim(probs))
    } else {
      # Continuous: add noise using sigma (chains x samples)
      if (is.null(sigma)) sigma <- matrix(1, nrow = n_ch, ncol = n_samp)
      # Expand sigma to chains x samples x units
      sigma_expanded <- array(NA_real_, dim = dim(yhat))
      for (ch in seq_len(n_ch)) {
        for (s in seq_len(n_samp)) {
          sigma_expanded[ch, s, ] <- sigma[ch, s]
        }
      }
      result <- yhat + stats::rnorm(length(yhat), 0, sigma_expanded)
    }
    
    # Result is chains x samples x units
    # Permute to units x samples x chains for consistency
    aperm(result, c(3, 2, 1))
  })

  # Keep the same iteration-by-unit matrix contract as the fitted chains.
  lapply(seq_len(n_chains), function(ch) {
    draws <- lapply(mi_list, function(x) {
      t(matrix(x[, , ch], nrow = n_units, ncol = n_samples))
    })
    names(draws) <- mi_vars
    draws
  })
}


# -----------------------------------------------------------------------------
# Helper: Predict Y(0) and Y(1) in external data using saved trees
# -----------------------------------------------------------------------------
predict_external_outcomes <- function(
  trees,
  newdata,
  mi_draws,
  mi_vars,
  propensity_bart,
  e_precomputed = NULL,
  scaled_center,
  scaled_scale,
  x_cols,
  lambda_z0 = 1,
  lambda_z1 = 1,
  w_max = Inf,
  n_samples,
  n_chains,
  n_cores = 1,
  verbose = FALSE
) {

  trees <- as.data.frame(trees)
  n_new <- nrow(newdata)
  
  # Get unique iterations and chains
  iterations <- sort(unique(trees$iteration))
  chains <- sort(unique(trees$chain))
  n_iter <- length(iterations)
  n_ch <- length(chains)

  # Get variable names for scaling (excluding e)
  scale_vars <- names(scaled_center)
  scale_vars <- setdiff(scale_vars, "e")

  # Function to predict for one (iteration, chain) combination
  pred_fun <- function(idx) {
    # Convert linear index to (iteration_idx, chain_idx)
    iter_idx <- ((idx - 1) %% n_iter) + 1
    chain_idx <- ((idx - 1) %/% n_iter) + 1
    
    iter <- iterations[iter_idx]
    ch <- chains[chain_idx]

    # Get base data
    x <- newdata

    # Inject MI imputations if available
    # mi_draws[[chain]][[variable]] is iteration by unit.
    if (!is.null(mi_draws) && length(mi_vars) > 0) {
      for (j in seq_along(mi_vars)) {
        x[[mi_vars[j]]] <- mi_draws[[chain_idx]][[mi_vars[j]]][iter_idx, ]
      }
    }

    # Use pre-computed propensity if available, otherwise compute per-draw
    if (!is.null(e_precomputed)) {
      e_newdata_qnorm <- e_precomputed
    } else {
      e_newdata <- colMeans(
        stats::pnorm(
          stats::predict(propensity_bart, newdata = x[, x_cols, drop = FALSE])
        )
      )
      e_newdata_qnorm <- stats::qnorm(pmax(pmin(e_newdata, 0.999), 0.001))
    }

    # Scale covariates
    x_scaled <- as.matrix(x[, scale_vars, drop = FALSE])
    x_scaled <- scale(x_scaled, center = scaled_center[scale_vars],
                      scale = scaled_scale[scale_vars])

    # Add propensity
    x_scaled <- cbind(x_scaled, e = e_newdata_qnorm)

    # Extract trees for this (iteration, chain) combination
    trees_s <- trees[trees$iteration == iter & trees$chain == ch, ]
    if (nrow(trees_s) == 0) {
      stop(
        "No trees found for iteration ", iter,
        " and chain ", ch,
        "."
      )
    }

    # Binary fits store y0co/y1co, ordinal fits store y0/y1.
    y0_label <- if ("y0co" %in% unique(trees_s$m)) "y0co" else "y0"
    y1_label <- if ("y1co" %in% unique(trees_s$m)) "y1co" else "y1"

    x_for_y <- x_scaled
    if (identical(y0_label, "y0") && identical(y1_label, "y1")) {
      trees_z0 <- trees_s[trees_s$m == "z0", ]
      trees_z1 <- trees_s[trees_s$m == "z1", ]
      if (nrow(trees_z0) == 0 || nrow(trees_z1) == 0) {
        stop(
          "Ordinal prediction requires 'z0' and 'z1' tree components for iteration ",
          iter, " chain ", ch, "."
        )
      }

      mu_z0 <- predict_one_sample_raw(trees_z0, x_scaled)
      mu_z1 <- predict_one_sample_raw(trees_z1, x_scaled)
      w0_hat <- round_floor(ibc(mu_z0, lambda_z0), y_max = w_max)
      w1_hat <- round_floor(ibc(mu_z1, lambda_z1), y_max = w_max)
      x_for_y <- cbind(x_scaled, w0 = w0_hat, w1 = w1_hat)
    }

    # Predict y0
    trees_y0 <- trees_s[trees_s$m == y0_label, ]
    if (nrow(trees_y0) == 0) {
      stop("Missing tree component '", y0_label, "' for iteration ", iter, " chain ", ch, ".")
    }
    y0 <- predict_one_sample(trees_y0, x_for_y)

    # Predict y1
    trees_y1 <- trees_s[trees_s$m == y1_label, ]
    if (nrow(trees_y1) == 0) {
      stop("Missing tree component '", y1_label, "' for iteration ", iter, " chain ", ch, ".")
    }
    y1 <- predict_one_sample(trees_y1, x_for_y)

    list(y0 = y0, y1 = y1, iter_idx = iter_idx, chain_idx = chain_idx)
  }

  # Total number of (iteration, chain) combinations
  n_total <- n_iter * n_ch

  if (n_cores > 1) {
    # Use future for cross-platform parallelization
    old_plan <- future::plan()
    if (inherits(old_plan, "sequential")) {
      future::plan(future::multisession, workers = n_cores)
      on.exit(future::plan(old_plan), add = TRUE)
    }
    
    # Increase globals size limit for large tree objects
    old_max_size <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 2 * 1024^3)  # 2 GiB
    on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)
    
    pred_list <- future.apply::future_lapply(
      seq_len(n_total),
      pred_fun,
      future.seed = TRUE
    )
  } else {
    pred_list <- lapply(seq_len(n_total), pred_fun)
  }

  # Combine into 3D arrays: [units, iterations, chains]
  y0_array <- array(NA_real_, dim = c(n_new, n_iter, n_ch))
  y1_array <- array(NA_real_, dim = c(n_new, n_iter, n_ch))
  
  for (res in pred_list) {
    y0_array[, res$iter_idx, res$chain_idx] <- res$y0
    y1_array[, res$iter_idx, res$chain_idx] <- res$y1
  }

  list(
    y0 = y0_array,
    y1 = y1_array
  )
}


# -----------------------------------------------------------------------------
# Helper: Compute generalizability overlap s = P(group|X) * P(in_source|X)
# -----------------------------------------------------------------------------
#' Compute Generalizability Overlap
#'
#' Estimates the selection score s = P(group|X, in_source) * P(in_source|X)
#' for assessing generalizability from source study to external population.
#' For binary fits, group is the complier stratum. For ordinal fits, group is
#' the affected set defined by W(0)-W(1)=1; this ordinal diagnostic is
#' experimental. Missing covariates are single-imputed using BART mean predictions.
#'
#' @param princebart_fit A fitted \code{prince_bart} object with saved trees.
#' @param newdata External population data. Covariates missing from source model
#'   will be auto-detected and single-imputed using BART.
#' @param weights Survey weights for external data.
#' @param verbose Print progress.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{pi_c}: P(group|X, in_source) for each external unit
#'     \item \code{pi_t}: P(in_source|X) for each external unit
#'     \item \code{pi_s}: Selection score s = pi_c * pi_t
#'     \item \code{e_s_tilde}: Standardized selection score (logit, then z-score)
#'     \item \code{e_s_tilde_source}: Standardized scores for source sample
#'   }
#'
#' @export
compute_generalizability_overlap <- function(
  princebart_fit,
  newdata,
  weights = NULL,
  verbose = FALSE
) {
  validate_fit_covariate_contract(princebart_fit, require_raw = TRUE)

  source_X <- get_fit_covariates(princebart_fit, type = "model")
  x_cols <- setdiff(colnames(source_X), "e")
  n_source <- nrow(source_X)
  n_new    <- nrow(newdata)

  if (is.null(weights)) {
    weights <- rep(1, n_new)
  }

  scaled_center <- princebart_fit$scaling$center
  scaled_scale  <- princebart_fit$scaling$scale
  source_X_unscaled <- as.data.frame(source_X[, x_cols, drop = FALSE])
  colnames(source_X_unscaled) <- x_cols

  # Auto-detect missing covariates
  missing_cols <- setdiff(x_cols, colnames(newdata))

  # Single-impute if newdata is incomplete
  if (length(missing_cols) > 0) {
    if (verbose) message("  Single-imputing missing covariates: "
      , paste(missing_cols, collapse = ", ")
    )

    common_cols <- intersect(x_cols, colnames(newdata))

    for (v in missing_cols) {
      # Use unscaled source data for imputation
      x_train <- source_X_unscaled[, common_cols, drop = FALSE]
      y_train <- source_X_unscaled[[v]]

      is_binary <- all(y_train %in% c(0, 1))

      mi_fit <- dbarts::bart2(
        x_train, y_train,
        newdata[, common_cols, drop = FALSE],
        verbose = FALSE
      )

      # yhat.test may be 3D (chains x samples x units) or 2D (samples x units)
      # Average over all posterior draws to get single imputation per unit
      yhat <- mi_fit$yhat.test
      if (length(dim(yhat)) == 3) {
        # 3D: chains x samples x units -> average over chains and samples
        if (is_binary) {
          newdata[[v]] <- apply(stats::pnorm(yhat), 3, mean)
        } else {
          newdata[[v]] <- apply(yhat, 3, mean)
        }
      } else {
        # 2D: samples x units -> colMeans
        if (is_binary) {
          newdata[[v]] <- colMeans(stats::pnorm(yhat))
        } else {
          newdata[[v]] <- colMeans(yhat)
        }
      }
    }
  }

  # Now newdata has all x_cols - stack for P(in_source|X) in original scale
  stacked_data <- rbind(
    data.frame(source_X_unscaled[, x_cols, drop = FALSE], source = 1),
    data.frame(newdata[, x_cols, drop = FALSE], source = 0)
  )

  if (verbose) message("  Fitting P(in_source|X)...")

  # Fit P(in_source|X) using all covariates
  pi_t_bart <- dbarts::bart2(
    stacked_data[, x_cols, drop = FALSE],
    stacked_data$source,
    verbose = FALSE
  )

  pi_t_all <- stats::fitted(pi_t_bart)
  pi_t_source <- pi_t_all[seq_len(n_source)]
  pi_t_new <- pi_t_all[seq(n_source + 1, n_source + n_new)]

  if (verbose) message("  Predicting P(complier|X, in_source)...")

  # Predict P(group|X):
  # - binary fit -> complier probability via "co" trees
  # - ordinal fit -> affected probability via imp draws
  trees <- as.data.frame(collect_fit_trees(princebart_fit))
  trees_co <- trees[trees$m == "co", ]

  # Fit instrument propensity e = P(Z|X) using unscaled covariates
  source_Z <- princebart_fit$data$Z
  if (verbose) message("  Fitting P(Z|X)...")

  propensity_bart <- dbarts::bart2(
    source_X_unscaled[, x_cols, drop = FALSE],
    source_Z,
    keepTrees = TRUE,
    verbose = FALSE
  )

  # Predict propensity for external data (now complete)
  # predict() returns matrix: rows = posterior samples, columns = units
  e_new <- colMeans(stats::pnorm(
    stats::predict(propensity_bart, newdata = newdata[, x_cols, drop = FALSE])
  ))
  e_new <- pmax(pmin(e_new, 0.999), 0.001)

  # Scale covariates for tree prediction (trees expect scaled X)
  # Use stored scaling parameters
  x_scaled_new <- as.matrix(newdata[, x_cols, drop = FALSE])
  x_scaled_new <- scale(x_scaled_new,
                        center = scaled_center,
                        scale = scaled_scale)
  x_scaled_new <- cbind(x_scaled_new, e = stats::qnorm(e_new))
  
  # Trees expect scaled model-space covariates + unscaled propensity column e
  source_x_scaled <- as.matrix(source_X[, x_cols, drop = FALSE])
  source_x_scaled <- scale(source_x_scaled,
                           center = scaled_center,
                           scale = scaled_scale)
  source_x_scaled <- cbind(source_x_scaled, e = source_X[, "e"])

  if (inherits(princebart_fit, "prince_bart_binary")) {
    # Binary fit path: use saved complier trees.
    if (nrow(trees_co) == 0) {
      stop("Binary overlap requires saved 'co' trees, but none were found")
    }
    samples <- unique(trees_co$iteration)

    pi_c_mat <- sapply(samples, function(s) {
      trees_s <- trees_co[trees_co$iteration == s, ]
      predict_one_sample(trees_s, x_scaled_new)
    })
    pi_c_mat <- as.matrix(pi_c_mat)
    if (nrow(pi_c_mat) != n_new) {
      pi_c_mat <- matrix(pi_c_mat, nrow = n_new)
    }
    pi_c_new <- rowMeans(pi_c_mat)

    pi_c_source_mat <- sapply(samples, function(s) {
      trees_s <- trees_co[trees_co$iteration == s, ]
      predict_one_sample(trees_s, source_x_scaled)
    })
    pi_c_source_mat <- as.matrix(pi_c_source_mat)
    if (nrow(pi_c_source_mat) != n_source) {
      pi_c_source_mat <- matrix(pi_c_source_mat, nrow = n_source)
    }
    pi_c_source <- rowMeans(pi_c_source_mat)
  } else if (inherits(princebart_fit, "prince_bart_ordinal")) {
    # Ordinal fit path: estimate affected probability from imp draws.
    if (verbose) {
      message("  Using experimental ordinal overlap diagnostic (affected-unit probability)")
    }
    pi_c_source <- Reduce(`+`, lapply(princebart_fit$chains,
      function(chain) colMeans((chain$w0 - chain$w1) == 1,
        na.rm = TRUE))) / length(princebart_fit$chains)

    pi_c_bart <- dbarts::bart2(
      source_X_unscaled[, x_cols, drop = FALSE],
      pi_c_source,
      keepTrees = TRUE,
      verbose = FALSE
    )

    pi_c_new <- colMeans(
      stats::pnorm(
        stats::predict(pi_c_bart, newdata = newdata[, x_cols, drop = FALSE])
      )
    )
    pi_c_new <- pmax(pmin(pi_c_new, 0.999), 0.001)
    pi_c_source <- pmax(pmin(pi_c_source, 0.999), 0.001)
  } else {
    stop("Unsupported prince_bart subclass for overlap computation")
  }

  # Compute selection scores
  pi_s_new <- pi_c_new * pi_t_new
  pi_s_source <- pi_c_source * pi_t_source

  # Standardize using source distribution (weighted)
  source_wt <- rep(1, n_source)  # Source weights if available
  e_s_source <- stats::qlogis(pmax(pmin(pi_s_source, 0.999), 0.001))
  e_s_new <- stats::qlogis(pmax(pmin(pi_s_new, 0.999), 0.001))

  m <- stats::weighted.mean(e_s_source, source_wt, na.rm = TRUE)
  m2 <- stats::weighted.mean(e_s_source^2, source_wt, na.rm = TRUE)
  s <- sqrt(m2 - m^2)

  e_s_tilde_source <- (e_s_source - m) / s
  e_s_tilde_new <- (e_s_new - m) / s

  list(
    pi_c = pi_c_new,
    pi_t = pi_t_new,
    pi_s = pi_s_new,
    e_s_tilde = e_s_tilde_new,
    e_s_tilde_source = e_s_tilde_source
  )
}


# -----------------------------------------------------------------------------
# Helper: Survey-weighted PATE with Bayesian bootstrap
# -----------------------------------------------------------------------------

#' Bayesian Bootstrap for Dirichlet weights
#' @keywords internal
bayesian_bootstrap <- function(n) {
  u <- stats::runif(n - 1)
  diff(c(0, sort(u), 1))
}


survey_pate <- function(tau, subpop, psu, weights, seed = NULL) {
  tau <- as.matrix(tau)
  tau[!subpop, ] <- NA

  n_psu <- length(unique(psu))
  n_draws <- ncol(tau)

  # PSU-level weights (mean within PSU)
  p_c <- tapply(weights, psu, mean)
  n_c <- tapply(!is.na(psu), psu, sum)

  if (!is.null(seed)) set.seed(seed)
  bb_list <- lapply(seq_len(n_draws), function(r) bayesian_bootstrap(n_psu))

  # Compute weighted mean for each posterior draw
  pate_draws <- sapply(seq_len(n_draws), function(k) {
    psu_tau <- tapply(tau[, k], psu, mean, na.rm = TRUE)
    stats::weighted.mean(psu_tau, bb_list[[k]] * p_c * n_c, na.rm = TRUE)
  })

  list(
    estimate = mean(pate_draws, na.rm = TRUE),
    sd = stats::sd(pate_draws, na.rm = TRUE),
    ci = stats::quantile(pate_draws, c(0.025, 0.975), na.rm = TRUE),
    draws = pate_draws
  )
}


# -----------------------------------------------------------------------------
# Sensitivity Function: Overlap Analysis
# -----------------------------------------------------------------------------

#' Generalizability Overlap Analysis for PATE Estimates
#'
#' Compute generalizability overlap scores and optionally trim observations
#' with low overlap to produce a trimmed PATE estimate.
#' For binary fits, overlap is based on complier similarity. For ordinal fits,
#' overlap uses an affected-unit analogue (W(0)-W(1)=1), which is experimental.
#'
#' @param object A `general_pate` object from [general_BART()]
#' @param threshold Numeric threshold for trimming based on standardized
#'   selection scores. Observations with `|e_s_tilde| > threshold` are
#'   excluded. Default is NULL (no trimming, just compute overlap).
#' @param overlap_value Value to replace `tau` with for trimmed observations.
#'   Either `"zero"` (set to 0) or `"NA"` (exclude completely). Default is `"zero"`.
#' @param verbose Logical; print progress messages. Default is FALSE.
#'
#' @return A list with components:
#' \describe{
#'   \item{overlap}{Data frame with overlap metrics: pi_c, pi_t, pi_s, e_s_tilde
#'   where pi_c is group probability (complier for binary; affected for ordinal)}
#'   \item{e_s_tilde_source}{Standardized selection scores for source data}
#'   \item{n_trimmed}{Number of observations trimmed (if threshold used)}
#'   \item{pate_trimmed}{Trimmed PATE estimate (if threshold used)}
#'   \item{ci_trimmed}{95 percent CI for trimmed PATE (if threshold used)}
#' }
#'
#' @seealso [general_BART()], [general_BART_transportability()]
#'
#' @examples
#' \dontrun{
#' fit <- general_BART(princebart_fit, newdata, ...)
#' overlap <- general_BART_overlap(fit, threshold = 2)
#' }
#'
#' @export
general_BART_overlap <- function(
  object
  , threshold = NULL
  , overlap_value = c("zero", "NA")
  , verbose = FALSE
) {
  if (!inherits(object, "general_pate")) {
    stop("object must be a 'general_pate' object from general_BART()")
  }

  overlap_value <- match.arg(overlap_value)

  if (verbose) message("Computing generalizability overlap...")

  # Use stored data from general_pate object
  newdata <- object$newdata
  source_X <- object$source_X
  source_Z <- object$source_Z
  trees <- object$trees
  scaled_center <- object$scaled_center
  scaled_scale <- object$scaled_scale
  fit_uptake_type <- object$fit_uptake_type
  source_group_prob <- object$source_group_prob
  # weights <- object$weights # removed unused variable

  if (is.null(scaled_center) || is.null(scaled_scale)) {
    stop(
      "general_BART_overlap() requires scaling metadata in general_pate object. ",
      "Please rerun general_BART() with the current package version."
    )
  }

  if (is.null(fit_uptake_type)) {
    trees_tmp <- as.data.frame(trees)
    fit_uptake_type <- if (any(trees_tmp$m == "co")) "binary" else "ordinal"
  }

  x_cols <- setdiff(colnames(source_X), "e")
  n_source <- nrow(source_X)
  n_new <- nrow(newdata)
  source_X_unscaled <- as.data.frame(source_X[, x_cols, drop = FALSE])
  colnames(source_X_unscaled) <- x_cols

  # Auto-detect and single-impute missing covariates
  missing_cols <- setdiff(x_cols, colnames(newdata))

  if (length(missing_cols) > 0) {
    if (verbose) message("  Single-imputing missing covariates: ",
                         paste(missing_cols, collapse = ", "))

    common_cols <- intersect(x_cols, colnames(newdata))

    for (v in missing_cols) {
      # Use unscaled source data for imputation model
      x_train <- source_X_unscaled[, common_cols, drop = FALSE]
      y_train <- source_X_unscaled[[v]]

      is_binary <- all(y_train %in% c(0, 1))

      mi_fit <- dbarts::bart2(
        x_train, y_train,
        newdata[, common_cols, drop = FALSE],
        verbose = FALSE
      )

      yhat <- mi_fit$yhat.test
      if (length(dim(yhat)) == 3) {
        if (is_binary) {
          newdata[[v]] <- apply(stats::pnorm(yhat), 3, mean)
        } else {
          newdata[[v]] <- apply(yhat, 3, mean)
        }
      } else {
        if (is_binary) {
          newdata[[v]] <- colMeans(stats::pnorm(yhat))
        } else {
          newdata[[v]] <- colMeans(yhat)
        }
      }
    }
  }

  # Stack source and target for P(in_source|X) - both in original scale
  stacked_data <- rbind(
    data.frame(source_X_unscaled[, x_cols, drop = FALSE], source = 1),
    data.frame(newdata[, x_cols, drop = FALSE], source = 0)
  )

  if (verbose) message("  Fitting P(in_source|X)...")

  pi_t_bart <- dbarts::bart2(
    stacked_data[, x_cols, drop = FALSE],
    stacked_data$source,
    verbose = FALSE
  )

  pi_t_all <- stats::fitted(pi_t_bart)
  pi_t_source <- pi_t_all[seq_len(n_source)]
  pi_t_new <- pi_t_all[seq(n_source + 1, n_source + n_new)]

  if (verbose) message("  Predicting P(group|X, in_source)...")

  # Predict P(group|X)
  trees_df <- as.data.frame(trees)
  trees_co <- trees_df[trees_df$m == "co", ]

  # Fit propensity P(Z|X) using unscaled source data
  if (verbose) message("  Fitting P(Z|X)...")

  propensity_bart <- dbarts::bart2(
    source_X_unscaled[, x_cols, drop = FALSE],
    source_Z,
    keepTrees = TRUE,
    verbose = FALSE
  )

  e_new <- colMeans(stats::pnorm(
    stats::predict(propensity_bart, newdata = newdata[, x_cols, drop = FALSE])
  ))
  e_new <- pmax(pmin(e_new, 0.999), 0.001)

  # Scale covariates for tree prediction
  x_scaled_new <- as.matrix(newdata[, x_cols, drop = FALSE])
  x_scaled_new <- scale(x_scaled_new, center = scaled_center, scale = scaled_scale)
  x_scaled_new <- cbind(x_scaled_new, e = stats::qnorm(e_new))

  source_x_scaled <- as.matrix(source_X[, x_cols, drop = FALSE])
  source_x_scaled <- scale(source_x_scaled, center = scaled_center, scale = scaled_scale)
  source_x_scaled <- cbind(source_x_scaled, e = source_X[, "e"])

  if (identical(fit_uptake_type, "binary")) {
    if (nrow(trees_co) == 0) {
      stop("Binary overlap requires saved 'co' trees, but none were found")
    }

    samples <- unique(trees_co$iteration)
    pi_c_mat <- sapply(samples, function(s) {
      trees_s <- trees_co[trees_co$iteration == s, ]
      predict_one_sample(trees_s, x_scaled_new)
    })
    pi_c_mat <- as.matrix(pi_c_mat)
    if (nrow(pi_c_mat) != n_new) {
      pi_c_mat <- matrix(pi_c_mat, nrow = n_new)
    }
    pi_c_new <- rowMeans(pi_c_mat)

    pi_c_source_mat <- sapply(samples, function(s) {
      trees_s <- trees_co[trees_co$iteration == s, ]
      predict_one_sample(trees_s, source_x_scaled)
    })
    pi_c_source_mat <- as.matrix(pi_c_source_mat)
    if (nrow(pi_c_source_mat) != n_source) {
      pi_c_source_mat <- matrix(pi_c_source_mat, nrow = n_source)
    }
    pi_c_source <- rowMeans(pi_c_source_mat)
  } else if (identical(fit_uptake_type, "ordinal")) {
    if (verbose) {
      message("  Using experimental ordinal overlap diagnostic (affected-unit probability)")
    }
    if (is.null(source_group_prob) || length(source_group_prob) != n_source) {
      stop(
        "Ordinal overlap requires source_group_prob in the general_pate object. ",
        "Please rerun general_BART() with the current package version."
      )
    }

    pi_c_source <- as.numeric(source_group_prob)
    pi_c_source <- pmax(pmin(pi_c_source, 0.999), 0.001)

    pi_c_bart <- dbarts::bart2(
      source_X_unscaled[, x_cols, drop = FALSE],
      pi_c_source,
      keepTrees = TRUE,
      verbose = FALSE
    )

    pi_c_new <- colMeans(
      stats::pnorm(
        stats::predict(pi_c_bart, newdata = newdata[, x_cols, drop = FALSE])
      )
    )
    pi_c_new <- pmax(pmin(pi_c_new, 0.999), 0.001)
  } else {
    stop("Unsupported or unknown fit_uptake_type in general_pate object")
  }

  # Compute selection scores
  pi_s_new <- pi_c_new * pi_t_new
  pi_s_source <- pi_c_source * pi_t_source

  # Standardize using source distribution
  e_s_source <- stats::qlogis(pmax(pmin(pi_s_source, 0.999), 0.001))
  e_s_new <- stats::qlogis(pmax(pmin(pi_s_new, 0.999), 0.001))

  m <- mean(e_s_source, na.rm = TRUE)
  s <- stats::sd(e_s_source, na.rm = TRUE)

  e_s_tilde_source <- (e_s_source - m) / s
  e_s_tilde_new <- (e_s_new - m) / s

  overlap_df <- data.frame(
    pi_c = pi_c_new,
    pi_t = pi_t_new,
    pi_s = pi_s_new,
    e_s_tilde = e_s_tilde_new
  )

  result <- list(
    overlap = overlap_df,
    e_s_tilde_source = e_s_tilde_source,
    n_trimmed = 0,
    pate_trimmed = NULL,
    ci_trimmed = NULL
  )

  # Apply trimming if threshold specified

  if (!is.null(threshold)) {
    if (verbose) message("  Applying overlap trimming (threshold = ", threshold, ")...")

    trim_mask <- abs(e_s_tilde_new) > threshold
    n_trimmed <- sum(trim_mask)

    if (verbose) message("  Trimmed ", n_trimmed, " of ", n_new, " observations")

    # Get tau from object and apply trimming
    tau <- object$tau
    n_iter <- dim(tau)[2]
    n_ch <- dim(tau)[3]
    tau_2d <- matrix(tau, nrow = n_new, ncol = n_iter * n_ch)

    if (overlap_value == "zero") {
      tau_2d[trim_mask, ] <- 0
    } else {
      tau_2d[trim_mask, ] <- NA
    }

    # Recompute PATE with trimmed tau
    pate_trimmed <- survey_pate(
      tau = tau_2d,
      subpop = object$subpop,
      psu = object$psu,
      weights = object$weights,
      seed = object$seed
    )

    result$n_trimmed <- n_trimmed
    result$pate_trimmed <- pate_trimmed$estimate
    result$ci_trimmed <- pate_trimmed$ci
    result$sd_trimmed <- pate_trimmed$sd
    result$draws_trimmed <- pate_trimmed$draws
  }

  if (verbose) message("Done.")
  result
}


# -----------------------------------------------------------------------------
# Sensitivity Function: Weight-Shift Bounds
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# S3 Methods for general_pate class
# -----------------------------------------------------------------------------

#' @export
print.general_pate <- function(x, ...) {
  cat("General BART PATE Estimate\n")
  cat("==========================\n")
  cat("PATE:    ", round(x$pate, 4), "\n")
  cat("95% CI:  [", round(x$ci[1], 4), ", ", round(x$ci[2], 4), "]\n", sep = "")
  cat("SD:      ", round(x$sd, 4), "\n")
  cat("N (subpop):", x$n_obs, "\n")

  cat("\nUse general_BART_overlap() for overlap analysis")
  cat("\nUse general_BART_transportability() for sensitivity analysis\n")

  invisible(x)
}


#' @export
summary.general_pate <- function(object, ...) {
  cat("General BART PATE Summary\n")
  cat("=========================\n\n")

  cat("Treatment Effect Estimate:\n")
  cat("  PATE:          ", round(object$pate, 4), "\n")
  cat("  Posterior SD:  ", round(object$sd, 4), "\n")
  cat("  95% CI:        [", round(object$ci[1], 4), ", ",
      round(object$ci[2], 4), "]\n", sep = "")
  cat("  N (subpop):    ", object$n_obs, "\n\n")

  cat("For sensitivity analyses, use:\n")
  cat("  general_BART_overlap(object)      - Generalizability overlap\n")
  cat("  general_BART_transportability(object)  - Weight-shift bounds\n")

  invisible(object)
}


#' Plot overlap diagnostics
#'
#' @param x Output from [general_BART_overlap()]
#' @param ... Additional arguments (ignored)
#'
#' @export
plot_overlap <- function(x, ...) {
  if (!is.list(x) || is.null(x$overlap) || is.null(x$e_s_tilde_source)) {
    stop("x must be output from general_BART_overlap()")
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  graphics::par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

  # Density plot of selection scores
  e_source <- x$e_s_tilde_source
  e_target <- x$overlap$e_s_tilde

  xlim <- range(c(e_source, e_target), na.rm = TRUE)

  d_source <- stats::density(e_source, na.rm = TRUE)
  d_target <- stats::density(e_target, na.rm = TRUE)

  graphics::plot(d_source, xlim = xlim,
                 main = "Generalizability Overlap",
                 xlab = "Standardized Selection Score",
                 ylab = "Density",
                 col = "blue", lwd = 2)
  graphics::lines(d_target, col = "red", lwd = 2, lty = 2)
  graphics::legend("topright",
                   legend = c("Source (compliers)", "Target population"),
                   col = c("blue", "red"),
                   lty = c(1, 2), lwd = 2)

  invisible(x)
}
