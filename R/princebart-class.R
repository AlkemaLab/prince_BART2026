#' @name prince_bart-class
#' @title S3 Methods for prince_bart Objects
#'
#' @description Print, summary, and coef methods for objects returned
#' by \code{prince_BART}. Behavior is specialized by modality: binary fits
#' report principal strata and treatment effects for compliers, while ordinal
#' fits report contrasts for the monotone compliance group and stratified by
#' baseline uptake level.
#'
#' @param x,object A prince_bart object (either prince_bart_binary or prince_bart_ordinal).
#' @param type Character; type of estimand: \code{"mixed"} (default) for mixed
#'   estimands averaging over principal stratum uncertainty, or \code{"sample"}
#'   for sample estimands using imputed principal strata (experimental).
#'   For ordinal fits, this parameter is accepted but may be ignored depending
#'   on implementation.
#' @param treated_only Logical; if \code{TRUE}, compute treatment effect only
#'   among the treated units (ATT). Default is \code{FALSE} (ATE). For ordinal
#'   fits, this parameter is accepted but may be ignored or handled differently.
#' @param adaptive_levels Logical; for ordinal fits, if \code{TRUE} (default),
#'   show level-specific effects up to a data-adaptive threshold based on
#'   cumulative affected-unit mass.
#' @param cumulative_mass Numeric in (0, 1]; for ordinal fits with adaptive
#'   grouping, levels are shown individually until this cumulative mass is
#'   reached. Default is \code{0.80}.
#' @param level_threshold Optional integer threshold \code{K} for ordinal fits.
#'   If supplied, levels \code{1..K} are shown individually and higher levels
#'   are pooled as \code{"Level > K"}. This overrides adaptive grouping.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' \code{print}: Invisibly returns the object.
#' \code{summary}: Prints and invisibly returns a list with posterior summaries
#'   (including uncertainty and diagnostics) appropriate to the modality.
#' \code{coef}: Named numeric vector of posterior mean estimates for the key
#'   estimands.
#'
#' @examples
#' \dontrun{
#' # Binary fit
#' fit_binary <- prince_BART(Y ~ X | Z | W, data = df)
#' summary(fit_binary)
#' coef(fit_binary)
#'
#' # Ordinal fit
#' fit_ordinal <- prince_BART(Y ~ X | Z | W, data = df, uptake_type = "ordinal")
#' summary(fit_ordinal)
#' coef(fit_ordinal)
#' }
NULL


#' @rdname prince_bart-class
#' @export
print.prince_bart <- function(x, ...) {
  cat("Principal Stratification BART Fit\n")
  cat("---------------------------------\n")

  if (length(x$chains) > 0) {
    first_draw <- x$chains[[1]][[1]]
    cat("Chains:      ", length(x$chains), "\n")
    cat("Iterations:  ", nrow(first_draw), "\n")
    cat("Units:       ", ncol(first_draw), "\n")
  }

  if (any(vapply(x$chains, function(chain) !is.null(chain$trees), logical(1)))) {
    cat("Trees:        saved\n")
  }

  uptake_label <- if (inherits(x, "prince_bart_binary")) {
    "Binary uptake (3 principal strata: complier, never-taker, always-taker)"
  } else {
    "Ordinal/count uptake"
  }
  cat("Uptake model: ", uptake_label, "\n")

  cat("\nUse summary() or coef() to extract contrasts and treatment effects.\n")
  invisible(x)
}


#' @rdname prince_bart-class
#' @export
summary.prince_bart <- function(object
  , type = c("mixed", "sample")
  , treated_only = FALSE
  , adaptive_levels = TRUE
  , cumulative_mass = 0.80
  , level_threshold = NULL
  , ...
) {
  if (inherits(object, "prince_bart_binary")) {
    summary_prince_bart_binary(object, type = type, treated_only = treated_only)
  } else if (inherits(object, "prince_bart_ordinal")) {
    summary_prince_bart_ordinal(
      object,
      type = type,
      treated_only = treated_only,
      adaptive_levels = adaptive_levels,
      cumulative_mass = cumulative_mass,
      level_threshold = level_threshold
    )
  } else {
    stop("Unknown prince_bart class: ", paste(class(object), collapse = ", "))
  }
}


#' @rdname prince_bart-class
#' @export
coef.prince_bart <- function(object
  , type = c("mixed", "sample")
  , treated_only = FALSE
  , adaptive_levels = TRUE
  , cumulative_mass = 0.80
  , level_threshold = NULL
  , ...
) {
  if (inherits(object, "prince_bart_binary")) {
    coef_prince_bart_binary(object, type = type, treated_only = treated_only)
  } else if (inherits(object, "prince_bart_ordinal")) {
    coef_prince_bart_ordinal(
      object,
      type = type,
      treated_only = treated_only,
      adaptive_levels = adaptive_levels,
      cumulative_mass = cumulative_mass,
      level_threshold = level_threshold
    )
  } else {
    stop("Unknown prince_bart class: ", paste(class(object), collapse = ", "))
  }
}


# =============================================================================
# Internal helper functions for modality-specific summary/coef
# =============================================================================


#' @keywords internal
.posterior_mean_vector <- function(summary_tbl) {
  summary_df <- as.data.frame(summary_tbl)
  if (nrow(summary_df) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }
  if (!all(c("variable", "mean") %in% names(summary_df))) {
    stop("Expected summary to contain 'variable' and 'mean' columns")
  }
  stats::setNames(
    as.numeric(summary_df$mean), as.character(summary_df$variable)
  )
}

#' @keywords internal
summary_prince_bart_binary <- function(object
  , type = c("mixed", "sample"), treated_only = FALSE
) {
  type <- match.arg(type)

  cat("Principal Stratification BART Summary (Binary Uptake)\n")
  cat("=====================================================\n\n")

  # Compute treatment effect using existing binary estimand functions
  effect_label <- if (treated_only) "ATT" else "ATE"
  type_label <- if (type == "mixed") "Mixed" else "Sample"

  if (type == "mixed") {
    if (treated_only) {
      effect <- matt_c(object)
    } else {
      effect <- mate_c(object)
    }
  } else {
    if (treated_only) {
      effect <- satt_c(object)
    } else {
      effect <- sate_c(object)
    }
  }

  # effect[[1]] = strata distribution, effect[[2]] = outcomes & treatment effect
  strata_summary <- effect[[1]]
  effect_summary <- effect[[2]]

  cat("Principal Strata Distribution:\n")
  print(as.data.frame(strata_summary))
  cat("\n")

  cat(sprintf("%s %s for Compliers:\n", type_label, effect_label))
  effect_summary_df <- as.data.frame(effect_summary)
  print(effect_summary_df)

  invisible(list(
    strata = strata_summary,
    effect = effect_summary,
    type = type,
    treated_only = treated_only
  ))
}

#' @keywords internal
coef_prince_bart_binary <- function(object
  , type = c("mixed", "sample"), treated_only = FALSE
) {
  type <- match.arg(type)

  if (type == "mixed") {
    if (treated_only) {
      result <- matt_c(object)[[2]]
    } else {
      result <- mate_c(object)[[2]]
    }
  } else {
    if (treated_only) {
      result <- satt_c(object)[[2]]
    } else {
      result <- sate_c(object)[[2]]
    }
  }

  result_df <- as.data.frame(result)
  .posterior_mean_vector(result_df)
}

#' @keywords internal
summary_prince_bart_ordinal <- function(
  object,
  type = c("mixed", "sample"),
  treated_only = FALSE,
  adaptive_levels = TRUE,
  cumulative_mass = 0.80,
  level_threshold = NULL
) {
  cat("Principal Stratification BART Summary (Ordinal Uptake)\n")
  cat("=====================================================\n\n")
  cat("Type and treated_only parameters are currently 
  interpreted in overall contrast computation.\n\n")

  results <- estimands_ordinal_mixed(
    object,
    adaptive_levels = adaptive_levels,
    cumulative_mass = cumulative_mass,
    level_threshold = level_threshold
  )

  grouping <- results$grouping

  cat("Treatment effects among affected units\n")
  cat("--------------------------------------\n\n")

  cat("Mixed ATE among affected units:\n")
  print(as.data.frame(results$overall))
  cat("\n")

  cat("By treatment level affected by the instrument:\n")
  print(as.data.frame(results$by_w0))

  cat("\n")
  if (grouping$rule == "manual") {
    cat(sprintf(
      "Levels 1..%d are shown individually by user threshold;
      higher levels are pooled when present.\n",
      grouping$threshold
    ))
  } else if (grouping$rule == "adaptive") {
    cat(sprintf(
      "Levels are shown individually until they cover %.0f%% of affected units;
       higher levels are pooled.\n",
      100 * grouping$cumulative_mass
    ))
  } else {
    cat("All observed affected levels are shown individually;
    no pooling was applied.\n")
  }

  if (any(is.na(results$overall$mean)) || any(is.na(results$by_w0$mean))) {
    cat("\nNote: NA rows indicate no posterior units 
    matched that contrast in some draws/strata.\n")
  }

  invisible(list(
    overall = results$overall,
    by_w0 = results$by_w0,
    grouping = grouping
  ))
}

#' @keywords internal
coef_prince_bart_ordinal <- function(
  object,
  type = c("mixed", "sample"),
  treated_only = FALSE,
  adaptive_levels = TRUE,
  cumulative_mass = 0.80,
  level_threshold = NULL
) {
  results <- estimands_ordinal_mixed(
    object,
    adaptive_levels = adaptive_levels,
    cumulative_mass = cumulative_mass,
    level_threshold = level_threshold
  )
  c(
    .posterior_mean_vector(results$overall),
    .posterior_mean_vector(results$by_w0)
  )
}
