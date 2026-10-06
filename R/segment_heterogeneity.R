# This is required for R CMD check to suppress visible binding notes for NSE variables
if (getRversion() >= "2.15.1") utils::globalVariables(
  c("segment_label", "is_overall", "ci_lower", "ci_upper", "d")
)
#' Effect Heterogeneity by Segments
#'
#' Partition the sample into covariate-defined segments and summarize how
#' estimated latent-group effects vary across segments, following the
#' heterogeneity analysis in the Prince BART paper.
#'
#' @param princebart A fitted object from \code{prince_BART}.
#' @param data Optional data.frame of covariates used for segmentation. If NULL,
#'   uses stored raw covariates \code{princebart$data$X_raw}.
#' @param vars Optional character vector of covariate names to include in the
#'   segmentation model. If NULL, uses all covariates in \code{data}.
#' @param rpart_control An \code{rpart.control} object for tree fitting.
#' @param min_compliers_bucket Minimum effective latent-group count for pruning
#'   (sum of posterior mean latent-group
#' membership probabilities within a segment).
#'   Default is 100.
#' @param ci_level Credible interval level
#' (e.g., 0.9 for 90% credible intervals).
#'   Default is 0.9.
#' @param plot Logical; if TRUE, return ggplot2 objects for segment effects
#'   (and the max--min contrast plot if \code{contrast = TRUE}).
#' @param contrast Logical; if TRUE, compute posterior draws for the difference
#'   between the segments with the highest and lowest estimated mean effects.
#'
#' @details
#' This function fits a shallow \code{rpart} regression tree to posterior mean
#' conditional latent-group effects and assigns each unit to
#' a terminal-node segment. Tree fitting is 
#' weighted by each unit's posterior mean
#' latent-group membership probability, so that splits emphasize regions of the
#' covariate space with more relevant units.
#'
#' Raw covariates stored in the fitted object are used by default so that factor
#' and ordered-factor variables can be represented more interpretably in segment
#' labels.
#'
#' For each segment, the function aggregates posterior draws of
#' conditional effects using latent-group membership weights to obtain a
#' segment-specific average effect for the relevant latent group.
#' For binary fits,
#' this reproduces segment-specific mixed effects among compliers.
#' For ordinal fits,
#' this targets segment-specific mixed effects among affected units defined by
#' \eqn{W(0)-W(1)=1}. These segment-level summaries are mixed (sample-based)
#' estimands: they average conditional effects over the empirical covariate
#' distribution within each segment rather than defining new population-level
#' causal parameters.
#'
#' The main segment summaries are returned in \code{res$effects} and
#' visualized in \code{res$plot$effect} (if \code{plot = TRUE}).
#'
#' Optionally, segments with small effective latent-group counts can be
#' pruned/merged for stability. If \code{contrast = TRUE}, heterogeneity is
#' summarized by the posterior distribution of the difference between the
#' segments with the largest and smallest estimated mean effects. Numerical
#' results are in \code{res$contrast$summary} and the histogram is in
#' \code{res$plot$diff} (if \code{plot = TRUE}).
#'
#' @return A list with components:
#' \item{tree}{Fitted (and optionally pruned) rpart tree.}
#' \item{segment_data}{Input data with \code{cate} (posterior mean unit-level
#'   conditional effect), \code{w} (posterior mean latent-group membership
#'   probability), and \code{segment}.}
#' \item{effects}{Segment-level effect summaries: posterior mean, sd, credible
#'   interval bounds, and \code{p_gt0}, plus \code{n} (segment size) and
#'   \code{n_group} (estimated weighted latent-group size).}
#' \item{draws}{List of posterior draw matrices by segment.}
#' \item{contrast}{If \code{contrast = TRUE}
#' , posterior comparison of the highest- vs
#'   lowest-effect segments, containing \code{$summary} and
#' posterior draws in \code{$draws}.}
#' \item{plot}{If \code{plot = TRUE}
#' , a list of ggplot objects including \code{$effect}
#'   (segment effects) and, when \code{contrast = TRUE}
#' , \code{$diff} (difference distribution).}
#'
#' @export

segment_heterogeneity <- function(
  princebart,
  data = NULL,
  vars = NULL,
  rpart_control = rpart::rpart.control(cp = 0, maxdepth = 3L),
  min_compliers_bucket = 100,
  ci_level = 0.9,
  plot = FALSE,
  contrast = FALSE
) {
  if (!inherits(princebart, "prince_bart")) {
    stop("princebart must be a fitted prince_bart object")
  }
  if (is.null(princebart$chains) || length(princebart$chains) == 0L) {
    stop("princebart must contain posterior chains")
  }

  if (is.null(data)) {
    data <- get_fit_covariates(princebart, type = "raw")
    if (is.null(data)) {
      stop("segment_heterogeneity() requires princebart$data$X_raw; please refit with current prince_BART()")
    }
  }
  data <- as.data.frame(data)

  # Make names syntactic to safely evaluate model formulas with factor dummies.
  orig_names <- names(data)
  safe_names <- make.names(orig_names, unique = TRUE)
  names(data) <- safe_names
  name_map <- stats::setNames(safe_names, orig_names)

  n_units <- ncol(princebart$chains[[1]][[1]])
  if (nrow(data) != n_units) {
    stop("data and princebart chains must have the same number of units")
  }

  if (is.null(vars)) {
    vars <- setdiff(names(data), c("cate", "w", "segment", "e"))
  } else {
    vars <- as.character(vars)
    if (!all(vars %in% names(name_map))) {
      stop("vars contains names not found in data")
    }
    vars <- unname(name_map[vars])
  }
  if (length(vars) == 0) {
    stop("vars must contain at least one covariate name")
  }
  if (!all(vars %in% names(data))) {
    stop("vars contains names not found in data")
  }

  alpha <- (1 - ci_level) / 2

  seg_q <- extract_segment_quantities(princebart)
  cate_draws <- seg_q$cate_draws
  group_prob_draws <- seg_q$group_prob_draws

  data$cate <- Reduce(`+`, lapply(cate_draws, colMeans)) / length(cate_draws)
  data$w <- Reduce(`+`, lapply(group_prob_draws, colMeans)) /
    length(group_prob_draws)

  rpart_data <- data[, unique(c("cate", vars, "w")), drop = FALSE]
  prep <- preprocess_segmentation_vars(rpart_data, vars)
  rpart_data <- prep$data
  ordered_meta <- prep$ordered_meta
  seg_formula <- stats::reformulate(vars, response = "cate")

  tree <- rpart::rpart(
    seg_formula,
    data = rpart_data,
    weights = rpart_data$w,
    control = rpart_control
  )

  if (min_compliers_bucket > 0) {
    tree <- prune_weighted_tree(rpart_data, tree, min_compliers_bucket)
  }

  data$segment <- assign_segments(tree, rpart_data, ordered_meta)

  segment_draws <- function(include) {
    if (!any(include)) {
      return(
        matrix(NA_real_, nrow = nrow(cate_draws[[1]]),
          ncol = length(cate_draws))
      )
    }
    estimates <- lapply(seq_along(cate_draws), function(i) {
      weight <- group_prob_draws[[i]][, include, drop = FALSE]
      numer <- rowSums(cate_draws[[i]][, include, drop = FALSE] * weight)
      denom <- rowSums(weight)
      result <- numer / denom
      result[denom == 0] <- NA_real_
      result
    })
    do.call(cbind, estimates)
  }

  summarize_draws <- function(draws) {
    v <- as.vector(draws)
    v <- v[is.finite(v)]
    if (length(v) == 0) {
      return(list(mean = NA_real_, sd = NA_real_, ci90 = c(NA_real_, NA_real_),
                  p_gt0 = NA_real_))
    }
    list(
      mean = mean(v),
      sd = stats::sd(v),
      ci = stats::quantile(v, c(alpha, 1 - alpha), names = FALSE),
      p_gt0 = mean(v > 0)
    )
  }

  effects <- list()
  draws <- list()

  add_segment_result <- function(seg_name, include) {
    est <- segment_draws(include)
    sumry <- summarize_draws(est)
    effects[[seg_name]] <<- data.frame(
      segment = seg_name,
      mean = sumry$mean,
      sd = sumry$sd,
      ci_lower = sumry$ci[1],
      ci_upper = sumry$ci[2],
      p_gt0 = sumry$p_gt0,
      n = sum(include),
      n_group = sum(data$w[include]),
      stringsAsFactors = FALSE
    )
    draws[[seg_name]] <<- est
  }

  add_segment_result("overall", rep(TRUE, nrow(data)))
  for (seg in levels(data$segment)) {
    add_segment_result(seg, data$segment == seg)
  }

  effects <- do.call(rbind, effects)
  class(effects) <- c("princebart_segment_effects", class(effects))

  contrast_res <- NULL
  if (contrast) {
    seg_levels <- levels(droplevels(data$segment))
    if (length(seg_levels) < 2L) {
      contrast_res <- list(
        segment_low = NA_character_,
        segment_high = NA_character_,
        draws = numeric(0),
        summary = data.frame(
          mean = NA_real_,
          sd = NA_real_,
          ci_lower = NA_real_,
          ci_upper = NA_real_,
          p_gt0 = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    } else {
      seg_means <- tapply(data$cate, data$segment, mean, na.rm = TRUE)
      seg_min <- names(which.min(seg_means))
      seg_max <- names(which.max(seg_means))
      est_min <- segment_draws(data$segment == seg_min)
      est_max <- segment_draws(data$segment == seg_max)
      diff_draws <- as.vector(est_max - est_min)
      diff_draws <- diff_draws[is.finite(diff_draws)]
      if (length(diff_draws) == 0) {
        diff_draws <- NA_real_
      }
      contrast_res <- list(
        segment_low = seg_min,
        segment_high = seg_max,
        draws = diff_draws,
        summary = data.frame(
          mean = mean(diff_draws, na.rm = TRUE),
          sd = stats::sd(diff_draws, na.rm = TRUE),
          ci_lower = stats::quantile(
            diff_draws, alpha, na.rm = TRUE, names = FALSE
          )
          , ci_upper = stats::quantile(
            diff_draws, 1 - alpha, na.rm = TRUE, names = FALSE
          )
          , p_gt0 = mean(diff_draws > 0, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  wrap_label <- function(x, width = 40) {
    vapply(
      x,
      function(s) paste(strwrap(s, width = width), collapse = "\n"),
      character(1)
    )
  }

  plot_res <- NULL
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 is required for plot = TRUE")
    }
    plot_df <- effects
    plot_df$segment_label <- as.character(plot_df$segment)
    plot_df$is_overall <- plot_df$segment == "overall"

    effect_plot <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        y = stats::reorder(wrap_label(segment_label, 40), mean),
        x = mean,
        color = !is_overall
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, color = "grey") +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
        width = 0.5
      ) +
      ggplot2::geom_point() +
      ggplot2::labs(
        x = paste("Subgroup-specific effect size estimate ("
          , ci_level * 100, "% CI)", sep = ""
        ),
        y = "",
        title = ""
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")

    plot_res <- list(effect = effect_plot)

    if (contrast && !is.null(contrast_res)) {
      diff_df <- data.frame(d = contrast_res$draws)
      diff_plot <- ggplot2::ggplot(diff_df, ggplot2::aes(x = d)) +
        ggplot2::geom_histogram(bins = 30, fill = "grey70", color = "white") +
        ggplot2::geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
        ggplot2::labs(
          x = "Difference in effect across subgroups",
          y = "Count",
          title = ""
        ) +
        ggplot2::theme_minimal()
      plot_res$diff <- diff_plot
    }
  }

  out <- list(
    tree = tree,
    segment_data = data,
    effects = effects,
    draws = draws,
    contrast = contrast_res,
    plot = plot_res
  )
  class(out) <- c("princebart_segment", "list")
  out
}

assign_segments <- function(tree, data = NULL, ordered_meta = list()) {
  party_tree <- partykit::as.party(tree)
  list_rules <- utils::getFromNamespace(".list.rules.party", "partykit")
  rules <- list_rules(party_tree)
  rules <- simplify_and_pretty_rules(rules)
  if (!is.null(data)) {
    rules <- prettify_binary_rules(rules, data, ordered_meta = ordered_meta)
  }
  rules[is.na(rules) | trimws(rules) == ""] <- "all units"
  nodes <- stats::predict(party_tree, type = "node")
  factor(nodes, labels = rules)
}

#' @export
print.princebart_segment_effects <- function(x, digits = 3, ...) {
  print.data.frame(x, digits = digits, row.names = FALSE, ...)
  invisible(x)
}

prune_weighted_tree <- function(data, tree, min_weighted_bucket = 100) {
  if (is.null(tree$cptable) || nrow(tree$cptable) == 0) {
    return(tree)
  }
  cp_seq <- sort(unique(tree$cptable[, "CP"]))
  for (cp in cp_seq) {
    tree <- rpart::prune(tree, cp = cp)
    seg_weights <- weighted_segment_sizes(data, tree)
    if (all(seg_weights > min_weighted_bucket, na.rm = TRUE)) {
      break
    }
  }
  tree
}

weighted_segment_sizes <- function(data, tree) {
  seg <- assign_segments(tree)
  tapply(data$w, seg, sum)
}

#' @keywords internal
extract_segment_quantities <- function(princebart) {
  if (inherits(princebart, "prince_bart_binary")) {
    group_prob_draws <- lapply(princebart$chains,
      function(chain) 1 - chain$p_n - chain$p_a)
    cate_draws <- lapply(princebart$chains,
      function(chain) chain$m_y1c - chain$m_y0c)
    group_label <- "complier"
  } else if (inherits(princebart, "prince_bart_ordinal")) {
    group_prob_draws <- lapply(princebart$chains,
      function(chain) 1 * ((chain$w0 - chain$w1) == 1))
    cate_draws <- lapply(princebart$chains,
      function(chain) chain$m_y1 - chain$m_y0)
    group_label <- "affected"
  } else {
    stop(
      "princebart must inherit from prince_bart_binary or prince_bart_ordinal"
    )
  }

  if (!all(vapply(seq_along(cate_draws), function(i)
      identical(dim(cate_draws[[i]]), dim(group_prob_draws[[i]])),
      logical(1)))) {
    stop("effect and group matrices must have identical dimensions")
  }

  list(
    cate_draws = cate_draws,
    group_prob_draws = group_prob_draws,
    group_label = group_label
  )
}

simplify_and_pretty_rules <- function(rules, digits = 1) {

  pretty_num <- function(x) {
    if (is.na(x)) return(NA_character_)
    if (abs(x) >= 1000) return(paste0(round(x/1000, digits), "K"))
    formatC(x, format = "f", digits = digits)
  }

  simplify_one <- function(rule) {
    parts <- strsplit(rule, " & ", fixed = TRUE)[[1]]

    # Parse: var op value  (value numeric)
    rx <- paste0(
      "^\\s*([A-Za-z.][A-Za-z0-9._]*)",
      "\\s*(<=|>=|<|>)",
      "\\s*([-+]?[0-9]*\\.?[0-9]+)\\s*$"
    )
    m <- regexec(rx, parts, perl = TRUE)
    g <- regmatches(parts, m)

    ok  <- lengths(g) == 4
    raw <- parts[!ok]
    g   <- g[ok]
    if (length(g) == 0) return(rule)

    var <- vapply(g, "[[", "", 2)
    op  <- vapply(g, "[[", "", 3)
    val <- as.numeric(vapply(g, "[[", "", 4))

    out <- character(0)
    for (v in unique(var)) {
      idx <- var == v

      # strongest lower bound: keep the maximum among >= or >
      lb <- idx & op %in% c(">=", ">")
      if (any(lb)) {
        best <- max(val[lb])
        out <- c(out, paste(v, ">=", pretty_num(best)))
      }

      # strongest upper bound: keep the minimum among < or <=
      ub <- idx & op %in% c("<", "<=")
      if (any(ub)) {
        best <- min(val[ub])
        out <- c(out, paste(v, "<", pretty_num(best)))
      }
    }

    paste(c(out, raw), collapse = " & ")
  }

  vapply(rules, simplify_one, character(1))
}

#' @keywords internal
find_dummy_groups <- function(binary_vars, min_prefix = 2L) {
  n <- length(binary_vars)
  if (n < 2L) return(list())

  lcp2 <- function(a, b) {
    ca <- strsplit(a, "")[[1]]
    cb <- strsplit(b, "")[[1]]
    len <- min(length(ca), length(cb))
    k <- 0L
    while (k < len && ca[k + 1L] == cb[k + 1L]) k <- k + 1L
    substr(a, 1L, k)
  }

  result <- list()
  for (i in seq_len(n - 1L)) {
    for (j in seq(i + 1L, n)) {
      pref <- lcp2(binary_vars[i], binary_vars[j])
      np   <- nchar(pref)
      if (np >= min_prefix) {
        vi <- binary_vars[i]
        vj <- binary_vars[j]
        if (nchar(vi) > np && is.null(result[[vi]])) {
          result[[vi]] <- c(base = pref, level = substr(vi, np + 1L, nchar(vi)))
        }
        if (nchar(vj) > np && is.null(result[[vj]])) {
          result[[vj]] <- c(base = pref, level = substr(vj, np + 1L, nchar(vj)))
        }
      }
    }
  }
  result
}

#' @keywords internal
prettify_binary_rules <- function(rules, data, ordered_meta = list()) {
  ordered_vars <- names(ordered_meta)

  # Identify binary (0/1) columns in the segmentation data.
  is_bin <- vapply(names(data), function(nm) {
    col <- data[[nm]]
    v <- unique(suppressWarnings(as.numeric(col)))
    v <- v[is.finite(v)]
    length(v) <= 2L && all(v %in% c(0, 1))
  }, logical(1))
  binary_vars <- names(data)[is_bin]

  # Detect factor dummy groups via longest common prefix.
  dummy_map <- if (length(binary_vars) > 0L) find_dummy_groups(binary_vars) else list()

  rx <- paste0(
    "^\\s*([A-Za-z.][A-Za-z0-9._]*)",
    "\\s*(<=|>=|<|>)",
    "\\s*([-+]?[0-9]*\\.?[0-9]+)\\s*$"
  )

  parse_part <- function(p) {
    m <- regmatches(p, regexec(rx, p, perl = TRUE))[[1]]
    if (length(m) == 4L) {
      list(var = m[2L], op = m[3L], val = as.numeric(m[4L]), raw = p)
    } else {
      NULL
    }
  }

  ordered_expr_from_conds <- function(var, conds) {
    lev <- ordered_meta[[var]]
    if (is.null(lev)) return(NULL)
    nlev <- length(lev)
    if (nlev == 0L) return(NULL)

    codes <- 0:(nlev - 1L)
    ok <- vapply(codes, function(code) {
      all(vapply(conds, function(cd) {
        if (cd$op == "<") return(code < cd$val)
        if (cd$op == "<=") return(code <= cd$val)
        if (cd$op == ">") return(code > cd$val)
        code >= cd$val
      }, logical(1)))
    }, logical(1))

    allowed <- codes[ok]
    if (length(allowed) == 0L) return(NULL)
    lo <- min(allowed)
    hi <- max(allowed)
    if (lo == 0L && hi == (nlev - 1L)) return(NULL)
    if (lo == 0L && hi == 0L && nlev >= 2L) {
      return(sprintf('%s < %s', var, lev[2L]))
    }
    if (lo == 0L) return(sprintf('%s <= %s', var, lev[hi + 1L]))
    if (hi == (nlev - 1L)) return(sprintf('%s > %s', var, lev[lo]))
    if (lo == hi) return(sprintf('%s == %s', var, lev[lo + 1L]))
    paste0(
      sprintf('%s > %s', var, lev[lo]),
      " & ",
      sprintf('%s <= %s', var, lev[hi + 1L])
    )
  }

  rewrite_part <- function(var, op, val_chr) {
    if (!(var %in% binary_vars)) return(paste(var, op, val_chr))
    if (abs(as.numeric(val_chr) - 0.5) > 1e-6) return(paste(var, op, val_chr))
    is_pos <- op %in% c(">=", ">")
    grp <- dummy_map[[var]]
    if (!is.null(grp)) {
      base  <- grp[["base"]]
      level <- grp[["level"]]
      if (is_pos) sprintf('%s == %s', base, level)
      else        sprintf('%s != %s', base, level)
    } else {
      if (is_pos) var
      else        paste("not", var)
    }
  }

  rewrite_rule <- function(rule) {
    parts <- strsplit(rule, " & ", fixed = TRUE)[[1]]
    parsed <- lapply(parts, parse_part)

    ordered_rewrites <- character(0)
    consumed_idx <- rep(FALSE, length(parts))
    if (length(ordered_vars) > 0L) {
      for (v in ordered_vars) {
        idx <- which(vapply(parsed, function(x) !is.null(x) && identical(x$var, v), logical(1)))
        if (length(idx) > 0L) {
          expr <- ordered_expr_from_conds(v, parsed[idx])
          if (!is.null(expr)) {
            ordered_rewrites <- c(ordered_rewrites, expr)
            consumed_idx[idx] <- TRUE
          }
        }
      }
    }

    remaining <- vapply(seq_along(parts), function(i) {
      if (consumed_idx[i]) return(NA_character_)
      p <- parsed[[i]]
      if (is.null(p)) return(parts[i])
      rewrite_part(p$var, p$op, as.character(p$val))
    }, character(1L))
    remaining <- remaining[!is.na(remaining)]
    paste(c(ordered_rewrites, remaining), collapse = " & ")
  }

  vapply(rules, rewrite_rule, character(1L))
}

#' @keywords internal
preprocess_segmentation_vars <- function(data, vars) {
  out <- data
  ordered_meta <- list()
  for (v in vars) {
    if (is.ordered(out[[v]])) {
      ordered_meta[[v]] <- levels(out[[v]])
      out[[v]] <- as.integer(out[[v]]) - 1L
    }
  }
  list(data = out, ordered_meta = ordered_meta)
}
