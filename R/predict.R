#' Predict from Saved Trees
#'
#' Generate predictions for new data using saved tree structures from
#' a prince_bart fit.
#'
#' @param trees A data.frame of tree structures from a prince_bart fit.
#' @param newdata A matrix of covariates for prediction (can be unscaled if
#'   scaling is provided), or a list of matrices (one per posterior sample).
#' @param scaling Optional list with `center` and `scale` named vectors for
#'   standardizing newdata before prediction. Typically from `prince_bart_fit$scaling`.
#'   If NULL, newdata is assumed to already be scaled.
#' @param n_cores Number of cores for parallel prediction (default: 1).
#'
#' @return A matrix of predicted probabilities with rows for observations
#'   and columns for posterior samples.
#'
#' @export
predict_trees <- function(trees, newdata, scaling = NULL, n_cores = 1) {
  trees <- as.data.frame(trees)
  samples <- unique(trees$sample)

  # Allow a single matrix or a list of matrices, one per sample
  x_list <- if (is.list(newdata) && !is.data.frame(newdata)) {
    newdata
  } else {
    rep(list(as.matrix(newdata)), length(samples))
  }
  stopifnot(length(x_list) == length(samples))

  # Apply scaling if provided
  if (!is.null(scaling) && !is.null(scaling$center) && !is.null(scaling$scale)) {
    x_list <- lapply(x_list, function(xmat) {
      xmat <- as.matrix(xmat)
      # Get column names excluding 'e' (propensity is not scaled)
      x_cols <- setdiff(colnames(xmat), "e")
      for (v in x_cols) {
        if (v %in% names(scaling$center) && v %in% names(scaling$scale)) {
          xmat[, v] <- (xmat[, v] - scaling$center[v]) / scaling$scale[v]
        }
      }
      xmat
    })
  }

  preds <- mapply(
    samples, x_list,
    SIMPLIFY = FALSE, USE.NAMES = FALSE,
    FUN = function(sidx, xmat) {
      xmat <- as.matrix(xmat)
      n_trees <- max(trees$tree)
      these <- trees[trees$sample == sidx, ]
      pre_mat <- sapply(seq_len(n_trees), function(i) {
        get_predictions_for_tree(these[these$tree == i, ], xmat)
      })
      apply(pre_mat, 1, sum)
    }
  )

  stats::pnorm(simplify2array(preds))
}
#' Get Predictions for a Single Tree
#'
#' Recursive function to traverse a single BART tree and generate predictions.
#'
#' @param tree A data.frame representing a single tree structure.
#' @param x A matrix of covariates.
#'
#' @return A numeric vector of predictions for each row in x.
#'
#' @keywords internal
get_predictions_for_tree <- function(tree, x) {
  tree <- as.data.frame(tree)
  if (nrow(tree) == 0) {
    stop("Encountered empty tree while traversing predictions")
  }

  predictions <- rep(NA_real_, nrow(x))

  get_predictions_recursive <- function(tree, indices) {
    if (tree$var[1] == -1) {
      predictions[indices] <<- tree$value[1]
      return(1)
    }

    split_var <- tree$var[1]
    if (!is.numeric(split_var) || is.na(split_var) || split_var < 1 || split_var > ncol(x)) {
      stop(
        "Tree split variable index out of bounds: ", split_var,
        " (ncol(x)=", ncol(x), ")."
      )
    }

    goes_left <- x[indices, split_var] <= tree$value[1]
    head_left <- tree[-1, ]
    n_nodes_left <- get_predictions_recursive(head_left, indices[goes_left])

    head_right <- tree[seq.int(2 + n_nodes_left, nrow(tree)), ]
    n_nodes_right <- get_predictions_recursive(head_right, indices[!goes_left])

    return(1 + n_nodes_left + n_nodes_right)
  }

  get_predictions_recursive(tree, seq_len(nrow(x)))
  predictions
}

# Predict from one saved posterior sample's trees.
predict_one_sample <- function(trees, x) {
  trees <- as.data.frame(trees)
  if (nrow(trees) == 0) {
    stop("No trees available for this posterior sample/model component")
  }
  n_trees <- max(trees$tree)

  preds <- sapply(seq_len(n_trees), function(i) {
    one_tree <- trees[trees$tree == i, ]
    if (nrow(one_tree) == 0) {
      stop("Missing tree id ", i, " in posterior sample component")
    }
    get_predictions_for_tree(one_tree, x)
  })

  stats::pnorm(rowSums(preds))
}

predict_one_sample_raw <- function(trees, x) {
  trees <- as.data.frame(trees)
  if (nrow(trees) == 0) {
    stop("No trees available for raw posterior sample/model component")
  }
  n_trees <- max(trees$tree)

  preds <- sapply(seq_len(n_trees), function(i) {
    one_tree <- trees[trees$tree == i, ]
    if (nrow(one_tree) == 0) {
      stop("Missing tree id ", i, " in raw posterior sample component")
    }
    get_predictions_for_tree(one_tree, x)
  })

  rowSums(preds)
}
