#' Extract a binary complier treatment effect
#'
#' Compute a posterior summary for one of four complier estimands from a
#' binary-uptake \code{prince_BART()} fit.
#'
#' @param fit A fitted binary-uptake \code{prince_bart} object.
#' @param value One of \code{"mate_c"}, \code{"matt_c"}, \code{"sate_c"},
#'   or \code{"satt_c"}.
#' @param induce_residual_corr Logical; for sample estimands, use correlated
#'   potential-outcome imputations. Defaults to \code{FALSE}. Ignored for mixed
#'   estimands.
#'
#' @details
#' Let \eqn{C_i} indicate that unit \eqn{i} is a complier,
#' \eqn{Z_i} be its assignment, \eqn{p_{ci}} its posterior probability of
#' compliance, and \eqn{\mu_{zi}^c = E[Y_i(z) \mid C_i = 1, X_i]}.
#' Each expression below is evaluated separately at every posterior draw:
#' \itemize{
#'   \item \code{"mate_c"}: mixed average treatment effect for compliers,
#'     \eqn{\sum_i p_{ci}(\mu_{1i}^c-\mu_{0i}^c) / \sum_i p_{ci}}.
#'   \item \code{"matt_c"}: mixed average treatment effect for treated
#'     compliers, \eqn{\sum_{i:Z_i=1} p_{ci}(\mu_{1i}^c-\mu_{0i}^c) /
#'     \sum_{i:Z_i=1} p_{ci}}.
#'   \item \code{"sate_c"}: sample average treatment effect for compliers,
#'     \eqn{\sum_{i:C_i=1} [Y_i(1)-Y_i(0)] / \sum_i C_i}, using imputed
#'     potential outcomes and sampled compliance classes.
#'   \item \code{"satt_c"}: sample average treatment effect for treated
#'     compliers, \eqn{\sum_{i:Z_i=1,C_i=1} [Y_i(1)-Y_i(0)] /
#'     \sum_{i:Z_i=1} C_i}, using imputed potential outcomes and sampled
#'     compliance classes.
#' }
#'
#' @return A list of two posterior summary tables. The first summarizes
#'   principal stratum shares; the second summarizes potential-outcome means
#'   and the requested treatment effect. Rows correspond to quantities and
#'   columns include posterior mean, median, standard deviation, and quantiles.
#'
#' @examples
#' \dontrun{
#' fit <- prince_BART(Y ~ x1 | Z | W, data = binary_data)
#' extract(fit, "mate_c")
#' extract(fit, "satt_c")
#' }
#' @export
extract <- function(fit, value, induce_residual_corr = FALSE) {
  choices <- c("mate_c", "matt_c", "sate_c", "satt_c")
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% choices) {
    stop("value must be one of: ", paste(choices, collapse = ", "),
      call. = FALSE)
  }
  if (!inherits(fit, "prince_bart_binary")) {
    stop("extract() currently requires a binary-uptake prince_bart fit",
      call. = FALSE)
  }
  if (!is.logical(induce_residual_corr) ||
      length(induce_residual_corr) != 1L || is.na(induce_residual_corr)) {
    stop("induce_residual_corr must be TRUE or FALSE", call. = FALSE)
  }

  switch(value,
    mate_c = mate_c(fit),
    matt_c = matt_c(fit),
    sate_c = sate_c(fit, induce_residual_corr),
    satt_c = satt_c(fit, induce_residual_corr)
  )
}
