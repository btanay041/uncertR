.new_unc_scores <- function(df, metric) {
  structure(list(scores = df, metric = metric), class = "unc_scores")
}

#' Per-sample uncertainty score from OOB predictions
#'
#' @param fit A `boot_unc` object.
#' @param metric `"sd"`, `"ci_width"` (percentile interval), or `"iqr"`.
#' @param level Interval level for `"ci_width"`.
#' @param min_oob Samples with fewer OOB predictions get `NA`.
#' @return An object of class `unc_scores` holding a data frame with
#'   `pred_mean`, `uncertainty`, `y`, `abs_error`, and `n_oob`.
#' @export
uncertainty_score <- function(fit, metric = c("ci_width", "sd", "iqr"),
                              level = 0.95, min_oob = 10) {
  stopifnot(inherits(fit, "boot_unc"))
  metric <- match.arg(metric)
  o <- fit$oob
  n_oob <- rowSums(!is.na(o))
  f <- switch(metric,
    sd = function(v) stats::sd(v, na.rm = TRUE),
    iqr = function(v) stats::IQR(v, na.rm = TRUE),
    ci_width = function(v) {
      q <- stats::quantile(v, c((1 - level) / 2, 1 - (1 - level) / 2),
                           na.rm = TRUE, names = FALSE)
      q[2] - q[1]
    })
  u <- apply(o, 1, f)
  u[n_oob < min_oob] <- NA_real_
  pm <- rowMeans(o, na.rm = TRUE)
  df <- data.frame(sample = rownames(o), pred_mean = pm, uncertainty = u,
                   y = unname(fit$y), abs_error = abs(pm - fit$y),
                   n_oob = n_oob, row.names = rownames(o))
  .new_unc_scores(df, metric)
}

#' @export
print.unc_scores <- function(x, ...) {
  cat("<unc_scores> metric =", x$metric, "|", nrow(x$scores), "samples\n")
  print(summary(x$scores$uncertainty))
  invisible(x)
}

#' @export
plot.unc_scores <- function(x, ...) {
  d <- x$scores
  graphics::plot(d$uncertainty, d$abs_error, pch = 19, col = "#00000066",
                 xlab = paste("Uncertainty (", x$metric, ")"),
                 ylab = "Absolute OOB error", ...)
  graphics::abline(stats::lm(abs_error ~ uncertainty, d), col = "red", lwd = 2)
  invisible(x)
}

#' Does the uncertainty score track real error?
#'
#' Run this *before* interpreting drivers. If uncertainty is not associated
#' with absolute out-of-bag error, explaining it is not meaningful.
#'
#' @param unc An `unc_scores` object.
#' @return List with Spearman correlation and mean error by uncertainty quartile.
#' @export
validate_uncertainty <- function(unc) {
  d <- unc$scores[stats::complete.cases(unc$scores), ]
  ct <- suppressWarnings(stats::cor.test(d$uncertainty, d$abs_error,
                                         method = "spearman"))
  q <- cut(d$uncertainty,
           stats::quantile(d$uncertainty, 0:4 / 4), include.lowest = TRUE,
           labels = paste0("Q", 1:4))
  by_q <- stats::aggregate(abs_error ~ q, d, mean)
  list(spearman_rho = unname(ct$estimate), p_value = ct$p.value,
       error_by_quartile = by_q)
}
