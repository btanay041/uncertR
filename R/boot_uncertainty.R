#' Bootstrap ensemble with out-of-bag predictions
#'
#' Resamples samples with replacement `B` times, fits an elastic-net model on
#' each bootstrap sample, and predicts only the samples that were *not* drawn
#' (out-of-bag, OOB). Every sample therefore accumulates a distribution of
#' honest (non-in-sample) predictions.
#'
#' @param expr Genes x samples matrix (bulk expression, e.g. log2 TPM), or a
#'   `SummarizedExperiment` (first assay is used).
#' @param y Numeric outcome vector (e.g. log IC50 or AUC), one per sample.
#' @param B Number of bootstrap replicates.
#' @param alpha glmnet mixing parameter (1 = lasso, 0 = ridge).
#' @param lambda_mode `"global"` selects lambda once by CV on all data (fast;
#'   slight optimism because OOB samples informed lambda). `"per_boot"` runs
#'   inner CV in every replicate (slower, cleaner).
#' @param nfolds Folds for cross-validation.
#' @param seed Optional integer for reproducibility.
#' @return An object of class `boot_unc`.
#' @export
boot_uncertainty <- function(expr, y, B = 200, alpha = 0.5,
                             lambda_mode = c("global", "per_boot"),
                             nfolds = 5, seed = NULL) {
  lambda_mode <- match.arg(lambda_mode)
  if (inherits(expr, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
      stop("Install SummarizedExperiment to pass SE objects.")
    expr <- SummarizedExperiment::assay(expr)
  }
  expr <- as.matrix(expr)
  if (!is.numeric(y)) stop("`y` must be numeric (continuous drug response).")
  if (ncol(expr) != length(y))
    stop("ncol(expr) must equal length(y) (expr is genes x samples).")
  keep <- !is.na(y)
  x <- t(expr[, keep, drop = FALSE])
  y <- y[keep]
  ids <- colnames(expr)[keep]
  if (is.null(ids)) ids <- paste0("S", seq_along(y))
  n <- nrow(x)

  if (!is.null(seed)) set.seed(seed)
  lam <- NULL
  if (lambda_mode == "global")
    lam <- glmnet::cv.glmnet(x, y, alpha = alpha, nfolds = nfolds)$lambda.min

  one <- function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    oob <- setdiff(seq_len(n), idx)
    pred <- rep(NA_real_, n)
    if (length(oob) == 0) return(pred)
    if (is.null(lam)) {
      fit <- glmnet::cv.glmnet(x[idx, , drop = FALSE], y[idx],
                               alpha = alpha, nfolds = nfolds)
      pred[oob] <- as.numeric(stats::predict(fit, x[oob, , drop = FALSE],
                                             s = "lambda.min"))
    } else {
      fit <- glmnet::glmnet(x[idx, , drop = FALSE], y[idx],
                            alpha = alpha, lambda = lam)
      pred[oob] <- as.numeric(stats::predict(fit, x[oob, , drop = FALSE]))
    }
    pred
  }
  res <- future.apply::future_lapply(seq_len(B), one, future.seed = TRUE)
  oob <- do.call(cbind, res)
  rownames(oob) <- ids

  structure(list(oob = oob, y = stats::setNames(y, ids), B = B,
                 params = list(alpha = alpha, lambda_mode = lambda_mode,
                               nfolds = nfolds, seed = seed),
                 call = match.call()),
            class = "boot_unc")
}

#' @export
print.boot_unc <- function(x, ...) {
  cat("<boot_unc>", nrow(x$oob), "samples,", x$B, "bootstrap replicates\n")
  cat("Median OOB predictions per sample:",
      stats::median(rowSums(!is.na(x$oob))), "\n")
  invisible(x)
}

#' @export
summary.boot_unc <- function(object, ...) {
  n_oob <- rowSums(!is.na(object$oob))
  cat("OOB predictions per sample:\n")
  print(summary(n_oob))
  invisible(n_oob)
}
