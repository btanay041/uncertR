.rss <- function(X, u) sum(stats::lm.fit(X, u)$residuals^2)

.partial_r2 <- function(X, u, assign) {
  full <- .rss(X, u)
  terms <- sort(unique(assign[assign > 0]))
  vapply(terms, function(t) {
    red <- .rss(X[, assign != t, drop = FALSE], u)
    (red - full) / red
  }, numeric(1))
}

#' Which covariates explain prediction uncertainty?
#'
#' Fits a joint linear model of (log) uncertainty on covariates and reports
#' the partial R^2 of each covariate with a permutation p-value and BH FDR.
#' Optionally removes the dependence of uncertainty on the predicted value
#' first, since uncertainty is often mechanically tied to how extreme a
#' prediction is.
#'
#' @param unc An `unc_scores` object.
#' @param covariates Data frame of covariates, rownames = sample IDs.
#' @param adjust_for `"prediction"` residualizes on a spline of the mean
#'   OOB prediction; `"none"` skips this.
#' @param log_transform Model log(uncertainty) (recommended; positive, skewed).
#' @param n_perm Number of permutations.
#' @param seed Optional seed.
#' @return An object of class `unc_drivers`.
#' @export
uncertainty_drivers <- function(unc, covariates,
                                adjust_for = c("prediction", "none"),
                                log_transform = TRUE, n_perm = 500,
                                seed = NULL) {
  adjust_for <- match.arg(adjust_for)
  if (!is.null(seed)) set.seed(seed)
  d <- unc$scores
  common <- intersect(rownames(d), rownames(covariates))
  if (length(common) < 20) stop("Fewer than 20 samples overlap with covariates.")
  d <- d[common, , drop = FALSE]
  cov <- as.data.frame(covariates[common, , drop = FALSE])
  ok <- stats::complete.cases(cbind(d[, c("uncertainty", "pred_mean")], cov)) &
    d$uncertainty > 0
  d <- d[ok, , drop = FALSE]; cov <- cov[ok, , drop = FALSE]
  if (ncol(cov) == 0) stop("No covariates supplied.")

  u <- if (log_transform) log(d$uncertainty) else d$uncertainty
  if (adjust_for == "prediction")
    u <- stats::residuals(stats::lm(u ~ splines::ns(d$pred_mean, df = 3)))

  # scale numerics so partial R2 / coefficients are comparable
  for (nm in names(cov)) if (is.numeric(cov[[nm]])) cov[[nm]] <- as.numeric(scale(cov[[nm]]))
  mm <- stats::model.matrix(~ ., data = cov)
  assign <- attr(mm, "assign")
  labs <- colnames(cov)[sort(unique(assign[assign > 0]))]

  obs <- .partial_r2(mm, u, assign)
  perm <- replicate(n_perm, .partial_r2(mm, sample(u), assign))
  if (is.null(dim(perm))) perm <- matrix(perm, nrow = length(obs))
  p <- (1 + rowSums(perm >= obs)) / (1 + n_perm)

  rho <- vapply(labs, function(v) {
    if (is.numeric(cov[[v]])) suppressWarnings(stats::cor(cov[[v]], u, method = "spearman")) else NA_real_
  }, numeric(1))

  tab <- data.frame(covariate = labs, partial_r2 = obs, spearman_rho = rho,
                    p_perm = p, fdr = stats::p.adjust(p, "BH"),
                    row.names = NULL)
  tab <- tab[order(-tab$partial_r2), ]
  full_r2 <- 1 - .rss(mm, u) / sum((u - mean(u))^2)

  structure(list(table = tab, model_r2 = full_r2, n = length(u),
                 adjust_for = adjust_for, metric = unc$metric),
            class = "unc_drivers")
}

#' @export
print.unc_drivers <- function(x, ...) {
  cat("<unc_drivers> n =", x$n, "| metric =", x$metric,
      "| adjusted for:", x$adjust_for, "\n")
  cat(sprintf("Joint model R^2 = %.3f\n\n", x$model_r2))
  print(utils::head(x$table, 10), row.names = FALSE, digits = 3)
  invisible(x)
}

#' @export
summary.unc_drivers <- function(object, fdr = 0.05, ...) {
  sig <- object$table[object$table$fdr <= fdr, ]
  cat(nrow(sig), "covariate(s) significant at FDR <=", fdr, "\n")
  print(sig, row.names = FALSE, digits = 3)
  invisible(sig)
}

#' @export
plot.unc_drivers <- function(x, fdr = 0.05, top = 15, ...) {
  t <- utils::head(x$table, top)
  t <- t[nrow(t):1, ]
  op <- graphics::par(mar = c(4, 10, 2, 1)); on.exit(graphics::par(op))
  graphics::barplot(t$partial_r2, names.arg = t$covariate, horiz = TRUE,
                    las = 1, col = ifelse(t$fdr <= fdr, "#c0392b", "grey70"),
                    xlab = "Partial R^2 (red = FDR significant)",
                    main = "Drivers of prediction uncertainty", ...)
  invisible(x)
}
