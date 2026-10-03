#' Assemble covariates that may drive uncertainty
#'
#' Combines technical covariates computed from `expr`, optional deconvolution
#' scores (patient tumors only; not meaningful for cell lines), and any
#' user-supplied covariates (e.g. tissue lineage, growth rate, mutation status).
#'
#' @param expr Genes x samples matrix (HGNC symbols, TPM-like scale for
#'   deconvolution).
#' @param methods Character vector of `immunedeconv` methods
#'   (e.g. `"estimate"`, `"quantiseq"`); `NULL` to skip.
#' @param extra Optional data frame of extra covariates; rownames = samples.
#' @param technical If `TRUE`, add log library size and number of detected
#'   genes. Only meaningful for raw/TPM RNA-seq counts; set `FALSE` for
#'   RMA-normalized microarray data such as GDSC.
#' @return Data frame, one row per sample.
#' @export
tumor_covariates <- function(expr, methods = NULL, extra = NULL,
                             technical = TRUE) {
  expr <- as.matrix(expr)
  out <- data.frame(row.names = colnames(expr))
  if (technical) {
    out$log_lib_size <- log10(colSums(expr) + 1)
    out$n_detected <- colSums(expr > 0)
  }
  for (m in methods) {
    if (!requireNamespace("immunedeconv", quietly = TRUE))
      stop("Install immunedeconv to use deconvolution methods.")
    res <- immunedeconv::deconvolute(expr, method = m)
    mat <- as.matrix(res[, -1])
    rownames(mat) <- paste(m, res[[1]], sep = "_")
    out <- cbind(out, t(mat)[rownames(out), , drop = FALSE])
  }
  if (!is.null(extra)) {
    extra <- as.data.frame(extra)
    out <- cbind(out, extra[rownames(out), , drop = FALSE])
  }
  out
}

#' Centered log-ratio transform for compositional fractions
#'
#' Cell-type fractions sum to a constant, so use CLR before regression.
#' @param fractions Samples x cell types matrix of non-negative fractions.
#' @param pseudo Pseudocount replacing zeros.
#' @export
clr_transform <- function(fractions, pseudo = 1e-3) {
  f <- as.matrix(fractions) + pseudo
  l <- log(f)
  l - rowMeans(l)
}

#' Collapse rare factor levels
#'
#' Tissue/lineage labels often have many tiny groups that make the driver
#' model unstable. Levels with fewer than `min_n` samples become `"Other"`.
#' @param x Character or factor vector.
#' @param min_n Minimum group size.
#' @export
collapse_rare_levels <- function(x, min_n = 15) {
  x <- as.character(x)
  tab <- table(x)
  x[x %in% names(tab)[tab < min_n]] <- "Other"
  factor(x)
}
