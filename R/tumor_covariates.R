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

  if (!is.null(methods) && length(methods) > 0) {

    if (!requireNamespace("immunedeconv", quietly = TRUE)) {
      stop(
        "The 'immunedeconv' package is required when 'methods' is specified. ",
        "Install it from GitHub with ",
        "remotes::install_github('omnideconv/immunedeconv').",
        call. = FALSE
      )
    }

    deconvolute <- getExportedValue(
      "immunedeconv",
      "deconvolute"
    )

    for (m in methods) {
      res <- deconvolute(expr, method = m)

      mat <- as.matrix(res[, -1, drop = FALSE])
      rownames(mat) <- paste(m, res[[1]], sep = "_")

      out <- cbind(
        out,
        t(mat)[rownames(out), , drop = FALSE]
      )
    }
  }

  if (!is.null(extra)) {
    extra <- as.data.frame(extra)
    out <- cbind(
      out,
      extra[rownames(out), , drop = FALSE]
    )
  }

  out
}
