set.seed(1)
n <- 80; p <- 100
expr <- matrix(rnorm(p * n), p, n, dimnames = list(paste0("g", 1:p), paste0("s", 1:n)))
y <- as.numeric(expr[1, ] * 2 + expr[2, ] + rnorm(n, sd = 0.5))

test_that("bootstrap returns OOB matrix with NAs for in-bag samples", {
  fit <- boot_uncertainty(expr, y, B = 20, seed = 1)
  expect_s3_class(fit, "boot_unc")
  expect_equal(dim(fit$oob), c(n, 20))
  expect_true(any(is.na(fit$oob)))
})

test_that("uncertainty scores and validation run", {
  fit <- boot_uncertainty(expr, y, B = 40, seed = 1)
  for (m in c("sd", "ci_width", "iqr")) {
    unc <- uncertainty_score(fit, metric = m, min_oob = 5)
    expect_true(all(unc$scores$uncertainty >= 0, na.rm = TRUE))
  }
  v <- validate_uncertainty(unc)
  expect_true(is.numeric(v$spearman_rho))
})

test_that("drivers recovers a planted covariate", {
  n <- 200
  purity <- runif(n); noise_cov <- rnorm(n)
  ids <- paste0("s", 1:n)
  df <- data.frame(sample = ids, pred_mean = rnorm(n),
                   uncertainty = exp(-2 * purity + rnorm(n, sd = 0.3)),
                   y = rnorm(n), abs_error = runif(n), n_oob = 50,
                   row.names = ids)
  unc <- uncertR:::.new_unc_scores(df, "sd")
  cov <- data.frame(purity = purity, noise = noise_cov, row.names = ids)
  drv <- uncertainty_drivers(unc, cov, n_perm = 200, seed = 1)
  expect_equal(drv$table$covariate[1], "purity")
  expect_lt(drv$table$fdr[1], 0.05)
})

test_that("clr rows sum to ~0", {
  f <- matrix(c(.2, .3, .5, .1, .1, .8), 2, byrow = TRUE)
  expect_equal(rowSums(clr_transform(f)), c(0, 0), tolerance = 1e-8)
})
