# uncertR

Bootstrap-based prediction uncertainty for transcriptomic drug-response models,
and a test of which biological/technical covariates drive it.

```r
fit <- boot_uncertainty(expr, ic50, B = 200)          # OOB bootstrap, glmnet
unc <- uncertainty_score(fit, metric = "ci_width")
validate_uncertainty(unc)                             # does uncertainty track error?
cov <- tumor_covariates(expr, extra = lineage_df)
drv <- uncertainty_drivers(unc, cov, adjust_for = "prediction")
plot(drv); summary(drv)
```

## Roadmap
- [ ] Survival outcomes (Cox / `Surv`)
- [ ] Alternative backends (ranger, xgboost)
- [ ] Conformal prediction intervals as a complementary score
- [ ] Grouped importance / elastic-net drivers for collinear covariates
- [ ] Benchmarks on GDSC / CCLE / TCGA, vignette, pkgdown, BiocCheck
