# uncertR

**Bootstrap-based prediction uncertainty for transcriptomic models**

`uncertR` provides tools for estimating **per-sample prediction uncertainty** using out-of-bag (OOB) bootstrap predictions and for investigating whether biological or technical covariates are associated with that uncertainty.

The package is designed for transcriptomic drug-response modeling, but the core uncertainty workflow can be applied to other regression problems.

## What does `uncertR` do?

A prediction can be accurate on average while still being unreliable for individual samples. `uncertR` addresses this by generating bootstrap ensembles of penalized regression models and using the variation among **out-of-bag predictions** as a measure of prediction uncertainty.

The workflow has three main stages:

1. **Estimate uncertainty** from OOB bootstrap predictions. Metrics include mean confidence interval width, standard deviation, interquartile range
2. **Validate uncertainty** by testing whether samples with greater uncertainty also tend to have larger prediction errors.
3. **Identify potential drivers** of uncertainty using biological and technical covariates.

## Installation

Install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("btanay041/uncertR")
```

Then load the package:

```r
library(uncertR)
```

## Basic workflow

### 1. Generate bootstrap predictions

Provide a gene-expression matrix and a continuous outcome such as drug-response measurements.

```r
fit <- boot_uncertainty(
  expr,
  ic50,
  B = 200,
  seed = 42
)
```

`boot_uncertainty()` fits an ensemble of penalized regression models and collects **out-of-bag predictions** for each sample.

* `expr`: gene × sample expression matrix
* `ic50`: numeric response vector corresponding to the samples
* `B`: number of bootstrap replicates
* `seed`: optional random seed for reproducibility

### 2. Calculate per-sample uncertainty

```r
unc <- uncertainty_score(
  fit,
  metric = "ci_width"
)
```

Three uncertainty metrics are available:

* `ci_width` — width of the empirical prediction interval
* `sd` — standard deviation of OOB predictions
* `iqr` — interquartile range of OOB predictions

These metrics quantify **prediction dispersion** across bootstrap models; they are not prediction-error measures themselves.

### 3. Validate uncertainty

```r
validate_uncertainty(unc)
```

This evaluates whether estimated uncertainty is associated with observed prediction error.

The validation compares uncertainty with the absolute error between the mean OOB prediction and the observed outcome.

### 4. Examine potential drivers of uncertainty

`uncertR` can combine uncertainty estimates with sample-level biological or technical covariates.

```r
cov <- tumor_covariates(
  expr,
  extra = lineage_df
)

drv <- uncertainty_drivers(
  unc,
  cov,
  adjust_for = "prediction"
)
```

Results can be inspected with:

```r
plot(drv)
summary(drv)
```

`uncertainty_drivers()` evaluates the contribution of individual covariates while accounting for the specified adjustment variables.

## Data format

### Expression data

Expression data should be provided as a numeric matrix:

```text
              Sample1  Sample2  Sample3
GENE1            ...      ...      ...
GENE2            ...      ...      ...
GENE3            ...      ...      ...
```

* Rows = genes
* Columns = samples
* Row names = gene identifiers
* Column names = sample identifiers

### Response data

The response should be a numeric vector whose names correspond to the expression matrix columns:

```r
ic50 <- c(
  Sample1 = 4.2,
  Sample2 = 5.1,
  Sample3 = 3.8
)
```

The expression and response data are matched by sample identifier.

## Why out-of-bag predictions?

Each bootstrap model is trained on a resampled subset of the data. Samples that are not selected for a particular bootstrap replicate are **out-of-bag** for that replicate and can be predicted without having been used to fit that model.

`uncertR` uses these OOB predictions to estimate prediction uncertainty, reducing the optimism that can occur when uncertainty is estimated from predictions made on the same observations used to fit each model.

## Main functions

| Function                 | Purpose                                                            |
| ------------------------ | ------------------------------------------------------------------ |
| `boot_uncertainty()`     | Generate OOB bootstrap predictions                                 |
| `uncertainty_score()`    | Calculate per-sample uncertainty                                   |
| `validate_uncertainty()` | Evaluate the relationship between uncertainty and prediction error |
| `tumor_covariates()`     | Construct biological/technical covariates                          |
| `uncertainty_drivers()`  | Test covariates associated with uncertainty                        |

## Reproducibility

Bootstrap procedures are stochastic. Set a seed when reproducibility is important:

```r
fit <- boot_uncertainty(
  expr,
  ic50,
  B = 200,
  seed = 42
)
```

For final analyses, use a sufficiently large number of bootstrap replicates. Smaller values of `B` can be useful for development and testing.


## Status

`uncertR` is currently under active development. Interfaces and functionality may change as additional validation and modeling approaches are added.

