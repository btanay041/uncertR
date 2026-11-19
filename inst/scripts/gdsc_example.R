## Worked example: GDSC2 expression (RMA, log2) + LN_IC50 for one drug.
## ADAPT the file paths / column names marked with <<< to match your files.

library(uncertR)
library(future)
plan(multisession, workers = 4)

## ---- 1. Load -------------------------------------------------------------
expr_raw <- readRDS("GDSC2_Expr (RMA Normalized and Log Transformed).rds")    
ic50_raw <- readRDS("GDSC2_Res.rds")                     


sym <- rownames(expr_raw)
mat <- as.matrix(expr_raw[, -1])
colnames(mat) <- sub("^DATA\\.", "", colnames(mat))
keep <- !is.na(sym) & sym != ""
mat <- mat[keep, ]; sym <- sym[keep]
# collapse duplicated symbols by averaging
expr <- rowsum(mat, sym) / as.vector(table(sym)[rownames(rowsum(mat, sym))])

## ---- 2. One drug at a time ----------------------------------------------
drug <- "Camptothecin_1003"

# Extract IC50 values for the selected drug
y <- ic50_raw[, drug]

# Cell-line IDs are the row names
names(y) <- rownames(ic50_raw)

# Keep only cell lines present in both datasets
common <- intersect(names(y), colnames(expr))

y <- y[common]
expr <- expr[, common, drop = FALSE]

message(
  length(y),
  " cell lines with expression + IC50 for ",
  drug
)

# Unsupervised variance filter: top 5000 genes
v <- apply(expr, 1, stats::var)

keep <- order(v, decreasing = TRUE)[
  seq_len(min(5000, nrow(expr)))
]

expr <- expr[keep, , drop = FALSE]

## ---- 3. Bootstrap uncertainty -------------------------------------------
fit <- boot_uncertainty(expr, y, B = 200, alpha = 0.5, seed = 42)
unc <- uncertainty_score(fit, metric = "sd")
plot(unc)
validate_uncertainty(unc)   # do NOT interpret drivers if this is weak

## ---- 4. Covariates -------------------------------------------------------
# Cell line annotations
ann <- read.csv("Cell_Lines_Details.csv")
rownames(ann) <- as.character(ann$COSMIC.identifier)

extra <- data.frame(
  tissue = collapse_rare_levels(
    ann[common, "GDSC.Tissue.descriptor.1"],
    15
  ),
  growth_properties = ann[common, "Growth.Properties"],
  row.names = common,
  stringsAsFactors = FALSE
)

# Convert categorical covariates to factors
extra$tissue <- factor(extra$tissue)
extra$growth_properties <- factor(extra$growth_properties)


# Expression-derived biological covariates

gene_set_score <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(expr)))
  }
  
  z <- t(scale(t(expr[genes, , drop = FALSE])))
  
  colMeans(z, na.rm = TRUE)
}


# Proliferation
extra$proliferation <- gene_set_score(
  expr,
  c("MKI67", "PCNA", "TOP2A", "MCM2", "CCNB1")
)

# DNA repair
extra$dna_repair <- gene_set_score(
  expr,
  c("BRCA1", "BRCA2", "RAD51", "ATM", "XRCC1")
)

# Apoptosis
extra$apoptosis <- gene_set_score(
  expr,
  c("BAX", "BAK1", "BCL2", "CASP3", "CASP8")
)

# Hypoxia
extra$hypoxia <- gene_set_score(
  expr,
  c("HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA")
)

# EMT
extra$emt <- gene_set_score(
  expr,
  c("VIM", "CDH2", "SNAI1", "SNAI2", "TWIST1")
)


# Build covariate matrix
cov <- tumor_covariates(expr,extra = extra,technical = FALSE)
## ---- 5. Drivers ----------------------------------------------------------
drv <- uncertainty_drivers(unc, cov, adjust_for = "prediction", n_perm = 1000)
print(drv); summary(drv); plot(drv)
