# REBEL: Repeated Measures Empirical Bayes Differential Expression Analysis Using Linear Mixed Models

**REBEL** is an R package for cell-type specific differential expression
analysis in repeated measures single-cell RNA-sequencing (scRNA-seq)
studies.

The package applies linear mixed models and empirical Bayes moderation
to identify differentially expressed genes while accounting for
within-subject correlation arising from paired, longitudinal, and other
repeated measures study designs.

**Website:** [REBEL website](https://ewynn610.github.io/REBEL/)

**Vignette:** - [Introduction to
REBEL](https://ewynn610.github.io/REBEL/articles/rebelVignette.html)

## Getting Started

### Install Dependencies

``` r

## CRAN packages
install.packages(
  c(
    "dplyr",
    "lme4",
    "lmerTest",
    "limma",
    "Matrix",
    "pbapply",
    "pbkrtest"
  )
)

## Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  c(
    "DESeq2",
    "SingleCellExperiment",
    "SummarizedExperiment"
  )
)
```

### Install REBEL from GitHub

**Note:** The package vignettes are available on the package website. If
you’d prefer to build the vignettes locally during installation, you can
set `build_vignettes = TRUE`, though this may take longer.

`{r} if (!requireNamespace("devtools", quietly = TRUE)) { install.packages("devtools") } devtools::install_github("https://github.com/ewynn610/REBEL", build_vignettes = FALSE)`
