# REBEL: Repeated measures Empirical Bayes differential Expression analysis using Linear mixed models
## Getting started 

REBEL is an R package for analyzing cell-type specific differential expression in longitudinal and other repeated measure scRNA-seq studies. It can be installed from github using the following code:

```{r}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("https://github.com/ewynn610/REBEL", build_vignettes = TRUE)
```

Once the package is installed, to access a vignette introducing the package workflow use:
```{r}
vignette("rebelVignette")
```
