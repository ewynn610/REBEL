# Estimate parameters for Kenward Rogers Degrees of freedom

Estimate parameters for Kenward Rogers Degrees of freedom

## Usage

``` r
rebelKRParams(RebelFitObj, parallel = FALSE, nCores = 1)
```

## Arguments

- RebelFitObj:

  Object created by running rebelLMM and rebelEmpBayes

- parallel:

  Logical value indicating whether to use parallelization via
  `mclapply`.

- nCores:

  Number of cores to use if `parallel` is `TRUE`.

## Value

A
[`RebelFit-class`](https://ewynn610.github.io/REBEL/reference/RebelFit-class.md)` object`.
