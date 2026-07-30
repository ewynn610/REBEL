# Fit linear mixed models on cell level of pseudo-bulk data

Fit linear mixed models on cell level of pseudo-bulk data

## Usage

``` r
rebelLMM(
  fixedEffects,
  normalizedCounts = NULL,
  colData = NULL,
  pseudoBulk,
  subjectVariable,
  sampleVariable = NULL,
  parallel = FALSE,
  cores = 2,
  outputFits = FALSE
)
```

## Arguments

- fixedEffects:

  One sided linear formula describing the fixed-effects on the right of
  the ~ operator. Terms should be separated by + operators. Terms should
  be variables in `colData`.

- normalizedCounts:

  Optional matrix-like object containing normalized counts. Used in
  conjunction with `colData` argument in place of object. Matrix should
  contain genes in the rows and cells (cell-level data) or samples
  (pseudo-bulk data) in the columns.

- colData:

  Optional dataframe object containing cell/sample level meta data which
  is used with `normalizedCounts` argument in place of object. Row names
  of dataframe should match column names of `normalizedCounts`.

- pseudoBulk:

  Logical value indicating whether a pseudo-bulk or cell-level analysis
  is being performed.

- subjectVariable:

  String denoting name of subject identifier in the `colData` to be used
  as the subject-level random effect in the LMM model.

- sampleVariable:

  String denoting name of sample identifier variable in the `colData` to
  be used as the sample-level random effect in the LMM model. For
  pseudo-bulk data where no sample-level random effect is used, should
  be `NULL`.

- parallel:

  Logical value indicating whether to use parallelization via
  `mclapply`.

- outputFits:

  Logical value indicating whether or not to include fit objects from
  [`lmerTest`](https://rdrr.io/pkg/lmerTest/man/lmerTest-package.html)
  in the output. Only necessary if user would like to inspect elements
  of the object. May use a large amount of memory if TRUE.

- nCores:

  Number of cores to use if `parallel` is `TRUE`.
