# Fit REBEL Models

Fit REBEL Models

## Usage

``` r
rebelFit(
  object,
  assay = "normcounts",
  fixedEffects,
  subjectVariable,
  sampleVariable = NULL,
  normalizedCounts = NULL,
  colData = NULL,
  pseudoBulk = TRUE,
  parallel = FALSE,
  nCores = 1,
  outputFits = FALSE,
  quiet = FALSE,
  REML = TRUE
)
```

## Arguments

- object:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object (cell-level data) or
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object (pseudo-bulk data) containing normalized count values in one of
  the object assays and cell/sample level meta data in the `colData`
  slot.

- assay:

  String specifying assay in object which holds normalized count values.

- fixedEffects:

  One sided linear formula describing the fixed-effects on the right of
  the ~ operator. Terms should be separated by + operators. Terms should
  be variables in `colData`.

- subjectVariable:

  String denoting name of subject identifier in the `colData` to be used
  as the subject-level random effect in the LMM model.

- sampleVariable:

  String denoting name of sample identifier variable in the `colData` to
  be used as the sample-level random effect in the LMM model. For
  pseudo-bulk data where no sample-level random effect is used, should
  be `NULL`.

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

- parallel:

  Logical value indicating whether to use parallelization via
  `mclapply`.

- nCores:

  Number of cores to use if `parallel` is `TRUE`.

- outputFits:

  Logical value indicating whether or not to include fit objects from
  [`lmerTest`](https://rdrr.io/pkg/lmerTest/man/lmerTest-package.html)
  in the output. Only necessary if user would like to inspect elements
  of the object. May use a large amount of memory if TRUE.

- quiet:

  Logical value indicating whether messages should be printed at each
  step.

- REML:

  Logical value indicating if LMM models should be fit using REML or
  regular ML.

## Value

A
[`RebelFit-class`](https://ewynn610.github.io/REBEL/reference/RebelFit-class.md)` object`.

## Examples

``` r
## Run cell-level analysis
## Read in data (scTransform normalized data saved in normcounts slot)
data("RecAM_sim_sce")

## Just use first 10 genes
RecAM_sim_sce_fil <- RecAM_sim_sce[1:10,]

## Fit models with time, group and time/group interaction effect
cell_fit <-rebelFit(object=RecAM_sim_sce_fil,
                                   fixedEffects = ~time*group,
                                   subjectVariable ="subjectID",
                                   sampleVariable = "sampleID",
                                   pseudoBulk = FALSE)
#> [1] "Fitting LMM Models"
#> Warning: the ‘mkReTrms’ function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainer to do so.
#> Warning: the ‘findbars’ function has moved to the reformulas package. Please update your imports, or ask an upstream package maintainer to do so.
#> [1] "Getting Empirical Bayes Estimates"
#> [1] "Getting Parameters for Kenward-Rogers method"

## Run pseudo-bulk analysis
## Read in data (DESeq2 VST normalized data saved in normcounts slot)
data("RecAM_sim_pb")

## Just use first 10 genes
RecAM_sim_pb_fil <- RecAM_sim_pb[1:10,]

## Fit models with time, group and time/group interaction effect
pb_fit <- rebelFit(object=RecAM_sim_pb_fil,
                                 fixedEffects = ~time*group,
                                 subjectVariable ="subjectID",
                                 pseudoBulk = TRUE)
#> [1] "Fitting LMM Models"
#> [1] "Getting Empirical Bayes Estimates"
#> [1] "Getting Parameters for Kenward-Rogers method"
```
