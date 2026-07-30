# Aggregate data for pseudo-bulk analysis

Aggregate data for pseudo-bulk analysis

## Usage

``` r
rebelAggregate(
  object,
  countsAssay = "counts",
  counts = NULL,
  colData = NULL,
  sampleVariable
)
```

## Arguments

- object:

  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object containing count values in one of the object assays and cell
  level meta data in the colData slot.

- countsAssay:

  String specifying assay in object which holds count values.

- counts:

  Optional matrix-like object containing counts. Used in conjunction
  with `colData` argument in place of object. Matrix should contain
  genes in the rows and cells in the columns.

- colData:

  Optional dataframe object containing cell level meta data which is
  used with counts argument in place of object. Row names of dataframe
  should match column names of `normalizedCounts`.

- sampleVariable:

  String denoting name of sample identifier variable in the `colData`
  which counts should be aggregated over.

## Value

[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
and includes pseudo-bulk counts in the `counts` assay. All colData
variables that do not vary within the `sampleVariable` are also included
in the `colData` slot.

## Author

Elizabeth Wynn

## Examples

``` r
## Read in single cell data
data("RecAM_sim_sce")

## Aggregate into pseudo-bulk data
pb_dat <- rebelAggregate(object=RecAM_sim_sce, countsAssay="counts", sampleVariable="sampleID")
```
