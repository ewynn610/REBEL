# Example Simulated Repeated Measures scRNA-Seq Data

Sample
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object with gene expression data simulated using the `RESCUE` R package.

## Usage

``` r
RecAM_sim_sce
```

## Format

[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object with the following data in object slots:

- `counts`:

  Matrix of raw counts with genes represented by rows and cells
  represented by columns.

- `normcounts`:

  Data normalized using the variance stabilizing transformation from the
  `sctransform` package.

- `colData`:

  sampleID

  :   Sample identifier

  subjectID

  :   Subject identifier

  time

  :   Timepoint identifier

  group

  :   Group identifier
