# The RebelFit class

The RebelFit class

## Slots

- `geneNames`:

  Vector containing gene names taken from the row names of the `object`
  or `normalizedCounts` provided in the `rebelFit` function call

- `fits`:

  List of
  [`lmerModLmerTest`](https://rdrr.io/pkg/lmerTest/man/lmerModLmerTest-class.html)
  fit objects with one object for each gene. Not included if
  `outputFits=FALSE`.

- `modelMatrix`:

  Design matrix for the model fixed effects.

- `coefficients`:

  Dataframe of model coefficients with columns representing the model
  fixed effects and each row representing an individual gene.

- `originalFitVar`:

  Dataframe of residual and random effect variance values from original
  model fits before empirical Bayes estimation. A logical value
  indicating whether the original fit was singular or not is also
  included.

- `miscFitInfo`:

  Miscellaneous information from model fits that is needed for
  hypothesis testing.

- `sampleVariable`:

  String denoting name of sample identifier variable in the `colData`
  which was used as the sample-level random effect in the LMM model. For
  pseudo-bulk data where no sample-level random effect is used, this
  value is `NULL`.

- `subjectVariable`:

  String denoting name of subject identifier variable in the `colData`
  which was used as the subject-level random effect in the LMM model.

- `pseudoBulk`:

  Logical value indicating whether a pseudo-bulk or cell-level analysis
  is being performed.

- `EBEstimates`:

  Dataframe with the empirical Bayes residual variance and random
  intercept variance for each gene.

- `EBParameters`:

  Dataframe with the parameters used to get empirical Bayes estimates.

- `vcovBetaList`:

  List of variance/covariance matrices with one matrix per gene.

- `vcovBetaAdjList`:

  List of adjusted variance/covariance matrices for each gene calculated
  using the Kenward-Roger's method. Additional parameters used for the
  Kenward-Roger's degrees of freedom estimation are also included. See
  `pbkrtest` function `vcovAdj` for more information.
