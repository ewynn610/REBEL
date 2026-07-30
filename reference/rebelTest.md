# Run hypothesis testing using models fit using REBEL

Run hypothesis testing using models fit using REBEL

## Usage

``` r
rebelTest(RebelFitObj, coef = NULL, contrast = NULL)
```

## Arguments

- RebelFitObj:

  A RebelFit object obtained using the `rebelFit` function.

- coef:

  Character string indicating which coefficient to summarize. If both
  coef and contrast are included, coef takes precedence.

- contrast:

  A numeric vector or matrix specifying the linear contrast(s) of fixed
  effects to test. If a vector, it must be the same length of the number
  of fixed effects and represents a one-dimensional contrast. If a
  matrix, each row represents a contrast, and the number of columns must
  match the number of fixed effects. Multi-row matrices are used for
  joint F-tests.

## Value

a dataframe with testing information, including Kenward-Rogers estimates
degrees of freedom as well as the Benjamini-Hochberg adjusted p-values,
for each gene. Genes appear in the order that they appear in the
dataset.

## Examples

``` r

## Read in cell-level data (scTransform normalized data saved in normcounts slot)
data("RecAM_sim_sce")

## Just use first 10 genes
RecAM_sim_sce_fil <- RecAM_sim_sce[1:10,]

## Fit models with time, group and time/group interaction effect
cell_fit <- rebelFit(object=RecAM_sim_sce_fil,
                     fixedEffects = ~time*group,
                     subjectVariable ="subjectID", sampleVariable = "sampleID",
                     pseudoBulk = FALSE)
#> [1] "Fitting LMM Models"
#> [1] "Getting Empirical Bayes Estimates"
#> [1] "Getting Parameters for Kenward-Rogers method"

## Run test on interaction coefficient
interaction_coef_test=rebelTest(cell_fit, coef="timetime1:groupgroup1")

## Run test on interaction coefficient using a contrast
interaction_contrast_test=rebelTest(cell_fit, contrast = c(0,0,0,1))

## Using coefficient and contrast to test interaction produce identical results
identical(interaction_coef_test, interaction_contrast_test)
#> [1] TRUE

##  multi-dimensional contrast (F-test)
contrast_mat <- rbind(c(0, 1, 0, 0),
                   c(0, 0, 1, 0),
                   c(0, 0, 0, 1))
## simultaneous test of all coefficients
joint_contrast_test=rebelTest(cell_fit, contrast = contrast_mat)
```
