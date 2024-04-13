


#' The RebelFit class
#'
#' @slot geneNames Vector containing gene names taken from the row names of the
#' \code{object} or \code{normalizedCounts} provided in the \code{rebelFit}
#' function call
#' @slot fits List of \code{\link[lmerTest]{lmerModLmerTest}} fit objects with one object for each gene.
#'  Not included if \code{outputFits=FALSE}.
#' @slot modelMatrix Design matrix for the model fixed effects.
#' @slot coefficients Dataframe of model coefficients with columns representing
#' the model fixed effects and each row representing an individual gene.
#' @slot originalFitVar Dataframe of residual and random effect variance values
#' from original model fits before empirical Bayes estimation. A logical value
#' indicating whether the original fit was singular or not is also included.
#' @slot miscFitInfo Miscellaneous information from model fits that is needed
#' for hypothesis testing.
#' @slot sampleVariable String denoting name of sample identifier variable in
#' the \code{colData} which was used as the sample-level random effect in the LMM
#' model. For pseudo-bulk data where no sample-level random effect is used,
#' this value is \code{NULL}.
#' @slot subjectVariable String denoting name of subject identifier variable in
#' the \code{colData} which was used as the subject-level random effect in the LMM
#' model.
#' @slot pseudoBulk Logical value indicating whether a pseudo-bulk or cell-level
#'  analysis is being performed.
#' @slot EBEstimates Dataframe with the empirical Bayes residual variance and
#' random intercept variance for each gene.
#' @slot EBParameters Dataframe with the parameters used to get empirical Bayes
#' estimates.
#' @slot vcovBetaList List of variance/covariance matrices with one matrix per
#' gene.
#' @slot vcovBetaAdjList List of adjusted variance/covariance matrices for each
#' gene calculated using the Kenward-Roger's method. Additional parameters used
#' for the Kenward-Roger's degrees of freedom estimation are also included. See
#' \code{pbkrtest} function \code{vcovAdj} for more information.
#'
#' @return
#' @export
#'
#' @examples
setClass("RebelFit",
         slots = c(
             geneNames="character",
             fits = "list",
             modelMatrix = "matrix",
             coefficients="matrix",
             originalFitVar = "data.frame",
             miscFitInfo = "list",
             sampleVariable = "character",
             subjectVariable = "character",
             pseudoBulk = "logical",
             EBEstimates = "data.frame",
             EBParameters= "data.frame",
             vcovBetaList = "list",
             vcovBetaAdjList = "list"
         )
)

## ADD in show method so not all data is output


#' Title
#'
#' @param RebelFitObj
#'
#' @return
#' @export
#'
#' @examples
getModelMatrix <- function(RebelFitObj) {
    methods::slot(RebelFitObj, "modelMatrix")
}

#' Title
#'
#' @param RebelFitObj
#'
#' @return
#' @export
#'
#' @examples
getEBEstimates <- function(RebelFitObj) {
    methods::slot(RebelFitObj, "EBEstimates")
}

#' Title
#'
#' @param RebelFitObj
#'
#' @return
#' @export
#'
#' @examples
getEBParameters <- function(RebelFitObj) {
    methods::slot(RebelFitObj, "EBParameters")
}

#' Title
#'
#' @param RebelFitObj
#'
#' @return
#' @export
#'
#' @examples
getOriginalFitVar<- function(RebelFitObj) {
    methods::slot(RebelFitObj, "originalFitVar")
}


#' Title
#'
#' @param RebelFitObj
#'
#' @return
#' @export
#'
#' @examples
getGeneNames<- function(RebelFitObj) {
    methods::slot(RebelFitObj, "geneNames")
}
