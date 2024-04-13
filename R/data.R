#' Example Simulated Repeated Measures scRNA-Seq Data
#'
#' Sample \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with
#' gene expression data simulated using the \code{RESCUE} R package.
#'
#' @format \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with
#' the following data in object slots:
#'
#' \describe{
#' \item{\code{counts}}{Matrix of raw counts with genes represented by rows and
#' cells represented by columns.}
#' \item{\code{normcounts}}{Data normalized using the variance stabilizing
#' transformation from the \code{sctransform} package.}
#'     \item{\code{colData}}{
#'         \describe{
#'             \item{sampleID}{Sample identifier}
#'             \item{subjectID}{Subject identifier}
#'             \item{time}{Timepoint identifier}
#'             \item{group}{Group identifier}
#'         }
#'     }
#' }
#'
"RecAM_sim_sce"

#' Example Simulated Repeated Measures Pseudo-bulk scRNA-Seq Data
#'
#' Sample \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#' containing pseudo-bulk data aggregated from scRNA-seq data simulated using
#' the \code{RESCUE} package.
#'
#'
#'
#'
#' @format \code{\link[SummarizedExperiment]{SummarizedExperiment}} object with
#' the following data in object slots:
#'
#' \describe{
#' \item{\code{counts}}{Matrix of raw pseudo-bulk counts with genes represented
#' by rows and samples represented by columns.}
#' \item{\code{normcounts}}{Pseudo-bulk normalized using the}
#' \code{\link[DESeq2]{VST} transformation from the \code{DESeq2} package
#'     \item{\code{colData}}{
#'         \describe{
#'             \item{sampleID}{Sample identifier}
#'             \item{subjectID}{Subject identifier}
#'             \item{time}{Timepoint identifier}
#'             \item{group}{Group identifier}
#'         }
#'     }
#' }
#'
"RecAM_sim_pb"
