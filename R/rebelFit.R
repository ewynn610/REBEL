#' Fit REBEL Models
#'
#' @param object \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' (cell-level data) or \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object (pseudo-bulk data) containing normalized count values in one of the
#' object assays and cell/sample level meta data in the \code{colData} slot.
#' @param assay String specifying assay in object which holds normalized count
#' values.
#' @param fixedEffects One sided linear formula describing the fixed-effects on the right of the ~ operator.
#'  Terms should be separated by + operators. Terms should be variables
#'  in \code{colData}.
#' @param subjectVariable String denoting name of subject identifier in
#' the \code{colData} to be used as the subject-level random effect in the
#' LMM model.
#' @param sampleVariable String denoting name of sample identifier variable in
#' the \code{colData} to be used as the sample-level random effect in the LMM model.
#' For pseudo-bulk data where no sample-level random effect is used, should be
#' \code{NULL}.
#' @param normalizedCounts Optional matrix-like object containing normalized
#' counts. Used in conjunction with \code{colData} argument in place of object. Matrix
#' should contain genes in the rows and cells (cell-level data) or samples
#' (pseudo-bulk data) in the columns.
#' @param colData Optional dataframe object containing cell/sample level meta
#' data which is used with \code{normalizedCounts} argument in place of object.
#' Row names of dataframe should match column names of \code{normalizedCounts}.
#' @param pseudoBulk Logical value indicating whether a pseudo-bulk or
#' cell-level analysis is being performed.
#' @param parallel Logical value indicating whether to use parallelization via
#' \code{mclapply}.
#' @param nCores Number of cores to use if \code{parallel} is \code{TRUE}.
#' @param outputFits Logical value indicating whether or not to include fit
#' objects from \code{\link[lmerTest]{lmerTest}} in the output. Only necessary
#' if user would like to inspect elements of the object. May use a large amount
#' of memory if TRUE.
#' @param quiet Logical value indicating whether messages should be printed at
#' each step.
#' @param REML Logical value indicating if LMM models should be fit using REML
#' or regular ML.
#'
#' @return A \code{\link[REBEL]{RebelFit-class} object}.
#' @export
#'
#' @examples
#' ## Run cell-level analysis
#' ## Read in data (scTransform normalized data saved in normcounts slot)
#' data("RecAM_sim_sce")
#'
#' ## Just use first 10 genes
#' RecAM_sim_sce_fil <- RecAM_sim_sce[1:10,]
#'
#' ## Fit models with time, group and time/group interaction effect
#' cell_fit <-rebelFit(object=RecAM_sim_sce_fil,
#'                                    fixedEffects = ~time*group,
#'                                    subjectVariable ="subjectID",
#'                                    sampleVariable = "sampleID",
#'                                    pseudoBulk = F)
#'
#' ## Run pseudo-bulk analysis
#' ## Read in data (DESeq2 VST normalized data saved in normcounts slot)
#' data("RecAM_sim_pb")
#'
#' ## Just use first 10 genes
#' RecAM_sim_pb_fil <- RecAM_sim_pb[1:10,]
#'
#' ## Fit models with time, group and time/group interaction effect
#' pb_fit <- rebelFit(object=RecAM_sim_pb_fil,
#'                                  fixedEffects = ~time*group,
#'                                  subjectVariable ="subjectID",
#'                                  pseudoBulk = T)
rebelFit <- function(object,assay="normcounts", fixedEffects,subjectVariable,
                     sampleVariable=NULL, normalizedCounts=NULL,
                     colData=NULL, pseudoBulk=T,parallel=F,
                     nCores=1, outputFits=F, quiet=F, REML=T){

    if(exists("object")){

        if(inherits(object,"SingleCellExperiment")|inherits(object,"SummarizedExperiment")){
            normalizedCounts=assay(object, assay)
            colData=colData(object)
        }else{
            stop("Object should be a SingleCellExperiment or SummarizedExperiment object.
                 Alternatively, you can provide counts and colData explicitly")
        }
    }
    if(!quiet) print("Fitting LMM Models")
    ## Fit LMM models
    RebelFitObj=rebelLMM(fixedEffects=fixedEffects,
                  normalizedCounts=normalizedCounts,
                  colData=colData,
                  pseudoBulk=pseudoBulk,
                  subjectVariable=subjectVariable,
                  sampleVariable = sampleVariable,
                  REML = REML,
                  parallel = parallel,
                  cores = nCores,
                  outputFits= outputFits)

    if(!quiet) print("Getting Empirical Bayes Estimates")
    ## Fit empirical Bayes
    RebelFitObj=rebelEmpBayes(RebelFitObj)

    if(!quiet) print("Getting Parameters for Kenward-Rogers method")
    ## Get adjusted Vcov for KR
    RebelFitObj=rebelKRParams(RebelFitObj, parallel = parallel, nCores = nCores)


}














