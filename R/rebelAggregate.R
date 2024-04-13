
#' Aggregate data for pseudo-bulk analysis
#'
#' @param object \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' containing count values in one of the object assays and cell level meta data
#' in the colData slot.
#' @param countsAssay String specifying assay in object which holds count values.
#' @param counts Optional matrix-like object containing counts. Used in
#' conjunction with \code{colData} argument in place of object. Matrix should
#' contain genes in the rows and cells in the columns.
#' @param colData Optional dataframe object containing cell level meta data
#' which is used with counts argument in place of object. Row names of dataframe
#'  should match column names of \code{normalizedCounts}.
#' @param sampleVariable String denoting name of sample identifier variable in
#' the \code{colData} which counts should be aggregated over.
#'
#' @return \code{\link[SummarizedExperiment]{SummarizedExperiment}} and includes
#' pseudo-bulk counts in the \code{counts} assay. All colData variables that do not
#' vary within the \code{sampleVariable} are also included in the \code{colData}
#' slot.
#' @export
#'
#' @author Elizabeth Wynn
#'
#' @examples
#' ## Read in single cell data
#' data("RecAM_sim_sce")
#'
#' ## Aggregate into pseudo-bulk data
#' pb_dat <- rebelAggregate(object=RecAM_sim_sce, countsAssay="counts", sampleVariable="sampleID")

rebelAggregate=function(object,countsAssay="counts", counts=NULL, colData=NULL, sampleVariable){
    if(exists("object")){
        if(inherits(object,"SingleCellExperiment")){
            counts=assay(object, countsAssay)
            colData=SingleCellExperiment::colData(object)
        }else{
            stop("Object should be a SingleCellExperiment object. Alternatively, you can provide counts and colData explicitly")
        }
    }


    ##Get summarized counts
    pbCounts=t(as.matrix(Matrix.utils::aggregate.Matrix(t(counts),
                                                      colData[,sampleVariable],
                                                      fun="sum")))

    ## Put back in original order
    pbCounts=pbCounts[,unique(colData[,sampleVariable])]


    ## Which varriables don't change within a sample?
    vars_keep=vapply(colnames(colData)[colnames(colData)!=sampleVariable], function(x){
        form=as.formula(paste0(x, "~", sampleVariable))
        df=aggregate(form, data = colData,
                     FUN = function(y) length(unique(y))==1)
        all(df[,x])

    }, FUN.VALUE = logical(1))

    ## Get rid of variables that vary within a sample
    colData=colData[!duplicated(colData[,sampleVariable]),
                    c(sampleVariable,names(vars_keep)[vars_keep])]
    rownames(colData)=NULL

    ## Make summarized experiment
    SummarizedExperiment::SummarizedExperiment(assays=list(counts=pbCounts), colData = colData)
}

