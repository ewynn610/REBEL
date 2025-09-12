
#' Title
#'
#' @param RebelFitObj
#'
#' @return
#' @export
#'
#' @examples
rebelEmpBayes=function(RebelFitObj){

    geneNames=getGeneNames(RebelFitObj)
    ## Get param/observation info
    modelMatrix=getModelMatrix(RebelFitObj)
    numParams=ncol(modelMatrix)
    numObs=nrow(modelMatrix)

    ## Get RE variables
    flist=methods::slot(RebelFitObj, "miscFitInfo")[["flist"]]
    subjectVariable=methods::slot(RebelFitObj, "subjectVariable")
    if (length(subjectVariable) == 0) {
      subjectVariable <- NULL
    } else {
      subjectVariable <- subjectVariable
    }
    sampleVariable=methods::slot(RebelFitObj, "sampleVariable")
    if (length(sampleVariable) == 0) {
      sampleVariable <- NULL
    } else {
      sampleVariable <- sampleVariable
    }
    pseudoBulk=methods::slot(RebelFitObj, "pseudoBulk")
    modInfo=getOriginalFitVar(RebelFitObj)


    ## Number of between/within variables
    paramTypeN=.findNumParamTypes(modelMatrix=data.frame(modelMatrix),flist=flist,
                                  subjectVariable=subjectVariable,
                                  sampleVariable=sampleVariable,
                                  pseudoBulk=pseudoBulk)


    ## Find empirical bayes DF
    
    ## Initialize btw DF
    btw_df=0
    if(!is.null(subjectVariable)){
      btw_df=dfSubjRE=length(unique(flist[[subjectVariable]]))-paramTypeN[["btw_subj"]]
    } else dfSubjRE=NULL
    if(!pseudoBulk){
        dfRes=dfSampRE=length(unique(flist[[sampleVariable]]))-paramTypeN[["btw_samp"]]
    }else{
      if(!is.null(sampleVariable)){
        dfSampRE=length(unique(flist[[sampleVariable]]))-paramTypeN[["btw_samp"]]
        btw_df=dfSampRE+btw_df
      }else dfSampRE=NULL
      
      dfRes=numObs-numParams-btw_df
    }
    priorDF=c(reVarSubj=dfSubjRE, reVarSamp=dfSampRE, resVar=dfRes)

    ## Get empirical bayes estimates
    empBayesInfo=lapply(names(priorDF), function(varType){
        vars=if(varType=="resVar")modInfo[,varType] else modInfo[,varType][!modInfo$singular]
        FFit=limma::fitFDist(vars, priorDF[varType])
        ests=data.frame(.squeezeVar(modInfo[,varType],priorDF[varType], FFit$scale,
                                    FFit$df2))
        rownames(ests)=geneNames
        colnames(ests)=varType
        params=data.frame(s0=FFit$scale, df0=FFit$df2, dfg=priorDF[varType])
        list(ests=ests, params=params)
    })
    empBayesEsts=data.frame(dplyr::bind_cols(lapply(empBayesInfo, function(x) x$ests)))
    empBayesParams=data.frame(dplyr::bind_rows(lapply(empBayesInfo, function(x) x$params)))
    rownames(empBayesParams)=names(priorDF)

    ## Save estimates/parameters in object
    methods::slot(RebelFitObj, "EBEstimates")=empBayesEsts
    methods::slot(RebelFitObj, "EBParameters")=empBayesParams

    RebelFitObj



}

## Find number of betweein/within parameters
.findNumParamTypes=function(modelMatrix, flist, subjectVariable, sampleVariable,
                            pseudoBulk){
    fixEffNames=colnames(modelMatrix)
    btw_subj_names <- fixEffNames  # default if no subjectVariable
    
    out <- c()
    
    if (!is.null(subjectVariable)) {
      btw_subj_log <- .findParamTypeLog(modelMatrix,flist,subjectVariable, fixEffNames)
      btw_subj_names <- names(btw_subj_log)[btw_subj_log]
      btw_subj_n <- sum(btw_subj_log)
      out["btw_subj"] <- btw_subj_n
      
    }
    
    # Add sample-level info if sampleVariable exists
    if (!is.null(sampleVariable)) {
      btw_samp_log <- .findParamTypeLog(modelMatrix, flist, sampleVariable, btw_subj_names)
      btw_samp_n <- sum(btw_samp_log)
      out["btw_samp"] <- btw_samp_n
    }
    return(out)
    
}

.findParamTypeLog=function(modelMatrix,flist, var, fixEffNames){
    modelMatrix[,var]=flist[[var]]

    btw_subj_log=vapply(fixEffNames, function(x){
        form=as.formula(paste0(x, "~", var))

        df=aggregate(form, data = modelMatrix,
                     FUN = function(y) length(unique(y))==1)
        all(df[,x])

    }, FUN.VALUE = logical(1))
}


## *NOTE: Internal Function taken from limma package
## Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma
##powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids
##Research 43(7), e47

## Have to use internal function rather than exported limma function because
## I want to estimate parameters from non-singular models and then
## estimate variance with full set of data
.squeezeVar=function (var, df, var.prior, df.prior)
{
    n <- length(var)
    isfin <- is.finite(df.prior)
    if (all(isfin))
        return((df * var + df.prior * var.prior)/(df + df.prior))
    if (length(var.prior) == n) {
        var.post <- var.prior
    }
    else {
        var.post <- rep_len(var.prior, length.out = n)
    }
    if (any(isfin)) {
        i <- which(isfin)
        if (length(df) > 1)
            df <- df[i]
        df.prior <- df.prior[i]
        var.post[i] <- (df * var[i] + df.prior * var.post[i])/(df +
                                                                   df.prior)
    }
    var.post
}















