
#' Title
#'
#' @param fixedEffects
#' @param normalizedCounts
#' @param colData
#' @param pseudoBulk
#' @param subjectVariable
#' @param sampleVariable
#' @param REML
#' @param parallel
#' @param cores
#' @param outputFits
#'
#' @return
#' @export
#'
#' @examples
rebelLMM <- function(fixedEffects, # Formula for fixed effects
                     normalizedCounts = NULL, # Matrix of transformed RNA-Seq counts where rows are genes and columns are samples
                     colData = NULL, # A data frame with meta data
                     pseudoBulk, # Logical value, pseudo-bulk or not (if not cell level, pseudo-bulk)
                     subjectVariable,
                     sampleVariable=NULL,
                     REML = TRUE, # Fit mixed models using REML or ML
                     parallel = FALSE,
                     cores = 2,
                     outputFits=FALSE

){  ## Insufficient Information


    if(is.null(normalizedCounts)==TRUE ) {
        stop("An expression matrix must be provided.")}

    if(is.null(fixedEffects)==TRUE ) {
        stop("A fixed effect formula must be provided.")}

    if(is.null(colData)==TRUE ) {
        stop("colData is missing.")}

    ## Inconsistent information

    if((ncol(normalizedCounts)==nrow(colData))==FALSE ) {
        stop("The expression matrix and sample data include differing numbers of samples.")}

    # Gene Names
    if(is.null(rownames(normalizedCounts))==TRUE){rownames(normalizedCounts)<-paste0("gene",seq(1,nrow(normalizedCounts)))}
    gene_names = as.character(rownames(normalizedCounts))


    # Make sure normalizedCounts is a matrix
    normalizedCounts <- as.matrix(normalizedCounts)

    ## Add random intercepts to formula
    if(pseudoBulk){
        formula=update(fixedEffects, paste0("expr~.+(1|",
                                            subjectVariable, ")"))

    }else{
        formula=update(fixedEffects, paste0("expr~.+(1|", sampleVariable,
                                            ")+(1|", subjectVariable, ")"))

    }

    ## Fit models
    if(parallel == FALSE){
        ret <- pbapply::pblapply(X = 1:nrow(normalizedCounts),FUN = function(i){
            .fitGeneMod(normalizedCounts[i,], gene_name=gene_names[i], colData,
                        formula, REML,subjectVariable = subjectVariable,
                        sampleVariable = sampleVariable,
                        pseudoBulk = pseudoBulk, outputFits = outputFits)
        })
    }else{
        ret = parallel::mclapply(1:nrow(normalizedCounts), mc.silent = TRUE,
                                 mc.cores = cores, function(i){
                                     .fitGeneMod(normalizedCounts[i,],
                                                 gene_name=gene_names[i],
                                                 colData,
                                                 formula, REML,
                                                 subjectVariable = subjectVariable,
                                                 sampleVariable = sampleVariable,
                                                 pseudoBulk = pseudoBulk,
                                                 outputFits = outputFits)
                                 })
    }

    ## Check if all models are NULL
    idx=which(!vapply(ret, function(x) is.null(x$modInfo), FUN.VALUE = logical(1)))[1]
    if(is.na(idx)){
        stop("All model fits are NULL. Check for error in function call.")
    }

    ## Get one fit object to extract model info from
    if(outputFits){
        fit=ret[[idx]]$fit
    }else{
        fit=.fitGeneMod(normalizedCounts[idx,], gene_name=gene_names[idx], colData,
                        formula, REML,subjectVariable = subjectVariable,
                        sampleVariable = sampleVariable,
                        pseudoBulk = pseudoBulk, outputFits = TRUE)
        fit=fit$fit
    }

    ## Get fit information
    flist=lme4::getME(fit, "flist")
    Ztlist=lme4::getME(fit, "Ztlist")
    modelMatrix=lme4::getME(fit, "X")
    fr=model.frame(fit)
    reTrms <- lme4:::mkReTrms(lme4::findbars(formula), fr)
    
    ## Collect misc info
    miscFitInfo=list(flist=flist, Ztlist=Ztlist, fr=fr, reTrms=reTrms)

    ## Gather coefficients
    coefficients=as.matrix(dplyr::bind_rows(
        lapply(ret, function(x) x$coefficients)))
    rownames(coefficients)=gene_names

    ## Gather model info
    modInfo=data.frame(dplyr::bind_rows(
        lapply(ret, function(x) x$modInfo)))
    rownames(modInfo)=gene_names

    RebelFitObj=methods::new("RebelFit",
                             geneNames=gene_names,
                             coefficients=coefficients,
                             modelMatrix=modelMatrix,
                             originalFitVar=modInfo,
                             miscFitInfo=miscFitInfo,
                             subjectVariable=subjectVariable,
                             pseudoBulk=pseudoBulk
                             )

    if(!pseudoBulk) methods::slot(RebelFitObj, "sampleVariable")=sampleVariable

    if(outputFits) methods::slot(RebelFitObj, "fits")=lapply(ret, function(x) x$fit)
    RebelFitObj
}





.fitGeneMod=function(geneExpr, gene_name, colData, formula, REML,
                     subjectVariable, sampleVariable, pseudoBulk, outputFits){

    ## Bind gene expression and meta data
    dat_sub <- data.frame(cbind(colData, data.frame(expr = as.numeric(geneExpr))))

    ## Fit model
    fit <- tryCatch({
        tmp1 <- suppressMessages(lmerTest::lmer(formula = formula,
                                                data = dat_sub,
                                                REML = REML))
    }, error = function(e) {
        ret_sub2 <- NULL
    })



    ## If model isn't null, return information
    if(!is.null(fit)){

        ## Is gene singular?
        singular=lme4::isSingular(fit)

        ## Variance info
        resVar=sigma(fit)^2
        reVarSubj=unlist(lme4::VarCorr(fit))[[subjectVariable]]
        modInfo=data.frame(gene=gene_name,
                            singular=singular,
                            resVar=resVar,
                            reVarSubj=reVarSubj)

        ## If cell level, add sample RE variance
        if(!pseudoBulk) modInfo$reVarSamp=reVar=unlist(lme4::VarCorr(fit))[[sampleVariable]]

        ## Misc info for adj variance calculation
        devFun=lme4::getME(fit, "devfun")

        ## Output info
        ret_sub=list(gene=gene_name,
                     modInfo=modInfo,
                     devFun=devFun,
                     coefficients=lme4::fixef(fit))
        if(outputFits) ret_sub$fit=fit



        ret_sub
    }else{
        ret_sub=list(gene=gene_name,
                     modInfo=NULL,
                     devfun=NULL,
                     coefficients=NULL)
        if(outputFits) ret_sub$fit=NULL
        ret_sub
    }

}
