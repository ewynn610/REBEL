
#' Fit linear mixed models on cell level of pseudo-bulk data
#'
#' @param fixedEffects One sided linear formula describing the fixed-effects on the right of the ~ operator.
#'  Terms should be separated by + operators. Terms should be variables
#'  in \code{colData}.
#' @param normalizedCounts Optional matrix-like object containing normalized
#' counts. Used in conjunction with \code{colData} argument in place of object. Matrix
#' should contain genes in the rows and cells (cell-level data) or samples
#' (pseudo-bulk data) in the columns.
#' @param colData Optional dataframe object containing cell/sample level meta
#' data which is used with \code{normalizedCounts} argument in place of object.
#' Row names of dataframe should match column names of \code{normalizedCounts}.
#' @param pseudoBulk Logical value indicating whether a pseudo-bulk or
#' cell-level analysis is being performed.
#' @param subjectVariable String denoting name of subject identifier in
#' the \code{colData} to be used as the subject-level random effect in the
#' LMM model.
#' @param sampleVariable String denoting name of sample identifier variable in
#' the \code{colData} to be used as the sample-level random effect in the LMM model.
#' For pseudo-bulk data where no sample-level random effect is used, should be
#' \code{NULL}.
#' @param parallel Logical value indicating whether to use parallelization via
#' \code{mclapply}.
#' @param nCores Number of cores to use if \code{parallel} is \code{TRUE}.
#' @param outputFits Logical value indicating whether or not to include fit
#' objects from \code{\link[lmerTest]{lmerTest}} in the output. Only necessary
#' if user would like to inspect elements of the object. May use a large amount
#' of memory if TRUE.
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
    random_effects <- c()
    
    if (!is.null(sampleVariable)) {
      random_effects <- c(random_effects, paste0("(1|", sampleVariable, ")"))
    }
    
    if (!is.null(subjectVariable)) {
        random_effects <- c(random_effects, paste0("(1|", subjectVariable, ")"))
    }
    
    re_string <- paste(random_effects, collapse = " + ")
    
    formula <- update(fixedEffects, paste("expr~ . +", re_string))
    

    ## Fit models
    if(parallel == FALSE){
        ret <- pbapply::pblapply(X = 1:nrow(normalizedCounts),FUN = function(i){
            .fitGeneMod(normalizedCounts[i,], gene_name=gene_names[i], colData,
                        formula,subjectVariable = subjectVariable,
                        sampleVariable = sampleVariable,
                        pseudoBulk = pseudoBulk, outputFits = outputFits)
        })
    }else{
        ret = parallel::mclapply(1:nrow(normalizedCounts), mc.silent = TRUE,
                                 mc.cores = cores, function(i){
                                     .fitGeneMod(normalizedCounts[i,],
                                                 gene_name=gene_names[i],
                                                 colData,
                                                 formula,
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
                        formula,subjectVariable = subjectVariable,
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
    miscFitInfo=list(flist=flist, Ztlist=Ztlist, fr=fr, reTrms=reTrms, formula=formula)

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
                             pseudoBulk=pseudoBulk
                             )

    if(!is.null(sampleVariable)) methods::slot(RebelFitObj, "sampleVariable")=sampleVariable
    if(!is.null(subjectVariable)) methods::slot(RebelFitObj, "subjectVariable")=subjectVariable
    

    if(outputFits) methods::slot(RebelFitObj, "fits")=lapply(ret, function(x) x$fit)
    RebelFitObj
}





.fitGeneMod=function(geneExpr, gene_name, colData, formula,
                     subjectVariable, sampleVariable, pseudoBulk, outputFits){

    ## Bind gene expression and meta data
    dat_sub <- data.frame(cbind(colData, data.frame(expr = as.numeric(geneExpr))))

    ## Fit model
    fit <- tryCatch({
        tmp1 <- suppressMessages(lmerTest::lmer(formula = formula,
                                                data = dat_sub,
                                                REML = T))
    }, error = function(e) {
        ret_sub2 <- NULL
    })



    ## If model isn't null, return information
    if(!is.null(fit)){

        ## Is gene singular?
        singular=lme4::isSingular(fit)

        ## Variance info
        resVar=sigma(fit)^2
        
        # Start with basic data.frame
        modInfo <- data.frame(
          gene = gene_name,
          singular = singular,
          resVar = resVar
        )
        
        # Conditionally add subject random effect variance
        if (!is.null(subjectVariable)) {
          reVarSubj <- unlist(lme4::VarCorr(fit))[[subjectVariable]]
          modInfo$reVarSubj <- reVarSubj
        }
        
        # Conditionally add sample random effect variance
        if (!is.null(sampleVariable)) {
          reVarSamp <- unlist(lme4::VarCorr(fit))[[sampleVariable]]
          modInfo$reVarSamp <- reVarSamp
        }

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
