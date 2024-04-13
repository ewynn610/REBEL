
#' Run hypothesis testing using models fit using REBEL
#'
#' @param RebelFitObj A RebelFit object obtained using the \code{rebelFit}
#' function.
#' @param coef 	Character string indicating which coefficient to summarize. If
#' both coef and contrast are included, coef takes precedence.
#' @param contrast A numeric vector the same length as the number of fixed
#' effects representing a one-dimensional contrast.
#'
#' @return a dataframe with testing information, including Kenward-Rogers
#' estimates for the standard error and degrees of freedom as well as the
#' Benjamini-Hochberg adjusted p-values, for each gene. Genes appear in the
#' order that they appear in the dataset.
#' @export
#'
#' @examples
#'
#' ## Read in cell-level data (scTransform normalized data saved in normcounts slot)
#' data("RecAM_sim_sce")
#'
#' ## Just use first 10 genes
#' RecAM_sim_sce_fil <- RecAM_sim_sce[1:10,]
#'
#' ## Fit models with time, group and time/group interaction effect
#' cell_fit <- rebelFit(object=RecAM_sim_sce_fil,
#'                      fixedEffects = ~time*group,
#'                      subjectVariable ="subjectID", sampleVariable = "sampleID",
#'                      pseudoBulk = F)
#'
#' ## Run test on interaction coefficient
#' interaction_coef_test=rebelTest(cell_fit, coef="timetime1:groupgroup1")
#'
#' ## Run test on interaction coefficient using a contrast
#' interaction_contrast_test=rebelTest(cell_fit, contrast = c(0,0,0,1))
#'
#' ## Using coefficient and contrast to test interaction produce identical results
#' identical(interaction_coef_test, interaction_contrast_test)
#'


rebelTest=function(RebelFitObj, coef=NULL, contrast=NULL){
    beta_mat=methods::slot(RebelFitObj, "coefficients")
    vcovBetaList=methods::slot(RebelFitObj, "vcovBetaList")
    vcovBetaAdjList=methods::slot(RebelFitObj, "vcovBetaAdjList")


    if(!is.null(coef)){
        if(!(coef %in% colnames(beta_mat))) stop("Coef must be one of the names of the model fixed effects")
        contrast=.contrastFromCoef(coef, beta_mat)
    }

    geneNames=getGeneNames(RebelFitObj)

    sumTab=dplyr::bind_rows(lapply(geneNames, function(x){
        Estimate=drop(contrast%*% beta_mat[x,])
        Std.Error=sqrt(Matrix::drop(t(contrast)%*%vcovBetaAdjList[[x]]%*%contrast))
        t.value=Estimate/Std.Error
        df=pbkrtest::Lb_ddf(L = contrast, V0 = vcovBetaList[[x]], Vadj = vcovBetaAdjList[[x]])
        p_val_raw=2*pt(-abs(t.value), df)
        data.frame(Estimate=Estimate, Std.Error=Std.Error,t.value=t.value,
                   df=df, p_val_raw=p_val_raw
                   )
    }))
    rownames(sumTab)=geneNames
    sumTab$p_val_adj=p.adjust(sumTab$p_val_raw, "BH")
    sumTab
}

.contrastFromCoef=function(coef, beta_mat){
    contrast=rep(0, length = ncol(beta_mat))
    idx_coef=which(colnames(beta_mat)==coef)
    contrast[idx_coef]=1
    contrast

}
