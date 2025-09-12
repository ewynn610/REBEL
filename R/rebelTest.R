
#' Run hypothesis testing using models fit using REBEL
#'
#' @param RebelFitObj A RebelFit object obtained using the \code{rebelFit}
#' function.
#' @param coef 	Character string indicating which coefficient to summarize. If
#' both coef and contrast are included, coef takes precedence.
#' @param contrast A numeric vector or matrix specifying the linear contrast(s)
#' of fixed effects to test. If a vector, it must be the same length of the number of
#' fixed effects and represents a one-dimensional contrast. If a matrix, each row represents
#' a contrast, and the number of columns must match the number of fixed effects. 
#' Multi-row matrices are used for joint F-tests.
#'
#' @return a dataframe with testing information, including Kenward-Rogers
#' estimates degrees of freedom as well as the
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
#'                      pseudoBulk = FALSE)
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
#' ##  multi-dimensional contrast (F-test)
#' contrast_mat <- rbind(c(0, 1, 0, 0),
#'                    c(0, 0, 1, 0),
#'                    c(0, 0, 0, 1))
#' ## simultaneous test of all coefficients
#' joint_contrast_test=rebelTest(cell_fit, contrast = contrast_mat)


rebelTest=function(RebelFitObj, coef=NULL, contrast=NULL){
    beta_mat=methods::slot(RebelFitObj, "coefficients")
    vcovBetaList=methods::slot(RebelFitObj, "vcovBetaList")
    vcovBetaAdjList=methods::slot(RebelFitObj, "vcovBetaAdjList")
    EBEstimates=methods::slot(RebelFitObj, "EBEstimates")
    if(!is.null(coef)){
        if(!(coef %in% colnames(beta_mat))) stop("Coef must be one of the names of the model fixed effects")
        contrast=.contrastFromCoef(coef, beta_mat)
    }
    joint_flag <- nrow(rbind(contrast)) > 1

    geneNames=getGeneNames(RebelFitObj)

    sumTab=if(joint_flag){
      beta_h=rep(0, ncol(beta_mat))
      
      dplyr::bind_rows(lapply(geneNames, function(x){

        test_res<-.KR_adjust(PhiA=vcovBetaAdjList[[x]], Phi =vcovBetaList[[x]], L = contrast,
                   beta = beta_mat[x,], betaH = beta_h )
        Fvalue = test_res$FstatU
        ndf = test_res$ndf
        ddf = test_res$ddf
        sigma = sqrt(EBEstimates[x, "resVar"])
        Fscale = test_res$F.scaling
       
        
        MS <- Fvalue * sigma^2
        Fvalue <- Fvalue * Fscale
        pvalue <- pf(q = Fvalue, df1 = ndf, df2 = ddf, lower.tail = FALSE)
        data.frame(Sum.sq = MS * ndf, Mean.sq = MS, Num.df = ndf, 
                   Den.df = ddf, F.value = Fvalue, p_val_raw = pvalue)
      }))
    }else if(!joint_flag){
      dplyr::bind_rows(lapply(geneNames, function(x){
        Estimate=drop(contrast%*% beta_mat[x,])
        Std.Error=sqrt(Matrix::drop(t(contrast)%*%vcovBetaAdjList[[x]]%*%contrast))
        t.value=Estimate/Std.Error
        df=pbkrtest::Lb_ddf(L = contrast, V0 = vcovBetaList[[x]], Vadj = vcovBetaAdjList[[x]])
        p_val_raw=2*pt(-abs(t.value), df)
        data.frame(Estimate=Estimate, Std.Error=Std.Error,t.value=t.value,
                   df=df, p_val_raw=p_val_raw
        )
      }))
    }
      
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

## Internal helper function from pbkrtest package
## Ulrich Halekoh, Søren Højsgaard (2014). A Kenward-Roger Approximation and Parametric Bootstrap
## Methods for Tests in Linear Mixed Models - The R Package pbkrtest. Journal of Statistical
## Software, 59(9), 1-30. URL https://www.jstatsoft.org/v59/i09/.

.KR_adjust<-function (PhiA, Phi, L, beta, betaH) 
{
  Theta <- t(L) %*% solve(L %*% Phi %*% t(L), L)
  P <- attr(PhiA, "P")
  W <- attr(PhiA, "W")
  A1 <- A2 <- 0
  ThetaPhi <- Theta %*% Phi
  n.ggamma <- length(P)
  for (ii in 1:n.ggamma) {
    for (jj in c(ii:n.ggamma)) {
      e <- ifelse(ii == jj, 1, 2)
      ui <- ThetaPhi %*% P[[ii]] %*% Phi
      uj <- ThetaPhi %*% P[[jj]] %*% Phi
      A1 <- A1 + e * W[ii, jj] * (.spur(ui) * .spur(uj))
      A2 <- A2 + e * W[ii, jj] * sum(ui * t(uj))
    }
  }
  q <- as.numeric(Matrix::rankMatrix(L))
  B <- (1/(2 * q)) * (A1 + 6 * A2)
  g <- ((q + 1) * A1 - (q + 4) * A2)/((q + 2) * A2)
  c1 <- g/(3 * q + 2 * (1 - g))
  c2 <- (q - g)/(3 * q + 2 * (1 - g))
  c3 <- (q + 2 - g)/(3 * q + 2 * (1 - g))
  EE <- 1 + (A2/q)
  VV <- (2/q) * (1 + B)
  EEstar <- 1/(1 - A2/q)
  VVstar <- (2/q) * ((1 + c1 * B)/((1 - c2 * B)^2 * (1 - c3 * 
                                                       B)))
  V0 <- 1 + c1 * B
  V1 <- 1 - c2 * B
  V2 <- 1 - c3 * B
  V0 <- ifelse(abs(V0) < 1e-10, 0, V0)
  rho <- 1/q * (.divZero(1 - A2/q, V1))^2 * V0/V2
  df2 <- 4 + (q + 2)/(q * rho - 1)
  F.scaling <- ifelse(abs(df2 - 2) < 0.01, 1, df2 * (1 - A2/q)/(df2 - 
                                                                  2))
  aux <- c(A1 = A1, A2 = A2, V0 = V0, V1 = V1, V2 = V2, rho = rho, 
           F.scaling = F.scaling)
  betaDiff <- cbind(beta - betaH)
  Wald <- as.numeric(t(betaDiff) %*% t(L) %*% solve(L %*% PhiA %*% 
                                                      t(L), L %*% betaDiff))
  WaldU <- as.numeric(t(betaDiff) %*% t(L) %*% solve(L %*% 
                                                       Phi %*% t(L), L %*% betaDiff))
  FstatU <- Wald/q
  pvalU <- pf(FstatU, df1 = q, df2 = df2, lower.tail = FALSE)
  Fstat <- F.scaling * FstatU
  pval <- pf(Fstat, df1 = q, df2 = df2, lower.tail = FALSE)
  stats <- list(ndf = q, ddf = df2, Fstat = Fstat, p.value = pval, 
                F.scaling = F.scaling, FstatU = FstatU, p.valueU = pvalU, 
                aux = aux)
  stats
}

.spur<-function (U) 
{
  sum(Matrix::diag(U))
}

.divZero<-function (x, y, tol = 1e-14) 
{
  if (abs(x) < tol & abs(y) < tol) 
    1
  else x/y
}
