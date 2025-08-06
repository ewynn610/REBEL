

#' Title
#'
#' @param RebelFitObj
#' @param parallel
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
rebelKRParams=function (RebelFitObj, parallel=FALSE, nCores=1) {

    EBEstimates=getEBEstimates(RebelFitObj)
    pseudoBulk=methods::slot(RebelFitObj, "pseudoBulk")
    flist=methods::slot(RebelFitObj, "miscFitInfo")[["flist"]]
    Ztlist=methods::slot(RebelFitObj, "miscFitInfo")[["Ztlist"]]
    fr=methods::slot(RebelFitObj, "miscFitInfo")[["fr"]]
    reTrms=methods::slot(RebelFitObj, "miscFitInfo")[["reTrms"]]
    subjectVariable=methods::slot(RebelFitObj, "subjectVariable")
    sampleVariable=methods::slot(RebelFitObj, "sampleVariable")
    modelMatrix=getModelMatrix(RebelFitObj)
    nObs=nrow(modelMatrix)
    gene_names=getGeneNames(RebelFitObj)




    devFuns=methods::slot(RebelFitObj, "miscFitInfo")[["devFuns"]]

    theta_list=lapply(gene_names,function(y){
        vapply(c(EBEstimates[y,"reVarSamp"], EBEstimates[y,"reVarSubj"]),
               function(x) .varCov2Theta(x, EBEstimates[y,"resVar"]),
               FUN.VALUE = numeric(1))
    })
    sigma_vals=sqrt(EBEstimates[, "resVar"])
    names(theta_list)=names(sigma_vals)=gene_names


    if (parallel) {
      vcovBetaAdj = parallel::mclapply(gene_names, mc.silent = TRUE, 
                                       mc.cores = nCores, function(gene) {
                                         ## Have to define devFun inside mclapply or won't work
                                         devFun=lme4::mkLmerDevfun(fr, X=modelMatrix, reTrms = reTrms)
                                         test=.VCovAdj_1Gene(theta=theta_list[[gene]], sigma=sigma_vals[gene], devFun=devFun,
                                                             Ztlist, sampleVariable=sampleVariable,
                                                             subjectVariable=subjectVariable,
                                                             nObs=nObs, modelMatrix=modelMatrix, pseudoBulk=pseudoBulk,
                                                             flist=flist)
                                       })
    }
    else {
      ## devFun can be outside of apply unless using parallel
      ## Can use same devfun to get all inverse cholesky factors - only theta changes
      devFun=lme4::mkLmerDevfun(fr, X=modelMatrix, reTrms = reTrms)
      vcovBetaAdj = pbapply::pblapply(gene_names, function(gene) {
        .VCovAdj_1Gene(theta = theta_list[[gene]], sigma = sigma_vals[gene], 
                       devFun = devFun, Ztlist, sampleVariable = sampleVariable, 
                       subjectVariable = subjectVariable, nObs = nObs, 
                       modelMatrix = modelMatrix, pseudoBulk = pseudoBulk, 
                       flist = flist)
      })
      ## Remake reterms because they get changed when devfun theta changes.
      ## Just for consistency
      methods::slot(RebelFitObj, "miscFitInfo")[["reTrms"]] <-
        lme4:::mkReTrms(lme4::findbars(methods::slot(RebelFitObj, "miscFitInfo")[["formula"]]), fr)
    }
    vcovBeta=lapply(vcovBetaAdj, function(x) x$vcovBeta)
    vcovBetaAdj=lapply(vcovBetaAdj, function(x) x$vcovBetaAdj)

    names(vcovBetaAdj)=names(vcovBeta)=gene_names


    methods::slot(RebelFitObj, "vcovBetaList")=vcovBeta
    methods::slot(RebelFitObj, "vcovBetaAdjList")=vcovBetaAdj
    RebelFitObj
}


.VCovAdj_1Gene=function(theta, sigma, devFun, Ztlist, sampleVariable,
                        subjectVariable, nObs, modelMatrix, pseudoBulk,
                        flist){
    vcovBeta=cov_beta=.get_covbeta(c(theta, sigma), devFun)
    G=.getGMats(Ztlist, sampleVariable, subjectVariable, nObs)
    if(!pseudoBulk){
        ggamma=.getGGamma(sigma, theta, sampleVariable, subjectVariable)
        SigmaG=NULL
    }else{
        SigmaG <- .rebelGetSigmaG(sigma, theta, sampleVariable, subjectVariable, G)
        ggamma=NULL
    }
    vcovAdj=.rebelVCovAdj_internal(Phi=vcovBeta, G=G, SigmaG=SigmaG,
                                ggamma=ggamma, modelMatrix=modelMatrix,
                                pseudoBulk=pseudoBulk,
                                subjectVariable=subjectVariable,
                                sampleVariable=sampleVariable,
                                flist=flist)

    list(vcovBeta=vcovBeta, vcovBetaAdj=vcovAdj)
}




## *NOTE: Internal Function taken from lme4 package
## Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects
##Models Using lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
.get_covbeta=function (varpar, devFun)
{
    nvarpar <- length(varpar)
    sigma <- varpar[nvarpar]
    theta <- varpar[-nvarpar]
    devFun(theta)
    df_envir <- environment(devFun)
    sigma^2 * Matrix::tcrossprod(df_envir$pp$RXi())
}


.varCov2Theta <- function(X, res_var){
    if(X==0) return(0)
    else{
        X <- t(chol(X/res_var))
        X_unpacked <- X[lower.tri(X, diag=TRUE)]
        return(X_unpacked)
    }

}




## Adapted from limma:::.get_SigmaG
.getGMats=function(Ztlist, sampleVariable, subjectVariable, nObs){
    G <- NULL

    ## Simplified this from pbkrtest code but now only works for random intercepts
    for (ss in c(sampleVariable,subjectVariable)) {
        ZZ <- Ztlist[[paste0(ss, ".(Intercept)")]]
        G <- c(G, list(Matrix::crossprod(ZZ, ZZ)))
    }

    ## Note: There is old code incorporating weights here that I might use later.
    G    <- c( G, list(Matrix::sparseMatrix(seq(1,nObs), seq(1,nObs), x=1 )))

    names(G)=c(sampleVariable, subjectVariable, "resid_var")
    G
}

## Adapted from limma:::.get_SigmaG
.getGGamma=function(sigma, theta, sampleVariable, subjectVariable){
    sc=sigma
    cnms=rep(list("(Intercept)"), length(theta))
    names(cnms)=c(sampleVariable,subjectVariable )
    nc=lengths(cnms)
    theta=theta
    GGamma <- lme4::mkVarCorr(sigma, cnms, nc, theta, names(cnms))
    ggamma <- NULL
    for (ii in c(sampleVariable, subjectVariable)) {
        Lii <- GGamma[[ii]]
        ggamma <- c(ggamma, Lii[lower.tri(Lii, diag = TRUE)])
    }
    ggamma <- c(ggamma, sigma^2)

    names(ggamma)=c(names(GGamma), "sigma2")
    return(ggamma)
}

## Adapted from limma:::.get_SigmaG
.rebelGetSigmaG=function (sigma, theta, sampleVariable, subjectVariable, G)
{


    ggamma=.getGGamma(sigma, theta, sampleVariable, subjectVariable)
    n.ggamma <- length(ggamma)

    Sigma <- ggamma[1] * G[[1]]
    for (ii in 2:n.ggamma) {
        Sigma <- Sigma + ggamma[ii] * G[[ii]]
    }
    #}

    SigmaG <- list(Sigma = Sigma, n.ggamma = n.ggamma)
    SigmaG



}

.rebelVCovAdj_internal=function ( Phi, G, SigmaG=NULL,ggamma=NULL, modelMatrix,
                               pseudoBulk, subjectVariable, sampleVariable, flist)
{

    if(!pseudoBulk){
        subj_vars=as.character(flist[[subjectVariable]])
        samp_vars=as.character(flist[[sampleVariable]])

        my_a=sum(ggamma)
        my_b=ggamma[[sampleVariable]]+ggamma[[subjectVariable]]
        my_c=ggamma[[subjectVariable]]

        SigmaInv_list=lapply(unique(subj_vars), function(subj){
            samps=unique(samp_vars[subj_vars==subj])
            samps_n=table(samp_vars)[samps]
            m=samps_n[1]
            ainv=.shermMorrInv(my_a,my_b,m)
            if(length(samps)>1){
                for (i in seq(2, length(samps))) {
                    n=samps_n[i]
                    short=i==2
                    my_sinv=.calc_sinv(my_a, my_b, my_c, ainv, n, short)
                    q1=.quad1(ainv, my_c, my_sinv, short)
                    q3=.quad3(my_sinv, my_c, ainv, short)
                    ainv=rbind(cbind(q1, t(q3)), cbind(q3, my_sinv))
                }
            }
            return(ainv)
        })
        SigmaInv=Matrix::bdiag(SigmaInv_list)
    }else{
        SigmaInv <- chol2inv(chol(Matrix::forceSymmetric(as(SigmaG$Sigma,
                                                            "matrix"))))
    }


    n.ggamma <- ifelse(!pseudoBulk, length(ggamma), SigmaG$n.ggamma)
    TT <- SigmaInv %*% modelMatrix
    HH_list <- HH <- OO <- vector("list", n.ggamma)


    if(!pseudoBulk){
        index_df=t(sapply(unique(subj_vars),  function(subj){
            maxval=max(which(subj_vars==subj))
            minval=min(which(subj_vars==subj))
            n_samps=length(unique(samp_vars[subj_vars==subj]))
            data.frame(min=minval, max=maxval, n_samps)
        }))
        for (ii in 1:n.ggamma) {

            ## This only works for sample/subject random intercepts
            gmat_name=names(G)[ii]
            if(gmat_name=="resid_var"){
                HH_list[[ii]]=SigmaInv_list
                HH=SigmaInv
            }else if(gmat_name==sampleVariable){
                HH_list[[ii]] <- lapply(1:length(SigmaInv_list), function(x){
                    if(index_df[[x,"n_samps"]]==1){
                        matrix(apply(SigmaInv_list[[x]], 2, sum),nrow(SigmaInv_list[[x]]),
                               nrow(SigmaInv_list[[x]]), byrow = TRUE)
                    }else{
                        min=index_df[[x,"min"]]
                        max=index_df[[x,"max"]]
                        Matrix::crossprod(G[[ii]][min:max,min:max],SigmaInv_list[[x]])
                    }
                })
                HH=Matrix::bdiag(HH_list[[ii]])
            }else{
                HH_list[[ii]]=lapply(SigmaInv_list,  function(x){
                    test=matrix(apply(x, 2, sum),nrow(x), nrow(x), byrow = TRUE)
                })
                HH=Matrix::bdiag(HH_list[[ii]])
            }
            OO[[ii]] <- HH %*% modelMatrix

        }
    } else{
        for (ii in 1:n.ggamma) {
            HH[[ii]] <- Matrix::crossprod(G[[ii]],SigmaInv)
            OO[[ii]] <- HH[[ii]]%*%modelMatrix
        }
    }

    ## Finding PP, QQ
    PP <- QQ <- NULL
    for (rr in 1:n.ggamma) {
        OrTrans <- Matrix::t(OO[[rr]])
        PP <- c(PP, list(Matrix::forceSymmetric(-1 * OrTrans %*% TT)))
        for (ss in rr:n.ggamma) {
            QQ <- c(QQ, list(OrTrans %*% SigmaInv %*% OO[[ss]]))
        }
    }

    ## Finding Ktrace
    Ktrace <- matrix(NA, nrow = n.ggamma, ncol = n.ggamma)
    if(!pseudoBulk){
        for (rr in 1:n.ggamma) {
            HrTrans <- lapply(HH_list[[rr]], Matrix::t)
            for (ss in rr:n.ggamma) {
                Ktrace[rr, ss] <- Ktrace[ss, rr] <- sum(
                    sapply(1:length(HrTrans), function(x){
                        sum(HrTrans[[x]]*HH_list[[ss]][[x]])
                    }))
            }
        }
    }else{
        for (rr in 1:n.ggamma) {
            HrTrans <- Matrix::t(HH[[rr]])
            for (ss in rr:n.ggamma) {
                Ktrace[rr, ss] <- Ktrace[ss, rr] <- sum(HrTrans *
                                                            HH[[ss]])
            }
        }
    }

    ## Finding IE2
    IE2 <- matrix(NA, nrow = n.ggamma, ncol = n.ggamma)
    for (ii in 1:n.ggamma) {
        Phi.P.ii <- Matrix::crossprod(Phi,PP[[ii]])
        for (jj in c(ii:n.ggamma)) {
            www <- .indexSymmat2vec(ii, jj, n.ggamma)
            IE2[ii, jj] <- IE2[jj, ii] <- Ktrace[ii, jj] - 2 *
                sum(Phi * QQ[[www]]) + sum(Phi.P.ii * (PP[[jj]] %*%
                                                           Phi))
        }
    }

    eigenIE2 <- eigen(IE2, only.values = TRUE)$values
     condi <- min(abs(eigenIE2))
    # WW <- if (condi > 3e-10)
    #     Matrix::forceSymmetric(2 * solve(IE2)) else Matrix::forceSymmetric(2 * ginv(IE2))
    WW <- tryCatch({
        if(condi > 3e-10){
            Matrix::forceSymmetric(2 * solve(IE2))
        }else Matrix::forceSymmetric(2 * ginv(IE2))
    }, error = function(e) {
        Matrix::forceSymmetric(2 * ginv(IE2))
    })
    UU <- matrix(0, nrow = ncol(modelMatrix), ncol = ncol(modelMatrix))
    for (ii in 1:(n.ggamma - 1)) {
        for (jj in c((ii + 1):n.ggamma)) {
            www <-.indexSymmat2vec(ii, jj, n.ggamma)
            UU <- UU + WW[ii, jj] * (QQ[[www]] - PP[[ii]] %*%
                                         Phi %*% PP[[jj]])
        }
    }
    UU <- UU + Matrix::t(UU)
    for (ii in 1:n.ggamma) {
        www <- .indexSymmat2vec(ii, ii, n.ggamma)
        UU <- UU + WW[ii, ii] * (QQ[[www]] - PP[[ii]] %*% Phi %*%
                                     PP[[ii]])
    }
    GGAMMA <- Phi %*% UU %*% Phi
    PhiA <- Phi + 2 * GGAMMA
    attr(PhiA, "P") <- PP
    attr(PhiA, "W") <- WW
    attr(PhiA, "condi") <- condi
    PhiA
}

## Functions to speed up matrix inverse
## See https://en.wikipedia.org/wiki/Invertible_matrix under "Blockwise inversion"
.shermMorrInv=function(my_a,my_b,n){
    val=-my_b/((my_a-my_b)*(n*my_b+my_a-my_b))
    mat=matrix(val, n,n)
    diag(mat)=val+1/(my_a-my_b)
    mat
}

.quad1=function(ainv, c1, sinv, short){
    if(short){
        ainv+sum(ainv[1,])^2*c1^2*ncol(sinv)*sum(sinv[,1])
    }else{
        mat=cbind(rowSums(ainv)*c1*sum(sinv[,1])*ncol(sinv)*c1)
        mat=Matrix::tcrossprod(mat[,rep(1,ncol(ainv))]%*%ainv)
        ainv+mat
    }

}

.quad3=function(sinv, c1, ainv, short){
    if(short){
        matrix(-sum(sinv[1,])*c1*sum(ainv[,1]),nrow=nrow(sinv),ncol=ncol(ainv))
    }else{
        mat=rbind(-sum(sinv[1,])*c1*colSums(ainv))
        mat[rep(1, nrow(sinv)),]
    }
}

.calc_sinv=function(d1, d2, c1, ainv,n, short){
    val=if(short){
        c1^2*sum(ainv[,1])*ncol(ainv)
    }else{
        sum(c1*colSums(ainv))*c1
    }
    d1_new=d1-val
    d2_new=d2-val
    .shermMorrInv(d1_new, d2_new,n)
}

## Internal helper function from pbkrtest package
## Ulrich Halekoh, Søren Højsgaard (2014). A Kenward-Roger Approximation and Parametric Bootstrap
## Methods for Tests in Linear Mixed Models - The R Package pbkrtest. Journal of Statistical
## Software, 59(9), 1-30. URL https://www.jstatsoft.org/v59/i09/.

.indexSymmat2vec = function (i, j, N)
{
    k <- if (i <= j) {
        (i - 1) * (N - i/2) + j
    }
    else {
        (j - 1) * (N - j/2) + i
    }
}
