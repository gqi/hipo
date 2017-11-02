#' @title HIPO eigendecomposition
#'
#' @description This function performs HIPO eigendecomposition and calculate the z statistics and p-values of HIPO components.
#'
#' @param sumstats.all a list of K that contains summary statistics for K traits
#' @param traitvec Vector of trait names, second output of preprocess().
#' @param coherit.mat Genetic covariance matrix.
#' @param ldscint.mat LD score regression intercept matrix.
#' @param HIPOD.num Number of HIPO components to calculate z stats and p-value for.
#' @param truncate If NULL, use the full coherit.mat and ldscint.mat; if a decimal, truncate the eigenvectors of ldscint.mat with eigenvalues < truncate*(max eigenvalue of ldscint.mat). Recommended truncate = 0.05.
#'
#' @return A list that contains
#'
#' \item{sumstats.all}{A data.frame containing the merged individual trait summary statistics and the z statistic (z.HIPODx) and p-value (pval.HIPODx) of HIPO components.}
#' \item{coherit.mat}{Heritability-coheritability matrix.}
#' \item{ldscint.mat}{Matrix of LD score regression intercepts.}
#' \item{eigenvalue}{Eigenvalues from HIPO eigendecomposition. Proportional to the average non-centrality parameter.}
#' \item{HIPOD.mat}{A matrix of which the columns are the HIPO components.}
#'
#' @export
hipoeigen = function(sumstats.all, traitvec, coherit.mat, ldscint.mat, HIPOD.num = 2, truncate = NULL){

    Nsqrt = sapply(traitvec, function(x) sqrt(median(sumstats.all[[paste0("N.",x)]])))
    Nadjust.mat = Nsqrt %*% t(Nsqrt)
    betahatcov.mat.scaled = ldscint.mat/Nadjust.mat*(prod(Nsqrt)^(2/length(Nsqrt)))

    if (is.null(truncate)){
        R = chol(betahatcov.mat.scaled)
        eigen.ratio.mat = eigen(solve(t(R)) %*% coherit.mat %*% solve(R))
        HIPOD.mat = solve(R) %*% eigen.ratio.mat$vectors # This HIPOD.mat is not orthogonal.
    } else{
        temp = eigen(betahatcov.mat.scaled)
        Uk = temp$vectors[,temp$values > truncate*max(temp$values)]
        coherit.Uk = t(Uk) %*% coherit.mat %*% Uk
        betahatcov.scaled.Uk = t(Uk) %*% betahatcov.mat.scaled %*% Uk

        R = chol(betahatcov.scaled.Uk)
        eigen.ratio.mat = eigen(solve(t(R)) %*% coherit.Uk %*% solve(R))
        HIPOD.mat = Uk %*% solve(R) %*% eigen.ratio.mat$vectors
    }

    # Calculate the z stats and p-values for genetic principal components
    for (i in 1:HIPOD.num){
        HIPOD = HIPOD.mat[,i]
        znum = (as.matrix(sumstats.all[,paste0('z.',names(sumstats))]) / sqrt(as.matrix(sumstats.all[,paste0('N.',names(sumstats))]))) %*% HIPOD # linear combination of betahat
        # ldscint.mat and SNP specific sample size are used to allow for different variance for different SNPs
        # Do not use betahatcov.mat.scaled
        zdenom = (matrix(rep(diag(ldscint.mat), each = nrow(sumstats.all)), nrow = nrow(sumstats.all))/as.matrix(sumstats.all[,paste0('N.',names(sumstats))])) %*% HIPOD^2
        for (k in 1:(length(sumstats)-1)){
            for (l in (k+1):length(sumstats))
                zdenom = zdenom + 2*HIPOD[k]*HIPOD[l]*ldscint.mat[k,l]/sqrt(sumstats.all[[paste0("N.",names(sumstats)[k])]]*sumstats.all[[paste0("N.",names(sumstats)[l])]])
        }
        zdenom = sqrt(zdenom)
        sumstats.all[,paste0("z.HIPOD",i)] = as.vector(znum / zdenom)
        sumstats.all[,paste0("pval.HIPOD",i)] = 2 * (1-pnorm(abs(sumstats.all[[paste0("z.HIPOD",i)]])))
    }

    colnames(HIPOD.mat) = paste0("D",1:ncol(HIPOD.mat))

    return(list(sumstats.all = sumstats.all,
                coherit.mat = coherit.mat,
                ldscint.mat = ldscint.mat,
                eigenvalue = eigen.ratio.mat$values,
                HIPOD.mat = HIPOD.mat))
}

