#' Heritability Informed Power Optimization (HIPO)
#'
#' @description Performs heritability informed power optimization (HIPO) to conduct powerful association testing across multiple traits.
#'
#' @param sumstats A list of K containing summary statistics for K traits.
#' Each element is a data frame that contains at least 6 columns：\code{rsid}, \code{A1} (effect allele), \code{A2}, \code{N} (sample size), \code{z} (z-statistic), \code{pval}.
#' Three optional columns can be provided: \code{chr} (chromosome number), \code{bp} (base pair; physical position), \code{freqA1} (allele frequency of A1);
#' if provided, SNP filtering will be applied: (1) remove MHC region (26-34Mb of chromosome 6) (2) remove variants with MAF < \code{maf.thr} (see below for \code{maf.thr}).
#' @param out.path The path where the LDSC intermediate outputs are stored.
#' @param HIPOD.num Number of HIPO components for which to calculate z statistics and p-values. Default 2.
#' @param ldsc.path The path to LDSC software.
#' @param python.path The path to Python, if you need to use a version other than the default one.
#' @param ldscore.path The path containing the LDSC LD score files. Default to be \code{file.path(ldsc.path,"eur_w_ld_chr/")}.
#' @param maf.thr MAF threshold for quality control. SNPs with MAF < \code{maf.thr} are removed. Default 0.05, constrained between 0 and 0.5. Only effective when \code{freqA1} is present in \code{sumstats}.
#' @param mergeallele Corresponds to the \code{--merge-allele} flag in LDSC, indicates whether to merge alleles to HapMap 3 SNPs. Default to be \code{TRUE}.
#' @param truncate Used only for high-dimensional phenotypes. If \code{NULL}, use the full \code{coherit.mat} and \code{ldscint.mat}; if a decimal, truncate the eigenvectors of ldscint.mat with eigenvalues < truncate*(max eigenvalue of ldscint.mat). Recommended 0.05.
#'
#' @details This function fits LD score regression to estimate the genetic covariance matrix and the covariance matrix of GWAS parameter estimates.
#' Suitable eigendecomposition is then performed to find the weights that lead to the largest average non-centrality parameter of the underlying test statistic (HIPO-D1),
#' as well as subsequent HIPO components. Z-statistics of HIPO components are computed for association testing.
#'
#' @return A list that contains
#'
#' \item{sumstats.all}{A data.frame containing the merged individual trait summary statistics and the z statistics (\code{z.HIPODx}) and p-values (\code{pval.HIPODx}) of HIPO components.}
#' \item{coherit.mat}{Genetic covariance matrix.}
#' \item{ldscint.mat}{Matrix of LD score regression intercepts.}
#' \item{eigenvalue}{Eigenvalues from HIPO eigendecomposition. Proportional to the average non-centrality parameter.}
#' \item{HIPOD.mat}{A matrix of which the columns are the HIPO components.}
#'
#' @references Qi, Guanghao, and Nilanjan Chatterjee. "Heritability Informed Power Optimization (HIPO) Leads to Enhanced Detection of Genetic Associations Across Multiple Traits." bioRxiv (2017): 218404.
#'
#' @import dplyr
#' @export
hipo = function(sumstats, out.path, HIPOD.num = 2, ldsc.path, python.path = NULL, ldscore.path = file.path(ldsc.path,"eur_w_ld_chr/"), maf.thr = 0.05, mergeallele = TRUE, truncate = NULL){
    temp = preprocess(sumstats, maf.thr)
    sumstats.all = temp[[1]]
    traitvec = temp[[2]]
    rm(temp)

    temp = covmat(sumstats.all, traitvec, out.path, ldsc.path = ldsc.path, python.path = python.path, ldscore.path = ldscore.path, mergeallele = mergeallele)
    coherit.mat = temp[[1]]
    ldscint.mat = temp[[2]]

    hipo.return = hipoeigen(sumstats.all, traitvec, coherit.mat, ldscint.mat, HIPOD.num = min(HIPOD.num,length(sumstats)), truncate = truncate)
    return(hipo.return)
}
