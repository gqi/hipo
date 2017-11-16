#' Heritability Informed Power Optimization (HIPO)
#'
#' @param sumstats A list of K containing summary statistics for K traits. Each element of the list contains 9 columns：rsid, A1 (effect allele), A2, N (sample size), z (z-statistic), pval, chr (chromosome number), bp (base pair; physical position), freqA1 (allele frequency of A1).
#' @param maf.thr MAF threshold in SNP QC. SNPs with MAF < maf.thr are removed. Default 0.05.
#' @param HIPOD.num Number of HIPO components to calculate z stats and p-value for.
#' @param ldsc.path The directory of the LDSC software.
#' @param python.path The path of Python, if need to use a version other than the default one.
#' @param ldscore.path The path containing the LDSC LD score files. Default to be file.path(ldsc.path,"eur_w_ld_chr/").
#' @param mergeallele Whether to merge allele to HapMap 3 SNPs when using LDSC. Default to be TRUE.
#' @param truncate If NULL, use the full coherit.mat and ldscint.mat; if a decimal, truncate the eigenvectors of ldscint.mat with eigenvalues < truncate*(max eigenvalue of ldscint.mat). Recommended truncate = 0.05.
#'
#' @import dplyr
#' @export
hipo = function(sumstats, out.path, maf.thr = 0.05, HIPOD.num = 2, ldsc.path = "~/documents/software/ldsc", python.path = NULL, ldscore.path = file.path(ldsc.path,"eur_w_ld_chr/"), mergeallele = TRUE, truncate = NULL){
    temp = preprocess(sumstats, maf.thr)
    sumstats.all = temp[[1]]
    traitvec = temp[[2]]
    rm(temp)

    temp = covmat(sumstats.all, traitvec, out.path, ldsc.path = ldsc.path, python.path = python.path, ldscore.path = ldscore.path, mergeallele = mergeallele)
    coherit.mat = temp[[1]]
    ldscint.mat = temp[[2]]

    hipo.return = hipoeigen(sumstats.all, traitvec, coherit.mat, ldscint.mat, HIPOD.num = HIPOD.num, truncate = truncate)
    return(hipo.return)
}
