#' Estimate covariance matrices
#'
#' Uses LD score regression to estimate the genetic covariance matrix and LDSC intercept matrix, which are needed for HIPO analysis.
#'
#' @param sumstats.all Summary statitics of all individual traits, first output of \code{preprocess()}.
#' @param traitvec Vector of trait names, second output of \code{preprocess()}.
#' @param out.path The path where the LDSC intermediate outputs are stored.
#' @param ldsc.path The path of the LDSC software.
#' @param python.path The path of Python, if you need to use a version other than the default one.
#' @param ldscore.path The path containing the LDSC LD score files. Default to be \code{file.path(ldsc.path,"eur_w_ld_chr/")}.
#' @param mergeallele Corresponds to the \code{--merge-allele} flag in LDSC, indicates whether to merge alleles to HapMap 3 SNPs. Default to be \code{TRUE}.
#'
#' @return A list containing
#'
#' \item{coherit.mat}{Genetic covariance matrix.}
#' \item{ldscint.mat}{LD score regression intercept matrix.}
#'
#' @import dplyr
#' @export
covmat = function(sumstats.all, traitvec, out.path, ldsc.path, python.path = NULL, ldscore.path = file.path(ldsc.path,"eur_w_ld_chr/"), mergeallele = TRUE){
    # Add anaconda python to path: for ldsc (no need to do this if the path is already set by the system).
    if (!is.null(python.path)){
        Sys.setenv(PATH = paste(python.path, Sys.getenv('PATH'), sep = ':'))
    }
    # Check if out.path exists
    if (!dir.exists(out.path)){
        dir.create(out.path)
    }

    # Further preprocess the data: randomly simulate signs for z.trait1 (LDSC requires median(z) is small).
    inds = sample(nrow(sumstats.all), round(nrow(sumstats.all))/2)
    tempA1 = sumstats.all$A1
    sumstats.all$A1[inds] = sumstats.all$A2[inds]
    sumstats.all$A2[inds] = tempA1[inds]
    for (trait in traitvec){
        sumstats.all[[paste0("z.",trait)]][inds] = -sumstats.all[[paste0("z.",trait)]][inds]
    }

    # Write files for LDSC
    for (trait in traitvec){
        sumstats.all %>%
            rename_(N = paste0("N.",trait), z = paste0("z.",trait), pval = paste0("pval.",trait)) %>%
            select(rsid, A1, A2, N, z, pval) %>%
            write_delim(file.path(out.path, paste0(trait,'_for_ldsc.txt')), delim = '\t')
    }

    # Call ldsc to estimate heritability, coheritability
    # Change summary level data to .sumstats format
    for (trait in traitvec){
        munge.sumstats.code = paste(file.path(ldsc.path,"munge_sumstats.py"),
                                    "--sumstats",
                                    file.path(out.path, paste0(trait,'_for_ldsc.txt')),
                                    "--out", file.path(out.path, paste0(trait,"_ldsc_format")))
        if (mergeallele==TRUE){
            munge.sumstats.code = paste(munge.sumstats.code, "--merge-alleles", file.path(ldsc.path,"w_hm3.snplist"))
        }
        system(munge.sumstats.code)
        print(paste("ldsc_sumstasts", trait))
    }

    # Fit LD score regression
    for (i in 1:(length(traitvec)-1)){
        traits = traitvec[i:(length(traitvec))]
        sumstats.files = paste(file.path(out.path, paste0(traits,"_ldsc_format.sumstats.gz")), collapse = ",")
        ldsc.code = paste(file.path(ldsc.path,"ldsc.py"),
                          "--rg", sumstats.files,
                          "--ref-ld-chr", ldscore.path,
                          "--w-ld-chr", ldscore.path,
                          "--out", file.path(out.path, paste0(paste(traits, collapse = "_"),"_ldsc_results")))
        system(ldsc.code)
        print(paste("ldsc", traits[1]))
    }

    # Retrieve heritability-coheritability matrix from the log files
    coherit.mat = matrix(NA, nrow = length(traitvec), ncol = length(traitvec))
    ldscint.mat = matrix(NA, nrow = length(traitvec), ncol = length(traitvec)) # Stores LDSC intercepts

    for (i in 1:(length(traitvec)-1)){
        traits = traitvec[i:length(traitvec)]
        logfile = readLines(file.path(out.path, paste0(paste(traits, collapse = "_"),"_ldsc_results.log")))
        # Get heritability information
        if (i==1){
            for (j in 1:length(traitvec)){
                if (j == 1){  # When j==1 the format is slightly different
                    ind = which(logfile == "Heritability of phenotype 1") # Two lines below the title is the heritability
                } else{
                    ind = which(logfile == paste0("Heritability of phenotype ", j, "/", length(traitvec)))
                }
                coherit.mat[j,j] = as.numeric(strsplit(logfile[ind+2], split = " ")[[1]][5])
                ldscint.mat[j,j] = as.numeric(strsplit(logfile[ind+5], split = " ")[[1]][2])
            }
        }

        # Get coheritability information
        gencov.inds = which(logfile == "Genetic Covariance")
        gencovs = sapply(strsplit(logfile[gencov.inds+2], split = " "), function(x) x[5])
        coherit.mat[i,(i+1):length(traitvec)] = as.numeric(gencovs)
        coherit.mat[(i+1):length(traitvec),i] = as.numeric(gencovs)

        ldscints = sapply(strsplit(logfile[gencov.inds+4], split = " "), function(x) x[2])
        ldscint.mat[i,(i+1):length(traitvec)] = as.numeric(ldscints)
        ldscint.mat[(i+1):length(traitvec),i] = as.numeric(ldscints)
    }

    return(list(coherit.mat = coherit.mat, ldscint.mat = ldscint.mat))
}
