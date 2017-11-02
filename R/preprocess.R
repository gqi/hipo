#' Preprocess Multi-trait Summary Statistics
#'
#' Preprocess summary statistics for HIPO analysis. Perform SNP filtering and allele merge.
#'
#' @param sumstats A list of K containing summary statistics for K traits. Each element of the list contains 9 columns：rsid, A1 (effect allele), A2, N (sample size), z (z-statistic), pval, chr (chromosome number), bp (base pair; physical position), freqA1 (allele frequency of A1).
#' @param maf.thr MAF threshold in SNP QC. SNPs with MAF < maf.thr are removed. Default 0.05.
#'
#' @export
preprocess = function(sumstats, maf.thr = 0.05){
    # If sumstats is not named, name it.
    if (is.null(names(sumstats))){
        names(sumstats) = paste0('trait', 1:length(sumstats))
    }

    # Remove SNPs with sample size < 0.67 * (90 percentile)
    for (trait in names(sumstats)){
        sumstats[[trait]] = sumstats[[trait]] %>% filter(!(chr==6 & bp>26e6 & bp<34e6) & (z^2<=80) & (N>0.67*quantile(N,0.9)))
    }

    # Remove missing data, make A1 represent minor allele, remove rare variants
    # Merge data set: only keep SNPs that are available for all traits.
    for (trait in names(sumstats)){
        names(sumstats[[trait]])[2:ncol(sumstats[[trait]])] = paste0(names(sumstats[[trait]])[2:ncol(sumstats[[trait]])],'.',trait)
    }
    sumstats.all = sumstats[[1]]
    for (i in 2:length(sumstats)){
        sumstats.all = sumstats.all %>% inner_join(sumstats[[i]], by = 'rsid')
    }

    # Remove missing data
    sumstats.all = sumstats.all[complete.cases(sumstats.all),]

    # Align the alleles of the first trait such that A1 is the minor allele.
    trait = names(sumstats)[1]
    sumstats.all = sumstats.all %>% rename_(chr = paste0("chr.",trait), bp = paste0("bp.",trait))

    tempA1 = sumstats.all[[paste0('A1.',trait)]]
    tempA2 = sumstats.all[[paste0('A2.',trait)]]
    tempz = sumstats.all[[paste0("z.",trait)]]
    maf = sumstats.all[[paste0("freqA1.",trait)]]
    # Retrieve indices that needs to be changed
    inds = sumstats.all[[paste0("freqA1.",trait)]] > 0.5
    tempA1[inds] = sumstats.all[[paste0('A2.',trait)]][inds]
    tempA2[inds] = sumstats.all[[paste0('A1.',trait)]][inds]
    tempz[inds] = -tempz[inds]
    maf[inds] = 1 - maf[inds]

    sumstats.all = sumstats.all %>% mutate(A1 = tempA1, A2 = tempA2, maf = maf)
    sumstats.all[[paste0("z.",trait)]] = tempz
    sumstats.all = sumstats.all %>% filter(maf>maf.thr)

    # Align all alleles to the first trait (A1 is minor allele)
    for (i in 2:length(sumstats)){
        trait = names(sumstats)[i]
        inds = (sumstats.all[[paste0("A1.",trait)]] != sumstats.all[["A1"]])
        sumstats.all[[paste0("z.",trait)]][inds] = -sumstats.all[[paste0("z.",trait)]][inds]
    }

    trait.spec = NULL
    for (trait in names(sumstats)){
        trait.spec = c(trait.spec, paste0(c("N.", "z.", "pval."), trait))
    }
    sumstats.all = sumstats.all[,c("rsid", "A1", "A2", "chr", "bp", "maf", trait.spec)]

    return(list(sumstats.all = sumstats.all, traitvec = names(sumstats)))
}
