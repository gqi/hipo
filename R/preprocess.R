#' Preprocess Multi-trait Summary Statistics
#'
#' @description Preprocesses summary statistics for HIPO analysis. Performs SNP filtering and allele merge.
#'
#' @param sumstats A list of K containing summary statistics for K traits.
#' Each element is a data frame that contains at least 6 columns：\code{rsid}, \code{A1} (effect allele), \code{A2}, \code{N} (sample size), \code{z} (z-statistic), \code{pval}.
#' Three optional columns can be provided: \code{chr} (chromosome number), \code{bp} (base pair; physical position), \code{freqA1} (allele frequency of A1);
#' if provided, SNP filtering will be applied: (1) remove MHC region (26-34Mb of chromosome 6) (2) remove variants with MAF < \code{maf.thr} (see below for \code{maf.thr}).
#' @param maf.thr MAF threshold for quality control. SNPs with MAF < \code{maf.thr} are removed. Default 0.05, constrained between 0 and 0.5. Only effective when \code{freqA1} is present in \code{sumstats}.
#'
#' @return A list that contains
#'
#' \item{sumstats.all}{A data.frame containing the summary statistics of all individual traits.}
#' \item{traitvec}{A vector containing the trait names.}
#'
#' @import dplyr
#' @export
preprocess = function(sumstats, maf.thr = 0.05){
    # If sumstats is not named, name it.
    if (is.null(names(sumstats))){
        names(sumstats) = paste0('trait', 1:length(sumstats))
    }

    # Remove SNPs with sample size < 0.67 * (90 percentile)
    for (trait in names(sumstats)){
        sumstats[[trait]] = sumstats[[trait]] %>% filter((z^2<=80) & (N>0.67*quantile(N,0.9)))
        if (("chr"%in%colnames(sumstats[[trait]]) & ("bp"%in%colnames(sumstats[[trait]])))){
            sumstats[[trait]] = sumstats[[trait]] %>% filter(!(chr==6 & bp>26e6 & bp<34e6))
        }
        if ("freqA1"%in%colnames(sumstats[[trait]])){
            sumstats[[trait]] = sumstats[[trait]] %>% filter((freqA1>maf.thr) & (freqA1<1-maf.thr))
        }
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

    # If freqA1 is provided, align the alleles of the first trait such that A1 is the minor allele
    trait = names(sumstats)[1]
    if (paste0("freqA1.",trait) %in% colnames(sumstats.all)){
        tempA1 = sumstats.all[[paste0('A1.',trait)]]
        tempA2 = sumstats.all[[paste0('A2.',trait)]]
        tempz = sumstats.all[[paste0("z.",trait)]]
        # Retrieve indices that needs to be changed
        inds = sumstats.all[[paste0("freqA1.",trait)]] > 0.5
        tempA1[inds] = sumstats.all[[paste0('A2.',trait)]][inds]
        tempA2[inds] = sumstats.all[[paste0('A1.',trait)]][inds]
        tempz[inds] = -tempz[inds]
        sumstats.all = sumstats.all %>% mutate(A1 = tempA1, A2 = tempA2)
        sumstats.all[[paste0("z.",trait)]] = tempz
    } else{
        sumstats.all = sumstats.all %>% rename_(A1 = paste0('A1.',trait), A2 = paste0('A2.',trait))
    }

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
    sumstats.all = sumstats.all[,c("rsid", "A1", "A2", trait.spec)]

    return(list(sumstats.all = sumstats.all, traitvec = names(sumstats)))
}
