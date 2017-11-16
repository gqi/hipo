# hipo

hipo is an R package that performs heritability informed power optimization (HIPO) for conducting multi-trait association analysis on summary level data.

### Installing the package.

Before using this package, you need to install the LDSC software that performs LD score regression to estimate heritability and genetic correlation. Click [here](https://github.com/bulik/ldsc) for LDSC. You must have Python and all the Python packages required by LDSC. LDSC recommends installing [Anaconda](https://store.continuum.io/cshop/anaconda/) which comes with all the required Python packages.

### Example 

Once you have installed and loaded the package, here's an example to guide you through the procedure.

1. Download and unzip the [Global Lipids Genetics Consortium (2013)](http://csg.sph.umich.edu/abecasis/public/lipids2013/) summary results from joint analysis of metabochip and GWAS data. Click on the link below to download
[LDL](http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_LDL.txt.gz) 
[HDL](http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_HDL.txt.gz)
[triglycerides](http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_TG.txt.gz)
[total cholesterol](http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_TC.txt.gz) 

2. Read the data
```{r}
rm(list=ls())
library(readr)
library(dplyr)
traitvec = c('LDL', 'HDL', 'TG', 'TC') # Following the order on the consortium website
sumstats = list()
for (trait in traitvec){
    sumstats[[trait]] = read_delim(paste0('blood_lipid_data/jointGwasMc_', trait, '.txt'), delim = '\t') %>% filter(rsid!='.') # Remove SNPs without an rs number.
}
```

3. Preprocess the data and transform them into data frames with 9 columns:
rsid: RS number of SNPs.
A1: effect allele.
A2: reference allele.
N: sample size.
z: z-statistic.
pval: p-value.
chr: chromosome.
bp: physical position.
freqA1: allele frequency of A1.
```{r}
for (trait in traitvec){
    sumstats[[trait]] = sumstats[[trait]] %>%
        mutate(position = strsplit(SNP_hg19, split = ':')) %>%
        mutate(chr = as.integer(gsub('chr', '',sapply(position, function(x) x[1]))), bp = as.integer(sapply(position, function(x) x[2]))) %>%
        mutate(z = beta/se, pval = 2*(1-pnorm(abs(z))) , A1 = toupper(A1), A2 = toupper(A2)) %>%
        select(rsid, A1, A2, N, z, pval, chr, bp, freqA1 = Freq.A1.1000G.EUR) %>%
    print(trait)
}
```

4. Implement HIPO
```{r}
library(hipo)
res = hipo(sumstats, out.path = "/Users/guanghaoqi/dropbox/High_D_trait/miscellaneous/test_data", maf.thr = 0.05, HIPOD.num = 4, ldsc.path = "~/documents/software/ldsc", python.path = '/Users/guanghaoqi/anaconda2/bin/', mergeallele = TRUE, truncate = NULL)

```

