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
```{r}
