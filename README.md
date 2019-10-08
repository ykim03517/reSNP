## reSNP
Github repository for robust and efficient single SNP association analysis 

## Install This Package from github

```{r}
library(devtools)
install_github("ykim03517/reSNP")  
```

## How To Use This Package 
+ The following R codes provide a tutorial of how to produce results from this package using simulated data sets. 

```{r}
library(reSNP)

n <- 1000       # sample size 
eff <- 1        # genotype effect
alleleF <- 0.1  # minor allele frequency
sig <- 0.05     # significance level for confidence interval of estimated genetic model
num <- 500      # bootstrap replications for confidence interval of estimated genetic model

set.seed(0901)

## When recessive genetic model is true
gmodel=0
rec=singleSNP(n, geno_eff=eff, c=gmodel, MAF=alleleF)
ans_standard=reSNP(rec, boot=TRUE, boot_Num=num, method_se="standard", alpha_CI=sig)   

## When additve genetic model is true
gmodel=0.5
add=singleSNP(n, geno_eff=eff, c=gmodel, MAF=alleleF)
ans_BC=reSNP(add, boot=TRUE, boot_Num=num, method_se="BC", alpha_CI=sig)               # bias-correction method

## When dominant genetic model is true
gmodel=1
dom=singleSNP(n, geno_eff=eff, c=gmodel, MAF=alleleF)
ans_pct=reSNP(dom, boot=TRUE, boot_Num=num, method_se="percentile", alpha_CI=sig)      # percentile method 

## Please refer to the document of "reSNP" function by ?reSNP for details of method_se.

```
