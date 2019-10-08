## reSNP
Github repository for robust and efficient single SNP association analysis 

## Install This Package from github

```{r}
library(devtools)
install_github("ykim03517/reSNP")  
```

## How To Use This Package 
```{r}
library(reSNP)
rec=singleSNP(1000, geno_eff=1, c=0, MAF=.1)
ans_BC=reSNP(rec, boot=TRUE, boot_Num=500, method_se="BC",
             alpha_CI=0.05)

add=singleSNP(1000, geno_eff=1, c=.5, MAF=.1)

## p-value, estimated genetic model, confidence interval
ans_pct=reSNP(add, boot=TRUE, boot_Num=500, method_se="percentile",
              alpha_CI=0.05)
ans_pct
```
