##' Simulating Single SNP Data
##' @title Simulating Single SNP Data
##' @param size  a numeric value of sample size.
##' @param geno_eff a numeric value of genotype effect on phenotype.
##' @param grid an interval of grid search (default=0.01; range: 0<grid<1). A lower value of grid gives denser intervals.
##' @param c a numeric value of heterozygous genotype effect (0<=c<=1).
##' @param MAF a numeric value of minor allele frequency (0<MAF<=0.5).
##' @param pv a numeric value of prevalence of disease (default=0.5).
##' @param cov_eff a numeric value of covariate effect on phenotype (default=0.5).
##' @return a data set which consists of phenotype, genotype under additive assumption of genetic model, genotype under true genetic model, genotype matrix expanded with grid interval, and covariate.
##' @author Yeonil Kim
##' @import RcppZiggurat
##' @import mvtnorm
##' @useDynLib reSNP
##' @importFrom Rcpp sourceCpp
##' @export

singleSNP <- function(size, geno_eff, grid = 0.01, c, MAF, pv = 0.5, cov_eff = 0.5) {
    
    seqC = seq(0, 1, grid)
    
    ## dimension
    n <- size
    z = zrnormVec(zrnorm(n))
    g = rbinom(n, 2, (MAF))
    
    # min(count) >= N * MAF^2. Also, this gurantees 3 groups
    if (min(table(g)) < (n * MAF^2) || length(unique(g)) != 3) {
        while (min(table(g)) < (n * MAF^2) || length(unique(g)) != 3) {
            g = rbinom(n, 2, (MAF))
        }
    }
    g = sample(sample(g))
    
    # create x based on genotype and c
    geno_true = trueG(g, c)
    geno_grid = gridG(g, seqC)
    
    # create y
    a = log(pv/(1 - pv))
    p = exp(a + geno_eff * geno_true + cov_eff * z)/(1 + exp(a + geno_eff * geno_true + cov_eff * z))
    
    # generate y
    y = rbinom(n, 1, p)
    
    return(list(pheno = y, geno_true = geno_true, covr = z, geno_adtv = g, geno_grid = geno_grid))
}
