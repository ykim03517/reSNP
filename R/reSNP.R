##' Testing association between a SNP and phenotype
##' @title testing association between a SNP and phenotype
##' @param data  a data set consists of phenotype, genotype under additive genetic model assumption, genotype expanded with grid interval, and covariate if applicable.
##' @param grid an interval of grid search (default=0.01; range: 0<grid<1). A lower value of grid gives denser intervals.
##' @param R a numeric value of the number of resampling replications.
##' @param testStat a character value of a statistical test statistic (default='score'). 'LRT' indicates likelihood ratio test; 'score' means score test.
##' @param out_type a character value of outcome type (default='D'). 'D' indicates dichotomous phenotype; 'C' represents continuouse phenotype.
##' @param boot a logical value indicating whether to use bootstrap for confidence interval of the estimated genetic model.
##' @param boot_Num a numeric value of bootstrap replications.
##' @param method_se a method to compute standard error (default='percentile'). Users must choose a method in c('standard', 'percentile', 'BC', 'BCa').
##' 'standard' calculates a standard error of the bootstrap samples. 'percentile' uses percentiles of the bootstrap distribution to define the percentile confidence limits. 'percentile' method corrects for non-normality of the estimator of interest.
##' 'BC' is a Bias-Corrected percentile method which corrects for narrowness bias of percentile CIs. The bias-correction factor isestimated using the inverse cdf of standard normal distribution. 'BCa' is a bias-corrected and accelerated
##' percentaile method which corrects for narrowness bias and for nonconstant standard error of the estimator of interest (acceleration). The acceleration factor estimates
##' the rate of change of the standard error of the estimator.
##' @param alpha_CI a significance level (default=0.05).
##' @return a list of p-value, estimated genetic model, a confidence interval of the estimated genetic model.
##' @author Yeonil Kim
##' @import RcppZiggurat
##' @import mvtnorm
##' @import lmtest
##' @import stats
##' @import RcppEigen
##' @useDynLib reSNP
##' @importFrom Rcpp sourceCpp
##' @export

reSNP <- function(data, grid = 0.01, R = 10^4, testStat = "score", out_type = "D", boot = FALSE, boot_Num = NULL, 
    method_se = NULL, alpha_CI = NULL) {
    
    # read data
    y = data[[1]]
    z = data[[3]]
    g = data[[4]]
    geno_grid = data[[5]]
    
    # grid intervals
    seqC = seq(0, 1, grid)
    
    # dimension
    n = length(y)
    m = length(seqC)
    
    # outcome type
    if (out_type == "C") {
        
        lm.null = lm(y ~ 1)
        
    } else if (out_type == "D") {
        
        lm.null = glm(y ~ 1, family = "binomial")
    }
    
    # score functions
    U_V = effScore(y, geno_grid, z, lm.null$fitted, boot = boot, boot_Num = boot_Num)
    
    # test statistics type
    
    if (testStat == "LRT") {
        
        for (j in seq_along(seqC)) {
            lm.full = glm(y ~ geno_grid[, j] + z, family = "binomial")
            LRT[j] = lrtest(lm.full, lm.null)[2, 4]
        }
        
        T_stat = max(LRT)
        est.loc = which.max(LRT)
        
    } else if (testStat == "score") {
        
        score_vec = colSums(U_V[[1]])^2/U_V[[2]]
        T_stat = max(score_vec)
        est.loc = which.max(score_vec)
        
    }
    
    # generating standard normal by RcppZiggurat package
    s = n * R
    G = arrayC(zrnormVec(zrnorm(s)), c(n, R))
    
    # matrix multiplication U %*% G
    mult = Multmat(t(U_V[[1]]), G)
    mult.sq = mult^2
    V = arrayC(rep(U_V[[2]], R), c(m, R))
    W = mult.sq/V
    W.star = colMax(W)
    
    cum.W = ecdf(W.star)
    pval = 1 - cum.W(T_stat)
    
    
    lm_est_gmodel = glm(y ~ geno_grid[, est.loc], family = "binomial")
    est_geno = summary(lm_est_gmodel)$coeff[2, 1]
    std_geno = summary(lm_est_gmodel)$coeff[2, 2]
    CI_OR = round(exp(c(est_geno - 1.96 * std_geno, est_geno + 1.96 * std_geno)), 2)
    names(CI_OR) = c("2.5%", "97.5%")
    
    OR_est_gmodel = exp(est_geno)
    
    if (!isFALSE(boot)) {
        
        out = CI_scoreBoot(data, U_V, boot_Num = boot_Num, method_se = method_se, alpha_CI = alpha_CI, 
            grid = grid)
        
        NOTE <- "odds ratio, p.value, estmated genetic model & confidence interval by bootstrap"
        structure(list(OR = OR_est_gmodel, CI_OR = CI_OR, p.value = pval, est_gene_model = out$est, CI_gene_model = out$CI, 
            method = out$method_se, note = NOTE))
        
    } else {
        
        NOTE <- "odds ratio, p.value, estmated genetic model"
        structure(list(OR = OR_est_gmodel, CI_OR = CI_OR, p.value = pval, est_gene_model = est, note = NOTE))
        
        
    }
}



