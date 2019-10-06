##' Score functions
##' @title Score functions
##' @param pheno  a vector of phenotype.
##' @param geno_grid a n x m matrix of genotype (n=the number of sample size; m=the number of interval between 0 and 1).
##' @param covr a vector of covariate.
##' @param y_hat a vector of fitted phenotype values under the null hypothesis.
##' @param boot a logical value indicating whether to use bootstrap for confidence interval of the estimated genetic model.
##' @param boot_Num a numeric value of bootstrap replications.
##' @return a list of score functions, variances, bootstrap samples of the score functions, and bootstrap samples of the variances.
##' @author Yeonil Kim
##' @useDynLib reSNP
##' @importFrom Rcpp sourceCpp
##' @export


effScore <- function(pheno, geno_grid, covr = NULL, y_hat, boot = FALSE, boot_Num = NULL) {
    if (sum(covr) == 0) {
        covr = NULL
    }
    
    if (!is.null(dim(geno_grid))) {
        m = dim(geno_grid)[2]
    } else {
        m = 1
    }
    
    n = length(pheno)
    
    if (!is.null(covr)) {
        pheno = arrayC(rep(pheno, m), c(n, m))
        covr = arrayC(rep(covr, m), c(n, m))
        y_hat = arrayC(rep(y_hat, m), c(n, m))
        yp = (pheno - y_hat)
        pp = y_hat * (1 - y_hat)
        b1 = yp * geno_grid
        a1 = yp
        d1 = yp * covr
        ba = -pp * geno_grid
        bd = -pp * geno_grid * covr
        da = -pp * covr
        aa = -pp
        bb = -pp * geno_grid^2
        dd = -pp * covr^2
        colSba = colSums(ba)
        colSbd = colSums(bd)
        colSaa = sum(aa[, 1])
        colSda = sum(da[, 1])
        colSdd = sum(dd[, 1])
        I.1 = t(do.call(cbind, lapply(1:m, function(i) {
            -c(colSba[i], colSbd[i])
        })))
        I.2.inv = solve(-arrayC(c(colSaa, colSda, colSda, colSdd), dim = c(2, 2)))
        U.nus = do.call(cbind, lapply(1:m, function(i) {
            cbind(a1[, i], d1[, i])
        }))
        U = do.call(cbind, lapply(1:m, function(i) {
            b1[, i] - t(I.1[i, 1:2] %*% I.2.inv %*% t(U.nus[, (2 * i - 1):(2 * i)]))
        }))
        V = colSums(U^2)
        
    } else {
        
        pheno = arrayC(rep(pheno, m), c(n, m))
        y_hat = arrayC(rep(y_hat, m), c(n, m))
        yp = (pheno - y_hat)
        pp = y_hat * (1 - y_hat)
        b1 = yp * geno_grid
        a1 = yp
        ba = -pp * geno_grid
        aa = -pp
        bb = -pp * geno_grid^2
        colSba = colSums(ba)
        colSaa = sum(aa[, 1])
        I.1 = t(do.call(cbind, lapply(1:m, function(i) {
            colSba[i]
        })))
        I.2.inv = colSaa
        U.nus = a1[, 1]
        U = b1 - t((I.1/I.2.inv) %*% U.nus)
        V = colSums(U^2)
    }
    
    U_boot = NULL
    V_boot = NULL
    
    if (boot == TRUE) {
        seq = 1:n
        ind = sample(seq, n * boot_Num, replace = TRUE)
        U_boot0 = U[ind, ]
        
        U_boot = aperm(`dim<-`(t(U_boot0), list(m, n, boot_Num)), c(2, 1, 3))
        V_boot = do.call(rbind, lapply(1:boot_Num, function(b) {
            colSums(U_boot[, , b]^2)
        }))
        
    }
    
    return(list(U = U, V = V, U_boot = U_boot, V_boot = V_boot))
}


