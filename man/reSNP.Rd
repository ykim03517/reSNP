% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reSNP.R
\name{reSNP}
\alias{reSNP}
\title{testing association between a SNP and phenotype}
\usage{
reSNP(data, grid = 0.01, R = 10^4, testStat = "score",
  out_type = "D", boot = FALSE, boot_Num = NULL, method_se = NULL,
  alpha_CI = NULL)
}
\arguments{
\item{data}{a data set consists of phenotype, genotype under additive genetic model assumption, genotype expanded with grid interval, and covariate if applicable.}

\item{grid}{an interval of grid search (default=0.01; range: 0<grid<1). A lower value of grid gives denser intervals.}

\item{R}{a numeric value of the number of resampling replications.}

\item{testStat}{a character value of a statistical test statistic (default='score'). 'LRT' indicates likelihood ratio test; 'score' means score test.}

\item{out_type}{a character value of outcome type (default='D'). 'D' indicates dichotomous phenotype; 'C' represents continuouse phenotype.}

\item{boot}{a logical value indicating whether to use bootstrap for confidence interval of the estimated genetic model.}

\item{boot_Num}{a numeric value of bootstrap replications.}

\item{method_se}{a method to compute standard error (default='percentile'). Users must choose a method in c('standard', 'percentile', 'BC', 'BCa').
'standard' calculates a standard error of the bootstrap samples. 'percentile' uses percentiles of the bootstrap distribution to define the percentile confidence limits. 'percentile' method corrects for non-normality of the estimator of interest.
'BC' is a Bias-Corrected percentile method which corrects for narrowness bias of percentile CIs. The bias-correction factor isestimated using the inverse cdf of standard normal distribution. 'BCa' is a bias-corrected and accelerated
percentaile method which corrects for narrowness bias and for nonconstant standard error of the estimator of interest (acceleration). The acceleration factor estimates
the rate of change of the standard error of the estimator.}

\item{alpha_CI}{a significance level (default=0.05).}
}
\value{
a list of p-value, estimated genetic model, a confidence interval of the estimated genetic model.
}
\description{
Testing association between a SNP and phenotype
}
\author{
Yeonil Kim
}
