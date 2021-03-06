% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwrpv.R
\name{gwrpv_batch}
\alias{gwrpv_batch}
\title{Batch computation of a list of pvalues of GWA regression beta statistics using a bernoulli-normal mixture distribution}
\usage{
gwrpv_batch(regresults, mua, siga, mub, sigb, pa, pb, logdelta = -16,
  lognearnorm = -5, logtopsum = 8)
}
\arguments{
\item{regresults}{a list of four lists.
\describe{
  \item{$beta}{the list of betas being tested}
  \item{$n0}{the list of major allele homozygotes}
  \item{$n1}{the list of major allele heterozygotes}
  \item{$n2}{the list of minor allele zygotes}
}}

\item{mua}{parameter of the mixture distribution, can be any real number}

\item{siga}{parameter of the mixture distribution, can be any real number}

\item{mub}{parameter of the mixture distribution, can be any real number}

\item{sigb}{parameter of the mixture distribution, can be any real number}

\item{pa}{parameter of the mixture distribution, a real number between zero and one with pa+pb=1}

\item{pb}{parameter of the mixture distribution, a real number between zero and one with pa+pb=1}

\item{logdelta}{must be in log base 10 format, with default value set to -16}

\item{lognearnorm}{must be in log base 10 format, with default value set to -5}

\item{logtopsum}{must be in log base 10 format, with default value set to 8}
}
\value{
gwrpv_batch returns a list of lists containing the lists:
\describe{
  \item{$pvalue}{p-value of a two-sided hypothesis test for a true coefficient of zero}
  \item{$skew}{skewness}
  \item{$kurt}{kurtosis of the coefficient estimate under assumed model}
  \item{$skiptype}{type of trimming/skip which took place (zero means no trimming)}
  \item{$totnobs}{total number of observations}
  \item{$loopruns}{number of sums in the main computation for each regression case}
}


                      .
}
\description{
Batch computation of a list of pvalues of GWA regression beta statistics using a bernoulli-normal mixture distribution
}
\examples{
beta <- c(6.05879, -6.05879, 2.72055, -2.72055, 1.93347,
         -1.93347, 0.88288, -0.88288, 4.28421, -4.28421)
n0 <- c(499, 499, 495, 495, 490, 490, 451, 451, 998, 998)
n1 <- c(1, 1, 5, 5, 10, 10, 48, 48, 2, 2)
n2 <- c(0, 0, 0, 0, 0, 0, 1, 1, 0, 0)
myregresults <- list(beta = beta, n0 = n0, n1 = n1, n2 = n2)
mua <- 13.87226
siga <- 2.58807
mub <- 4.62829
sigb <- 2.51803
pa <- 0.96544
pb <- 1.0 - pa
gwrpv_batch(myregresults,mua,siga,mub,sigb,pa,pb)
# store results in a user-defined variable g
g <- gwrpv_batch(myregresults,mua,siga,mub,sigb,pa,pb,logdelta=-16,lognearnorm=-4,logtopsum=8)
g$pvalue
}

