# gwrpv v1.0 AUTHORS: Gregory Connor, Michael O'Neill : GC authored original code in RATS, MO'N translated into R

#' gwrpvr: A package for calculating Genome-Wide Regression P-Values (gwrpv) in R
#'
#' Computes the sample probability value (p-value) for the estimated coefficient
#' from a standard genome-wide univariate regression.
#' It computes the exact finite-sample p-value under the assumption that
#' the measured phenotype (the dependent variable in the regression)
#' has a known Bernoulli-normal mixture distribution.
#'
#' The gwrpvr package provides two functions:
#' gwrpv and gwrpv_batch.
#'
#'
#' @docType package
#' @name gwrpvr
NULL

#' regresults: sample data
#'
#' A sample dataset of input regression results based on machine-level accurate cumulative normal values.
#' Rather than just typing in a few digits of the 2.5%, .5% and 2.5%x10^(-6) normal p-value cutoffs,
#' the norminverse function in RATS was used to create sample-case betas which are exact
#'
#' @name regresults
#' @format csv format file with 4 variables (beta, n0, n1, n2) and 120 rows
NULL

#' Genome-Wide Regression P-Value (gwrpv) in R
#'
#' Computes the sample probability value (p-value) for the estimated coefficient
#' from a standard genome-wide univariate regression.
#' It computes the exact finite-sample p-value under the assumption that
#' the measured phenotype (the dependent variable in the regression)
#' has a known Bernoulli-normal mixture distribution.
#'
#' @param beta the beta being tested
#' @param n0 number of major allele homozygotes
#' @param n1 number of major allele heterozygotes
#' @param n2 number of minor allele zygotes
#' @param mua parameter of the mixture distribution, can be any real number
#' @param siga parameter of the mixture distribution, can be any real number
#' @param mub parameter of the mixture distribution, can be any real number
#' @param sigb parameter of the mixture distribution, can be any real number
#' @param pa parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param pb parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param logdelta must be in log base 10 format, with default value set to -16
#' @param lognearnorm must be in log base 10 format, with default value set to -5
#' @param logtopsum must be in log base 10 format, with default value set to 8
#'
#' @return gwrpv returns a list containing:
#' \describe{
#'   \item{$pvalue}{p-value of a two-sided hypothesis test for a true coefficient of zero}
#'   \item{$skew}{skewness}
#'   \item{$kurt}{kurtosis of the coefficient estimate under assumed model}
#'   \item{$skiptype}{type of trimming/skip which took place (zero means no trimming)}
#'   \item{$totnobs}{total number of observations}
#'   \item{$loopruns}{number of sums in the main computation for each regression case}
#' }
#'
#'
#'                       .
#'
#' @examples
#' beta <- 6.05879
#' n0 <- 499
#' n1 <- 1
#' n2 <- 0
#' mua <- 13.87226
#' siga <- 2.58807
#' mub <- 4.62829
#' sigb <- 2.51803
#' pa <- 0.96544
#' pb <- 0.03456    # alternatively: pb <- 1.0 - pa
#' gwrpv(beta,n0,n1,n2,mua,siga,mub,sigb,pa,pb)
#'
#' # note default values have been used for the trim parameters above
#' # in the following example we explicitly set the trim parameters
#' #
#' g <- gwrpv(beta,n0,n1,n2,mua,siga,mub,sigb,pa,pb,logdelta=-16,lognearnorm=-5,logtopsum=8)
#' g$pvalue

gwrpv <- function(beta,n0,n1,n2,mua,siga,mub,sigb,pa,pb,logdelta=-16,lognearnorm=-5,logtopsum=8) {

  # process the parameters...

  # beta is the coefficient and n0, n1 and n2 are the number of each type of SNP observations
  # default values are set for logdelta=-16,lognearnorm=-5,logtopsum=8
  # to override these defaults the user must explicitly reset these values when calling the function
  # pa+pb should sum to 1.0
  if (pa + pb != 1) stop("Parameter input error: pa+pb!=1")

  # remove the log base 10 scale
  # unless the user has zero'ed the parameters out
  if (logtopsum != 0) {
      topsum <- 10^logtopsum
  }
  if ( (logdelta != 0) && (logdelta > -50) ) {
    downtrim <- (1/6) * (10^logdelta)
  } else {
    downtrim <- 0
  }
  if ( (lognearnorm != 0) && (lognearnorm > -50) ) {
    nearnorm <- 10^lognearnorm
  } else {
    nearnorm <- 0
  }

  # subtract the expected value from the dependent variable
  meany <- pa * mua + pb * mub
  mua <- mua - meany
  mub <- mub - meany

  # dimension output vectors to match numtests & initialise to zero's
  pvalue <- 0
  loopruns <- 0
  totnobs <- 0
  skew <- 0
  kurt <- 0
  skiptype <- 0

  # Set the total number of observations in the regression equal to the sum of the three types of observations and save in the output vector totnobs
  totnobs <- n0 + n1 + n2

  # can we skip the full iteration over n0, n1 and n2?  and calculate skewness, kurtosis, sigbeta
  ctn <- close_to_normal(totnobs, n0, n1, n2, pa, pb, mua, mub, siga, sigb, beta, nearnorm)
  skew <- ctn$skewbeta
  kurt <- ctn$kurtbeta
  sigbeta <- ctn$sigbeta
  skipiter <- ctn$skipiter


  # create proportions of the three types of SNP observations
  prop <- rep(0, 3)
  prop[1] <- n0 / totnobs
  prop[2] <- n1 / totnobs
  prop[3] <- n2 / totnobs

  # create a zero mean explanatory variable from the SNP data set
  rawx <- rep(0, 3)
  rawx[1] <- 0
  rawx[2] <- 1
  rawx[3] <- 2
  meanrawx <- prop[1] * rawx[1] + prop[2] * rawx[2] + prop[3] * rawx[3]
  x <- rep(0, 3)
  x[1] <- rawx[1] - meanrawx
  x[2] <- rawx[2] - meanrawx
  x[3] <- rawx[3] - meanrawx

  # find the sum of squares of the explanatory variable
  sumsqx <- (n0 * (x[1]^2) + n1 * (x[2]^2) + n2 * (x[3]^2))

  # initialise pvalue at zero
  if (skipiter < 2) {
    pvalue <- 0
  }

  # find the upper and lower limits in order to drop the outcomes which have cumulative probability close to zero set the near-zero cumulative
  # probability limit for trimming

  # initialize index limits at zero
  low <- rep(0, 3)
  high <- rep(0, 3)
  high[1] <- n0
  high[2] <- n1
  high[3] <- n2

  hl <- highlow(downtrim, high[1], pa, pb)
  high[1] <- hl[1]
  low[1] <- hl[2]

  if (skipiter == 1) {
  	high[1] <- 0
  	low[1] <- 0
	}


  hl <- highlow(downtrim, high[2], pa, pb)
  high[2] <- hl[1]
  low[2] <- hl[2]

  hl <- highlow(downtrim, high[3], pa, pb)
  high[3] <- hl[1]
  low[3] <- hl[2]

  # check that the number of loop runs is not excessive
  loopruns <- (high[1] - low[1] + 1) * (high[2] - low[2] + 1) * (high[3] - low[3] + 1)

  # if the loop runs IS excessive, don't calculate and set
  # the pvalues of these tests as -999 to flag this special case
  if (loopruns > topsum) {
    pvalue <- -999
    skipiter <- 2
  }

  if (loopruns <= topsum)
  {
    # assign zero values in case any loops are null (zero cases)
    n0a <- 0
    n1a <- 0
    n2a <- 0
    vary <- pa * (mua^2 + siga^2) + pb * (mub^2 + sigb^2) - (pa * mua + pb * mub)^2
    # loop over all (untrimmed) outcomes of the three binomial random variables
    pvalue <- loop_calc_pvalue(low[1], high[1], low[2], high[2], low[3], high[3], n0a, n1a, n2a, n0, n1, n2,
                               pa, pb, x, mua, mub, sumsqx, siga, sigb, vary, beta, skipiter, pvalue)

    # change to two-tailed pvalue units
    if (beta > 0) {
      pvalue <- 1 - pvalue
    }

  }  # END of skipped runs: if(loopruns[i] <= topsum){
  skiptype <- skipiter

  # return a list of outputs
  list(pvalue = pvalue, skew = skew, kurt = kurt, skiptype = skiptype, totnobs = totnobs, loopruns = loopruns)
}  #END gwrpv() function



#' Batch computation of a list of pvalues of GWA regression beta statistics using a bernoulli-normal mixture distribution
#'
#' @param regresults a list of four lists.
#' \describe{
#'   \item{$beta}{the list of betas being tested}
#'   \item{$n0}{the list of major allele homozygotes}
#'   \item{$n1}{the list of major allele heterozygotes}
#'   \item{$n2}{the list of minor allele zygotes}
#' }
#' @param mua parameter of the mixture distribution, can be any real number
#' @param siga parameter of the mixture distribution, can be any real number
#' @param mub parameter of the mixture distribution, can be any real number
#' @param sigb parameter of the mixture distribution, can be any real number
#' @param pa parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param pb parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param logdelta must be in log base 10 format, with default value set to -16
#' @param lognearnorm must be in log base 10 format, with default value set to -5
#' @param logtopsum must be in log base 10 format, with default value set to 8
#'
#' @return gwrpv_batch returns a list of lists containing the lists:
#' \describe{
#'   \item{$pvalue}{p-value of a two-sided hypothesis test for a true coefficient of zero}
#'   \item{$skew}{skewness}
#'   \item{$kurt}{kurtosis of the coefficient estimate under assumed model}
#'   \item{$skiptype}{type of trimming/skip which took place (zero means no trimming)}
#'   \item{$totnobs}{total number of observations}
#'   \item{$loopruns}{number of sums in the main computation for each regression case}
#' }
#'
#'
#'                       .
#'
#' @examples
#' beta <- c(6.05879, -6.05879, 2.72055, -2.72055, 1.93347,
#'          -1.93347, 0.88288, -0.88288, 4.28421, -4.28421)
#' n0 <- c(499, 499, 495, 495, 490, 490, 451, 451, 998, 998)
#' n1 <- c(1, 1, 5, 5, 10, 10, 48, 48, 2, 2)
#' n2 <- c(0, 0, 0, 0, 0, 0, 1, 1, 0, 0)
#' myregresults <- list(beta = beta, n0 = n0, n1 = n1, n2 = n2)
#' mua <- 13.87226
#' siga <- 2.58807
#' mub <- 4.62829
#' sigb <- 2.51803
#' pa <- 0.96544
#' pb <- 1.0 - pa
#' gwrpv_batch(myregresults,mua,siga,mub,sigb,pa,pb)
#' # store results in a user-defined variable g
#' g <- gwrpv_batch(myregresults,mua,siga,mub,sigb,pa,pb,logdelta=-16,lognearnorm=-4,logtopsum=8)
#' g$pvalue

gwrpv_batch <- function(regresults,mua,siga,mub,sigb,pa,pb,logdelta=-16,lognearnorm=-5,logtopsum=8) {

  # process the parameters...

  # beta is the coefficient and n0, n1 and n2 are the number of each type of SNP observations
  # default values are set for logdelta=-16,lognearnorm=-5,logtopsum=8 and numtests=1
  # to override these defaults the user must explicitly reset these values when calling the function
  # pa+pb should sum to 1.0
  if (pa + pb != 1) stop("Parameter input error: pa+pb!=1")

  # remove the log base 10 scale
  # unless the user has zero'ed the parameters out
  if (logtopsum != 0) {
    topsum <- 10^logtopsum
  }
  if ( (logdelta != 0) && (logdelta > -50) ) {
    downtrim <- (1/6) * (10^logdelta)
  } else {
    downtrim <- 0
  }
  if ( (lognearnorm != 0) && (lognearnorm > -50) ) {
    nearnorm <- 10^lognearnorm
  } else {
    nearnorm <- 0
  }

  # remove the expected value from the dependent variable
  meany <- pa * mua + pb * mub
  mua <- mua - meany
  mub <- mub - meany

  # read in the coefficient and the number of each type of SNP observations
  beta <- regresults$beta
  n0 <- regresults$n0
  n1 <- regresults$n1
  n2 <- regresults$n2

  numtests <- length(beta)
  # dimension output vectors to match numtests & initialise to zero's
  pvalue <- rep(0, numtests)
  loopruns <- rep(0, numtests)
  totnobs <- rep(0, numtests)
  skew <- rep(0, numtests)
  kurt <- rep(0, numtests)
  skiptype <- rep(0, numtests)


  # loop over the number, "numtests", of regressions being tested
  # "i" is an index counting which test number we are looking at now
  #cat("looping over each regression being tested...")
  for (i in 1:numtests) {

    #cat(i, ", ")
    # Set the total number of observations in the regression equal to the sum of the three types of observations and save in the output vector totnobs
    totnobs[i] <- n0[i] + n1[i] + n2[i]

    # can we skip the full iteration over n0, n1 and n2?  and calculate skewness, kurtosis, sigbeta
    ctn <- close_to_normal(totnobs[i], n0[i], n1[i], n2[i], pa, pb, mua, mub, siga, sigb, beta[i], nearnorm)
    skew[i] <- ctn$skewbeta
    kurt[i] <- ctn$kurtbeta
    sigbeta <- ctn$sigbeta
    skipiter <- ctn$skipiter


    # create proportions of the three types of SNP observations
    prop <- rep(0, 3)
    prop[1] <- n0[i] / totnobs[i]
    prop[2] <- n1[i] / totnobs[i]
    prop[3] <- n2[i] / totnobs[i]

    # create a zero mean explanatory variable from the SNP data set
    rawx <- rep(0, 3)
    rawx[1] <- 0
    rawx[2] <- 1
    rawx[3] <- 2
    meanrawx <- prop[1] * rawx[1] + prop[2] * rawx[2] + prop[3] * rawx[3]
    x <- rep(0, 3)
    x[1] <- rawx[1] - meanrawx
    x[2] <- rawx[2] - meanrawx
    x[3] <- rawx[3] - meanrawx

    # find the sum of squares of the explanatory variable
    sumsqx <- (n0[i] * (x[1]^2) + n1[i] * (x[2]^2) + n2[i] * (x[3]^2))

    # initialise pvalue at zero
    if (skipiter < 2) {
      pvalue[i] <- 0
    }

    # find the upper and lower limits in order to drop the outcomes which have cumulative probability close to zero set the near-zero cumulative
    # probability limit for trimming

    # initialize index limits at zero
    low <- rep(0, 3)
    high <- rep(0, 3)
    high[1] <- n0[i]
    high[2] <- n1[i]
    high[3] <- n2[i]

    hl <- highlow(downtrim, high[1], pa, pb)
    high[1] <- hl[1]
    low[1] <- hl[2]

  	if (skipiter == 1) {
	    high[1] <- 0
	    low[1] <- 0
	  }

    hl <- highlow(downtrim, high[2], pa, pb)
    high[2] <- hl[1]
    low[2] <- hl[2]

    hl <- highlow(downtrim, high[3], pa, pb)
    high[3] <- hl[1]
    low[3] <- hl[2]

    # check that the number of loop runs is not excessive
    loopruns[i] <- (high[1] - low[1] + 1) * (high[2] - low[2] + 1) * (high[3] - low[3] + 1)

    # if the loop runs IS excessive, don't calculate and set
    # the pvalues of these tests as -999 to flag this special case
    if (loopruns[i] > topsum) {
      pvalue[i] <- -999
      skipiter <- 2
    }

    if (loopruns[i] <= topsum)
    {
      # assign zero values in case any loops are null (zero cases)
      n0a <- 0
      n1a <- 0
      n2a <- 0
      vary <- pa * (mua^2 + siga^2) + pb * (mub^2 + sigb^2) - (pa * mua + pb * mub)^2
      # loop over all (untrimmed) outcomes of the three binomial random variables
      pvalue[i] <- loop_calc_pvalue(low[1], high[1], low[2], high[2], low[3], high[3], n0a, n1a, n2a, n0[i], n1[i], n2[i],
                                    pa, pb, x, mua, mub, sumsqx, siga, sigb, vary, beta[i], skipiter, pvalue[i])

      # change to two-tailed pvalue units
      if (beta[i] > 0) {
        pvalue[i] <- 1 - pvalue[i]
      }

    }  # END of skipped runs: if(loopruns[i] <= topsum){
    skiptype[i] <- skipiter
  }  # END loop over the number of regressions being tested: for(i in 1:numtests){

  list(pvalue = pvalue, skew = skew, kurt = kurt, skiptype = skiptype, totnobs = totnobs, loopruns = loopruns)

}  #END gwrpv_batch() <-function


# Helper Functions ---------------------------------------------------------------

#' This is a CLT-linked run-time control.
#'
#' If the number of observations is large enough that a normality approximation
#' holds for the y average across the major homozygote subsample, then the code skips the
#' time-consuming loop over n0, n1 and n2 and and uses the normal approximation for the average y
#' for the major homozygote subsample. The remaining loop is only over n1 and n2.
#' The only new input/output variables are input lognearnorm (the magnitude of maximum allowed
#' tolerance (in log 10 format) for the sum of squared deviation of skewness and
#' kurtosis from their normal values and output stopiter
#' (a zero if the code does not mandate a stop to the iterative estimation and a one if it does).
#' The input variable lognearnorm has a default value set so that users only have to enter it if they want to over-ride the default value.
#'
#' @param totnobs the sum of n0, n1, n2
#' @param n0 the major allele homozygotes
#' @param n1 the major allele heterozygotes
#' @param n2 the minor allele zygotes
#' @param pa parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param pb parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param mua parameter of the mixture distribution, can be any real number
#' @param mub parameter of the mixture distribution, can be any real number
#' @param siga parameter of the mixture distribution, can be any real number
#' @param sigb parameter of the mixture distribution, can be any real number
#' @param beta the beta from the regression being tested
#' @param nearnorm must be in log base 10 format, with default value set to -5
#'
#'
#' @return list(skewbeta = skewbeta, kurtbeta = kurtbeta, sigbeta = sigbeta, skipiter = skipiter)
#'
#'
close_to_normal <- function(totnobs, n0, n1, n2, pa, pb, mua, mub, siga, sigb, beta, nearnorm) {
    # compute skewness and kurtosis of the beta coefficient
    skipiter <- 0

	  # This is a CLT-linked run-time control.
	  # If the number of observations is large enough that a normality approximation
	  # holds for the y average across the major homozygote subsample, then the code skips the
	  # time-consuming loop over n0, n1 and n2 and and uses the normal approximation for the average y
	  # for the major homozygote subsample. The remaining loop is only over n1 and n2.The only new input/output variables are input lognearnorm (the magnitude of maximum allowed tolerance (in log 10 format) for
    # the sum of squared deviation of skewness and kurtosis from their normal values and output stopiter (a zero if the code does not mandate a stop to
    # the iterative estimation and a one if it does).  The input variable lognearnorm has a default value set so that users only have to enter it if they
    # want to over-ride the default value.

    # Step 1 Create standardized exogenous and endogenous variables

    mx <- (n1 + 2 * n2)/totnobs
    varx <- (n1 + 4 * n2)/totnobs - mx^2
    sigx <- (varx)^(1/2)
    x0 <- -mx/sigx
    x1 <- (1 - mx)/sigx
    x2 <- (2 - mx)/sigx
    muy <- pa * mua + pb * mub
    sigy <- (pa * (mua^2 + siga^2) + pb * (mub^2 + sigb^2) - muy^2)^(0.5)
    muastn <- (mua - muy)/sigy
    mubstn <- (mub - muy)/sigy
    sigastn <- siga/sigy
    sigbstn <- sigb/sigy

    # Step 2 Skewness and kurtosis of the endogenous variable, of the estimated regression coefficient beta,
	  # and over the average y over the major homozygote subsample.

    skewy <- pa * (muastn^3 + 3 * muastn * sigastn^2) + pb * (mubstn^3 + 3 * mubstn * sigbstn^2)
    kurty <- pa * (muastn^4 + 6 * muastn^2 * sigastn^2 + 3 * sigastn^4) + pb * (mubstn^4 + 6 * mubstn^2 * sigbstn^2 + 3 * sigbstn^4)
    skewbeta <- ((1/totnobs))^((3/2)) * (n0 * x0^3 + n1 * x1^3 + n2 * x2^3) * skewy
    kurtbeta <- ((1/totnobs))^2 * ((n0 * x0^4 + n1 * x1^4 + n2 * x2^4) * kurty +
                                     3 * (n0 * (n0 - 1) * x0^4 + n1 * (n1 - 1) * x1^4 + n2 * (n2 - 1) * x2^4) +
                                     6 * (n0 * n1 * x0^2 * x1^2 + n0 * n2 * x0^2 * x2^2 + n1 * n2 * x1^2 * x2^2))
	  skewmjh <- (1/n0)^(.5) * skewy
	  kurtmjh <- 3 + (1/n0) * (kurty - 3)


    # If the distribution of the y average over the major homozygote subsample is close to normal
	  # at tolerance level nearnorm, then the code will use normal-based approximate and skip the loop over n0.
    sigbeta <- sigy/(sigx * totnobs^0.5)
    # central limit theorem z statistic
    cltzstat <- beta/sigbeta
    # cat('beta sigbeta sigy sigx cltzstat: ',beta,sigbeta,sigy,sigx,cltzstat,'\n')

    sqdif <- (skewmjh)^2 + (kurtmjh - 3)^2

    if (sqdif < nearnorm) {
        skipiter <- 1
    } else if (sqdif >= nearnorm) {
        skipiter <- 0
    }

    return(list(skewbeta = skewbeta, kurtbeta = kurtbeta, sigbeta = sigbeta, skipiter = skipiter))
}

#' highlow()
#'
#' If possible, trim the upper and lower bounds
#'
#'
#' @param downtrim lower bound
#' @param n upper bound
#' @param pa parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param pb parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#'
#' @return c(lhigh, llow))  # return the new upper and lower bounds
#'
#'
highlow <- function(downtrim, n, pa, pb) {
    # instantiate local variables which will be returned by highlow function
    llow <- 0  # initial lower bound
    lhigh <- n  # initial upper bound

    # instantiate local utility variables used to calculate return variables
    lcumprob <- 0  #cumulative prob
    rhscumprob <- 0

    for (index in 0:n) {
        remain <- n - index

        # can we improve the lower bound?
        lnprobpart1 <- lgamma(n + 1) - lgamma(index + 1) - lgamma(remain + 1)
        lnprobpart2 <- index * log(pa) + remain * log(pb)
        lcumprob <- lcumprob + exp(lnprobpart1 + lnprobpart2)

        if (lcumprob < downtrim) {
            llow <- index  # we can improve the lower bound, set to new value
        }

        # can we improve the upper bound?
        lnprobpart1 <- lgamma(n + 1) - lgamma(index + 1) - lgamma(remain + 1)
        lnprobpart2 <- index * log(pb) + remain * log(pa)
        rhscumprob <- rhscumprob + exp(lnprobpart1 + lnprobpart2)

        if (rhscumprob < downtrim) {
            lhigh <- remain  # we can improve the upper bound, set to new value
        }

    }
    return(c(lhigh, llow))  # return the new upper and lower bounds
}

#' calc_pvalue()
#'
#' calculate the pvalue : called from loop_calc_pvalue()
#'
#' @param n0a outer loop index
#' @param n1a middle loop index
#' @param n2a inner loop index
#' @param n0 the major allele homozygotes
#' @param n1 the major allele heterozygotes
#' @param n2 the minor allele zygotes
#' @param pa parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param pb parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param x a zero mean explanatory variable from the SNP data set
#' @param mua parameter of the mixture distribution, can be any real number
#' @param mub parameter of the mixture distribution, can be any real number
#' @param sumsqx sum of the squares of x
#' @param siga parameter of the mixture distribution, can be any real number
#' @param sigb parameter of the mixture distribution, can be any real number
#' @param vary vary <- pa*(mua^2+siga^2)+pb*(mub^2+sigb^2)-(pa*mua+pb*mub)^2
#' @param beta the beta from the regression being tested
#' @param skipiter flag to determine if we can skip some calculations
#' @param pvalue the input pvalue prior to calculating new improved pvalue
#'
#' @return pvalue
#'
#'
#'
calc_pvalue <- function(n0a, n1a, n2a, n0, n1, n2, pa, pb, x, mua, mub, sumsqx, siga, sigb, vary, beta, skipiter, pvalue) {
    n0b <- n0 - n0a
    n1b <- n1 - n1a
    n2b <- n2 - n2a


    lnprobpart1 <- lgamma(n0 + 1) - lgamma(n0a + 1) - lgamma(n0b + 1) + lgamma(n1 + 1) - lgamma(n1a + 1) - lgamma(n1b + 1) +
                   lgamma(n2 + 1) - lgamma(n2a + 1) - lgamma(n2b + 1)
    lnprobpart2 <- (n0a + n1a + n2a) * log(pa) + (n0b + n1b + n2b) * log(pb)
    lnprob <- lnprobpart1 + lnprobpart2

    # compute the conditional mean and volatility of beta for a given (n0a,n1a,n2a)
    condmeanbeta <- ((n0a * x[1] + n1a * x[2] + n2a * x[3]) * mua + (n0b * x[1] + n1b * x[2] + n2b * x[3]) * mub)/sumsqx
    condvarbeta <- (((n0a * (x[1]^2) + n1a * (x[2]^2) + n2a * (x[3]^2)) * (siga^2) + (n0b * (x[1]^2) + n1b * (x[2]^2) + n2b * (x[3]^2)) * (sigb^2)))/(sumsqx^2)

	  if (skipiter == 1) {
	    lnprobpart1 <-  lgamma(n1 + 1) - lgamma(n1a + 1) - lgamma(n1b + 1) +
	                    lgamma(n2 + 1) - lgamma(n2a + 1) - lgamma(n2b + 1)
      lnprobpart2 <- (n1a + n2a) * log(pa) + (n1b + n2b) * log(pb)
      lnprob <- lnprobpart1 + lnprobpart2
	    # compute the conditional mean and volatility of beta for a given (n1a,n2a)
      condmeanbeta <- ((n1a * x[2] + n2a * x[3]) * mua + (n1b * x[2] + n2b * x[3]) * mub)/sumsqx
      condvarbeta <- (((n1a * (x[2]^2) + n2a * (x[3]^2)) * (siga^2) + (n1b * (x[2]^2) + n2b * (x[3]^2)) * (sigb^2)))/(sumsqx^2)
	    condvarbeta <- condvarbeta + n0 * x[1]^2 * vary/sumsqx^2
	  }

  	condvolbeta <- condvarbeta^0.5
    zstat <- (beta - condmeanbeta)/condvolbeta
    pvalue <- pvalue + exp(lnprob) * pnorm(zstat)
    return(pvalue)
}

#' loop_calc_pvalue()
#'
#' calls calc_pvalue()
#'
#' @param lowone lower bound outer loop
#' @param highone upper bound outer loop
#' @param lowtwo lower bound middle loop
#' @param hightwo upper bound middle loop
#' @param lowthree lower bound inner loop
#' @param highthree upper bound inner loop
#' @param n0a outer loop index
#' @param n1a middle loop index
#' @param n2a inner loop index
#' @param n0 the major allele homozygotes
#' @param n1 the major allele heterozygotes
#' @param n2 the minor allele zygotes
#' @param pa parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param pb parameter of the mixture distribution, a real number between zero and one with pa+pb=1
#' @param x a zero mean explanatory variable from the SNP data set
#' @param mua parameter of the mixture distribution, can be any real number
#' @param mub parameter of the mixture distribution, can be any real number
#' @param sumsqx sum of the squares of x
#' @param siga parameter of the mixture distribution, can be any real number
#' @param sigb parameter of the mixture distribution, can be any real number
#' @param vary vary <- pa*(mua^2+siga^2)+pb*(mub^2+sigb^2)-(pa*mua+pb*mub)^2
#' @param beta the beta from the regression being tested
#' @param skipiter flag to determine if we can skip some calculations
#' @param pvalue the input pvalue prior to calculating new improved pvalue
#'
#' @return pvalue
#'
#'
loop_calc_pvalue <- function(lowone, highone, lowtwo, hightwo, lowthree, highthree,
                             n0a, n1a, n2a, n0, n1, n2, pa, pb, x, mua, mub, sumsqx, siga, sigb, vary, beta, skipiter, pvalue) {
    for (n0a in lowone:highone) {
        for (n1a in lowtwo:hightwo) {
            for (n2a in lowthree:highthree) {
                pvalue <- calc_pvalue(n0a, n1a, n2a, n0, n1, n2, pa, pb, x, mua, mub, sumsqx, siga, sigb, vary, beta, skipiter, pvalue)
            }  #end for n2a
        }  #end for n1a
    }  #end for n0a
    return(pvalue)
}






