#' @title Estimation of mediation effects via Monte Carlo simulation
#' @author Enrico Perinelli
#' @description Based on Preacher and Selig (2008) app, with enhancements: 
#' (a) accepts standard errors directly, 
#' (b) allows setting a seed, 
#' (c) accepts covariance between a and b from Mplus TECH3 output, 
#' (d) allows automation for multiple mediations (e.g., cross-cultural analyses).
#' @param a = path a (unstandardized)
#' @param b = path b (unstandardized)
#' @param se_a = standard error of path a
#' @param se_b = standard error of path b
#' @param cov_ab_number = number before D of the covariance between path a and b (tech3 output in Mplus)
#' @param D = number (and sign) after D of the the covariance between path a and b (tech3 output in Mplus); example: -3
#' @param rep = repetition (default = 20000)
#' @param my_seed = set your seed (default = 111)
#' @param conf = set confidence intervals (default = 95)
#' @examples
#' med_mc(a= -0.367, b= 0.324, se_a = 0.130, se_b = 0.031, cov_ab_number = -3.7, D = -04, rep = 20000, my_seed = 111, conf = 95)


med_mc <- function(a, b, se_a, se_b, cov_ab_number = 0, D = 0, rep = 20000, my_seed = 111, conf = 95){

  require(MASS)

  # paths
  a=as.numeric(c(a))
  b=as.numeric(c(b))

  # asymptotic sampling variance for a and b
  se_a <- as.numeric(c(se_a))^2
  se_b <- as.numeric(c(se_b))^2

  # cov between parameters
  cov_ab_number <- as.numeric(c(cov_ab_number))
  D <- as.numeric(c(D))
  cov_ab_D <- 10^(D)
  cov_ab <- cov_ab_number*cov_ab_D
  
  # set repetitions and confidence intervals
  #rep=20000 (now default)
  #conf=95   (now default)
  
  # build asymptotic covariance matrix
  pest=c(a,b)
  acov <- matrix(c(
    se_a, cov_ab,
    cov_ab, se_b
    ),2,2)

  # set seed
  set.seed(my_seed)

  # run Monte Carlo simulation
  mcmc <- mvrnorm(rep, pest, acov, empirical = FALSE)
  ab <- mcmc[,1]*mcmc[,2]
  
  # compute confidence interval
  low=(1-conf/100)/2
  upp=((1-conf/100)/2)+(conf/100)
  LL=quantile(ab,low)
  UL=quantile(ab,upp)
  LL4 <- round(LL, 4)
  UL4 <- round(UL, 4)

  # indirect effect
  ind_effect <- a * b

  # plot
  hist(ab, breaks='FD', col='skyblue',
       xlab = paste0('Indirect effect = ', round(ind_effect, 4),
                     '; ', conf, '% CI: [', round(LL, 4), ', ', round(UL, 4), ']'),
       main='Distribution of Indirect Effect')

  # output
  cat(
    paste0('Indirect effect = ', round(ind_effect, 4),
           '; ', conf, '% Confidence Interval = [', LL4, ', ', UL4, ']')
  )
}
