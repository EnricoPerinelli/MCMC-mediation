#' @title Computing mediation
#' @author Enrico Perinelli
#' @description This function is based on Preacher and Selig 2008 app, with some amelioration, such as (a) put standard errors, without needing to square it,
#' (b) setting seed, (c) Include the number and the "D" from the TECH3 Mplus output
#' @param a = path a (unstandardized)
#' @param b = path b (unstandardized)
#' @param se_a = standard error of path a
#' @param se_b = standard error of path b
#' @param cov_ab_number = number before D of the covariance between path a and b (tech3 output in Mplus)
#' @param D = number (and sign) after D of the the covariance between path a and b (tech3 output in Mplus); example: -3
#' @param rep = repetition (default = 20000)
#' @param my_seed = set your seed (default = 111)
#' @examples
#' med_mcmc(a= -0.367, b= 0.324, se_a = 0.130, se_b = 0.031, cov_ab_number = -0.368762, D = -03, rep = 20000, my_seed = 111)


med_mcmc <- function(a, b, se_a, se_b, cov_ab_number, D, rep = 20000, my_seed){

  require(MASS)

  # paths
  a=as.numeric(c(a))
  b=as.numeric(c(b))

  #asymptotic sampling variance for a and b
  se_a <- as.numeric(c(se_a))^2
  se_b <- as.numeric(c(se_b))^2

  # cov between parameters
  cov_ab_number <- as.numeric(c(cov_ab_number))
  D <- as.numeric(c(D))
  cov_ab_D <- 10^(D)
  cov_ab <- cov_ab_number*cov_ab_D

  #rep=20000
  conf=95
  pest=c(a,b)
  acov <- matrix(c(
    se_a, cov_ab,
    cov_ab, se_b
    ),2,2)

  my_seed <- 111
  set.seed(my_seed)

  mcmc <- mvrnorm(rep,pest,acov,empirical=FALSE)
  ab <- mcmc[,1]*mcmc[,2]
  low=(1-conf/100)/2
  upp=((1-conf/100)/2)+(conf/100)
  LL=quantile(ab,low)
  UL=quantile(ab,upp)
  LL4=format(LL,digits=4)
  UL4=format(UL,digits=4)
################################################
# The number of columns in the histogram can   #
# be changed by replacing 'FD' below with      #
# an integer value.                            #
################################################

ind_effect <- a*b
#print(ind_effect)

hist(ab,breaks='FD',col='skyblue',xlab=paste('Indirect effect = ', ind_effect,';', conf,'% Confidence Interval ','LL',LL4,'  UL',UL4),
     main='Distribution of Indirect Effect')
#print(ab)


cat(
  paste0('Indirect effect = ',
         ind_effect,'; ',
         conf,'% Confidence Interval = ',
         'LL ',LL4,', UL ',UL4
         )
)
}
