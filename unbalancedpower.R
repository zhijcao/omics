#' power analysis
#'
#' @param k group number
#' @param n sample sizes for each group
#' @param a type I error
#' @param gm group means
#' @param gsd group stds
#'
#' @return total sample size, lambda, delta, alpha, power
#' @export
#'
#' @examples
#' unbalancedANOVA(k=3, n=c(4,3,5), a=0.05, gm=c(3,5,6), gsd=c(0.4,0.32,0.7))
#' \dontrun{
#' unbalancedANOVA(k=3, n=c(4,3,5), a=0.05, gm=c(3,5,6), gsd=c(0.4,0.32,0.7))
#' }
unbalancedANOVA <- function(k=k, n=n, a=0.05, gm=gm, gsd=gsd){

  sig2 <- sum((n-1)*gsd^2)/sum(n-1)
  lambda <- sum(n*(gm-mean(gm))^2)/sig2 #non-centrality parameter
  delta <- sqrt((sum((gm-mean(gm))^2)/k)/sig2) #effective size
  power <- 1-pf(qf(1-a, k-1, sum(n-1)), k-1, sum(n-1), lambda)
  return (c(totalSample=sum(n),lambda=lambda,delta=delta, alpha=a, power=power))
}
