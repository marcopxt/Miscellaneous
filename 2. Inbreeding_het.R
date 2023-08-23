
#' Calculating inbreeding rates from a SNP markers matrix
#'
#' @param W matrix with the SNPs coded as 0,1,2 for AA,Aa,aa
#'
#' @return the inbreeding rate for the population
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#' @export

inbr_rate <- function(W){
  het=1-abs(W-1)
  fi=rowSums(het)/(ncol(W))
  inbreeding=1-fi
  return(inbreeding)
}
