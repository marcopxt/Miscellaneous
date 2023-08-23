#----------------
# Relationship Matrix calculations
#
# marco antonio peixoto
# marco.peixotom@gmail.com


#' Function to measure the f (inbreeding coefficient), based on
#' Xiang et al. 2016
#' 
#' @param W matrix of SNP for all genotypes coded 0, 1, 2.
#'
#' @return vector of inbreeding coefficient per individual
#' 
#' @export

getInbreCoef = function(W){
  M <- W
  M[M == 2] <- 0
  h = rowSums(M) 
  f = 1-(h/ncol(M)) # f = 1-h/N
  return(f)
}

#' Function to Building the matrices A (van Raden), DVit (Vitezica), and DSu (Su)
#' 
#'  
#' @param M matrix of SNPs coded as 0, 1, and 2, for BB, Bb, and bb, respectively
#' @param type the type of the matrix to be created 
#' @param thresh threshold for missing SNP values filtering
#' @param maf.limite control for minor allele frequency in the SNP dataset
#'
#' @returns A matrix, NSNPs by NSNPs for the specified kernel
#'
#' @export

relMatrix = function(M, type, maf.limite = 0.01, thresh = 0.8){
  
  cat("The method will return a", type, "matrix \n")
  
  W = M
  
  cat("Initial data: \n")
  cat("\tNumber of Individuals:", nrow(W), "\n")
  cat("\tNumber of Markers:", ncol(W), "\n")
  
  # Filtering for missing values
  missing <- apply(W, 2, function(x) sum(is.na(x))/nrow(W))
  missing.who <- missing <= thresh
  
  if (any(missing.who)) {
    idx.rm <- which(missing.who)
    W <- W[, idx.rm, drop = FALSE]
  } else {
    cat("\tNo SNPs with missing data, missing threshold of = ", 
        thresh, "\n")
  }
  
  # Filtering for maf
  W = maf_filter(W, maf.limite)
  
  cat("Final data: \n")
  cat("\tNumber of Markers:", ncol(W), "\n")
  
  # Allele frequency
  af = apply(W, 2, function(x) mean(x, na.rm = T)/2)
  
  # Minor allele frequency
  maf = ifelse(af < 0.5, af, 1-af)
  
  # FreqP matrix
  P <- matrix(rep(af,nrow(W)),ncol=ncol(W),byrow=T)
  
  
  if (type %in% c("A")){
    freq_p = maf;freq_q = 1-maf;Z = W
    
    twopq <- 2*(t(freq_p)%*%(freq_q))
    
    Z = Z - 2*P
    Z[is.na(Z)] <- 0
    
    MT <- (tcrossprod(Z,Z))/as.numeric(twopq)
    return(MT)
    
    
  }else if(type %in% c("DVit")){
    fq_mp = P; fq_mq = 1-P;freq_p = maf;freq_q = 1-maf;Z = W
    
    Z[is.na(Z)]<-0
    Z = (Z == 0)*-2*(fq_mp^2) + (Z==1)*2*fq_mq*fq_mp + (Z==2)*-2*(fq_mq^2)
    
    varD<-sum((2*(freq_p*freq_q))^2) #snp variance
    
    MT <- tcrossprod(Z,Z)/varD
    return(MT)
    
  }else if(type %in% c("DSu")) {
    fq_mp = P; fq_mq = 1-P;freq_p = maf;freq_q = 1-maf;Z = W
    
    twopq = 2*fq_mp*fq_mq
    
    Z[Z == 2] <- 0
    Z = Z - twopq
    Z[is.na(Z)]<-0
    
    MT <- tcrossprod(Z,Z)/sum(twopq[1,]*(1-twopq[1,]))
    
  }
}


#' Function to calculate the minor allele frequency of markers
#'
#'
#' @param SNP_Mat matrix of SNPs coded as 0,1,2 to A1A1, A1A2, and A2A2, respectively
#' @param thresh threshold for missing SNP values filtering
#'
#' @returns return the SNP matrix filtered
#' @export

maf_filter = function(SNP_Mat, thresh){
  MAF <- apply(SNP_Mat, 2, function(x) {
    AF <- mean(x, na.rm = T)/2
    MAF <- ifelse(AF > 0.5, 1 - AF, AF)
  })
  snps.low <- MAF < thresh
  if (any(snps.low)) {
    idx.rm <- which(snps.low)
    SNP_Mat <- SNP_Mat[, -idx.rm, drop = FALSE]
  }
  return(SNP_Mat)
}
