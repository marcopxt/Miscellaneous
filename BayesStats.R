#########################################
#
# Package: BayesianFAST
#
# File: fastStats.R
# Contains: fastStats
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2025
# Last update: Mar-2025
#
# License: GPL-3
#
##########################################

#' Calculates some statistics from a FA model outcomes generate from BGLR package
#'
#' @description
#'
#' @param Model Model outcomes from BGLR implementing a FA structure for genotypic effect.
#' 
#' 
#' @return A data frame with the overall performance and RMSE for the genotypes (and an index), reliability for the environment,
#' and variance explained per factor (in absolute values).
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @examples
#' \dontrun{
#' }
#'
#' @references \emph{Pérez, Paulino, and Gustavo de Los Campos. "Genome-wide regression and prediction with the BGLR statistical package." Genetics 198.2 (2014): 483-495.}
#'
#' @export



fastStats = function(Model = NULL){
  
  require(MASS)
  
  #-----  Bayesian model 
  nF <- Model$ETA[[1]]$Cov$nF
  uFAhat <- Model$ETAHat
  G0 <- Model$ETA[[1]]$Cov$Omega
  EnvLoadings <- Model$ETA[[1]]$Cov$W 
  psi <- Model$ETA[[1]]$Cov$PSI
  
  
  #----- Getting environmental loadings and genotypic scores from the Bayesian model
  
  loadings_pinv <- MASS::ginv(t(EnvLoadings)) 
  GenScores <- uFAhat %*% loadings_pinv
  colnames(GenScores) <- paste0('fa', 1:nF)
  
  #----- 1. Stats for the genotypes
  #----- RMSD
  RMSDj <- sapply(1:ncol(uFAhat), function(j) {
    (uFAhat[,j] - EnvLoadings[j,1] * GenScores[,1])^2
  })
  
  RMSD = apply(RMSDj, 1, FUN = function(.a) sum(.a)/ncol(uFAhat))
  
  #----- Overall performance
  if(nF == 1){
    OP = data.frame(Value =(mean(EnvLoadings[,1]) * GenScores[,1]))
  }else if(nF > 1){
    OP = data.frame(Value = t(apply(EnvLoadings[,1:nF], 2, mean) %*% t(GenScores[,1:nF])))                     #****
  }
  
  #----- index for OP and RMSD
  Index = 2*((OP$Value-mean(OP$Value))/sqrt(var(OP$Value))) - ((RMSD-mean(RMSD))/sqrt(var(RMSD)))
  
  StatsBGLR = data.frame(GenID = rownames(Model$ETAHat),
                         RMSD = RMSD,
                         OP = OP$Value,
                         Index = Index)
  
  rownames(StatsBGLR) <- NULL
  
  #----- 2. Reliability for each environment
  Reliability = data.frame(EnvID = names(Model$resCov$S0),
                           Reliability = apply(uFAhat, 2, FUN = function(.a) 1-(mean(.a^2)/mean(diag(G0))))) #G0 is the varCov among env matrix
  rownames(Reliability) <- NULL
  
  #----- 3. variance explained per factor
  var_rj_percent <- 100*(apply(EnvLoadings, 2, function(.a)  .a^2)) / (colSums(EnvLoadings^2) + psi)
  rownames(var_rj_percent) <- names(Model$resCov$S0)
  
  modOutput = list(Reliability = Reliability,
                   StatsBGLR = StatsBGLR,
                   VarExp = var_rj_percent)
  # outcomes
  return(modOutput)
  
}

