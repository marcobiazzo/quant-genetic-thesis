library(Matrix)
library(lme4)
library(pedigreemm)
library(rlist)
library(stringr)
library(docstring)

fit_me = function(formula, 
                  data, 
                  label, 
                  pedigree = NULL, 
                  n = 1, 
                  optimizer = "nloptwrap")
{
  #' Fit mixed effect model
  #'
  #' Fits a mixed effect model with a non block covariance matrix
  #' @param formula A language type formula for the model to fit
  #' @param data A dataframe with columns the labels included in the formula 
  #' @param label Given that the dataset "data" is attached, label is the column 
  #' relative to the pedigree's label
  #' @param pedigree A pedigree having as labels the column "label" of the dataframe
  #' @param n The number of terms inside the mixed effect term relative to label
  #' @param optimizer The wrapper to use in the non-linear optimization step
  #' @return Returns a fitted lme4 object
  
  f = lFormula(formula,
               data = data, 
               control = lmerControl(optimizer = optimizer, 
                                     restart_edge = TRUE, 
                                     boundary.tol = 0.01))
  
  # Extract the relationship matrix A from the pedigree
  pos = unique(label)[order(unique(label))]
  A = if(is.null(pedigree)) diag(length(unique(label)))else getA(pedigree)[pos,pos]
  
  I = as(diag(n), "dgCMatrix" )
  nrest = dim(f$reTrms$Zt)[1]-dim(A)[1]*n
  
  Lt = chol(A)
  Mt = kronecker(Lt,I)
  
  # Modify the upper part of the matrix Zt
  f$reTrms$Zt[1:(n*(dim(A)[1])), ] = as(Mt%*%(f$reTrms$Zt[1:(n*(dim(A)[1])), ]), 
                                        "dgCMatrix")
  
  #Optimization steps
  devfun <- do.call(mkLmerDevfun,f)
  opt <- optimizeLmer(devfun,
                      control = lmerControl(optimizer = optimizer, 
                                                        restart_edge = TRUE, 
                                                        boundary.tol = 0.01,
                                                        calc.derivs = FALSE))
  
  fit <- mkMerMod(environment(devfun),
                  opt,
                  f$reTrms,fr = f$fr)
  
  return(fit)
}
