# The first implementation of this function wants to do the same thing pedigreem did

library(Matrix)
library(lme4)
library(pedigreemm)
library(rlist)
library(stringr)
library(docstring)

# I pass "label" to the function, which is a column of the dataframe data and is the id number 
# ---which corresponds to the pedigree's label. Has to be pointed out that the label is relative 
# ---to one individual
# ACHTUNG: In this code we need to put the random effect relative to the label first! And in no other 
# ---random effects
# label has to be a string with the name of the column

fit_me = function(formula, data, label, pedigree = NULL, n = 1, optimizer = "nloptwrap")
{
  #' Fit mixed effect model
  #'
  #' Fits a mixed effect model with a non block covariance matrix
  #' @param formula A language type formula for the model to fit
  #' @param data A dataframe with columns the labels included in the formula 
  #' @param label Given that the dataset "data" is attached, label is the column relative to the pedigree?s label
  #' @param pedigree A pedigree having as labels the column "label" of the dataframe
  #' @param n The number of terms inside the mixed effect term relative to label
  #' @param optimizer The wrapper to use in the non-linear optimization step
  #' @return Returns a list with the fitted vector u* and the fitted lme4 object
  #' 
  
  
  # For now it is necessary to attach the dataframe before calling the function
  # ---and then inserting in label an attached column
  #if (typeof(label)!="character")
  #  print("The slot label must be a character")
  # comment(data) = label può tornare utile
  # attach(data)
  
  # Create a formula
  f = lFormula(formula, data = data, control = lmerControl(optimizer = optimizer, restart_edge = TRUE, boundary.tol = 0.01))# (per eliminare label provare: form$reTrms$flist[[1]])
  
  # Extract the relationship matrix A from the pedigree (re components are ordered lexicographically)
  pos = unique(label)[order(unique(label))]
  A = if(is.null(pedigree)) diag(length(unique(label)))else getA(pedigree)[pos,pos]
  
  # Volevo automatizzare la conta dei termini ad effetto misto ma risulta lungo, per adesso salto
  # Now we need to take the matrix constructed by lme4 and modify it to allow a full cov matrix
  # 1)First we construct the covariance matrix of the vector of random effect
  
  #terms = as.list(attr(terms(f),"variables"))
  
  # 1.1)We need to know how many terms there are in the random effect term of the discriminator
  #for (l in terms) # Attenzione l è proprio l'oggetto interno della lista, che è un'oggetto di tipo language
  #{
  #  if (length(grep(label,l))!=0)# Questo comando dovrebbe convertire l in lista e poi in character, qui dentro devo contare i termini del random effect(conto il numero di +)
  #    form = strsplit(toString(l))
  #}(Anche Ztlist potrebbe aiutare!)
  
  # This is the structure of the covariance matrix of a single random effect vector
  I = as(diag(n), "dgCMatrix" )
  nrest = dim(f$reTrms$Zt)[1]-dim(A)[1]*n
  
  # Compute the kronecker product 
  # M = if(n==1) A else kronecker(A, I)
  
  Lt = chol(A)
  
  Mt = kronecker(Lt,I)
  
  # Modify the upper part of the matrix Ztransposed
  f$reTrms$Zt[1:(n*(dim(A)[1])), ] = as(Mt%*%(f$reTrms$Zt[1:(n*(dim(A)[1])), ]), "dgCMatrix")
  
  #Optimization steps
  dfun <- do.call(mkLmerDevfun,f)
  opt <- optimizeLmer(dfun, control = lmerControl(optimizer = optimizer, restart_edge = TRUE, boundary.tol = 0.01))
  fit <- mkMerMod(environment(dfun), opt, f$reTrms,fr = f$fr)
  
  # The output of the function is a list made of the vector u* and the fitted lme4 object 
  # um = as.vector(t(Lt)%*%fit@u)
  
  return(fit)
}

# Alcune funzioni che potrebbero essere d'aiuto per trovare nuovi alg d'ottimizzazione

optimizeLmer(devfun,
             optimizer    = formals(lmerControl)$optimizer,
             restart_edge = formals(lmerControl)$restart_edge,
             boundary.tol = formals(lmerControl)$boundary.tol,
             start = NULL, verbose = 0L,
             control = list(), ...)

lmerControl(optimizer = "nloptwrap",
            
            restart_edge = TRUE,
            boundary.tol = 1e-5,
            calc.derivs = TRUE,
            use.last.params = FALSE,
            sparseX = FALSE,
            standardize.X = FALSE,
            ## input checking options
            check.nobs.vs.rankZ = "ignore",
            check.nobs.vs.nlev = "stop",
            check.nlev.gtreq.5 = "ignore",
            check.nlev.gtr.1 = "stop",
            check.nobs.vs.nRE= "stop",
            check.rankX = c("message+drop.cols", "silent.drop.cols", "warn+drop.cols",
                            "stop.deficient", "ignore"),
            check.scaleX = c("warning","stop","silent.rescale",
                             "message+rescale","warn+rescale","ignore"),
            check.formula.LHS = "stop",
            ## convergence checking options
            check.conv.grad     = .makeCC("warning", tol = 2e-3, relTol = NULL),
            check.conv.singular = .makeCC(action = "message", tol = formals(isSingular)$tol),
            check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
            ## optimizer args
            optCtrl = list(),
            mod.type = "lmer"
)
