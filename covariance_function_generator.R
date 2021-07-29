# genetic covariance function
creator_gen_cov_function = function(alpha_matrix, type = "fourier")
{
  gen_func = function(t1,t2)
  {
    if(length(t1)==1){
      result = basis_matrix(t1,type = type, r = dim(alpha_matrix)[1])%*%
        alpha_matrix%*%
        (basis_matrix(t2,type = type, r = dim(alpha_matrix)[1]))
      result = as.vector(diag(result))
      
    }else{
    result = basis_matrix(t1,type = type, r = dim(alpha_matrix)[1])%*%
           alpha_matrix%*%
           t(basis_matrix(t2,type = type, r = dim(alpha_matrix)[1]))
    result = as.vector(diag(result))
    }
    return(result)
  }
  return(gen_func)
}