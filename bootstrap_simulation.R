library(MASS)
library(ggplot2)
library(ggpubr)
library(Metrics)
library(lme4)
library(pedigreemm)
library(mvnfast)
library(plotly)

TRFUN25PUP4 <- read.delim("TRFUN25PUP4.DAT", header=FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")

df = TRFUN25PUP4

# set the number in each random effect vector
n_of_terms_gen = 3
n_of_terms_nongen = 3

iterations = 400

# remapping of x in [-1,1]
df$x = df$x*(1/13)-1
head(df)

# matrix of the fourier basis
Z = basis_matrix(df$x, r = 3, type = "legendre")
df = cbind(df,Z)
names(df)[6:8] = c("z1", "z2", "z3")
head(df)

# pedigree costruction from the dataset df
posFirstUniqueId = which((!duplicated(df$id))==TRUE)

sr = c( rep(NA, 100),
        rep(NA, 9900),
        rep(NA, 1500))
sr[unique(df$id)] = (df$sire)[posFirstUniqueId]

dm = c( rep(NA, 100),
        rep(NA, 9900),
        rep(NA, 1500))
dm[unique(df$id)] = (df$dam)[posFirstUniqueId]

lbl = c(1:11500)
pedigree = pedigree(sr, dm, lbl)

# get the matrix A 
position = unique(df$id)[order(unique(df$id))]
A = getA(pedigree)[position,position]
dim(A)

df$id2 = df$id

# create the matrices from which to generate data 

# small genetic covariance matrix 
small_mat_gen = matrix(c(400, 40, 0,
                         40, 200, 40,
                         0, 40, 200), nrow = 3, byrow = T)
small_mat_gen

# small environmental covariance 
small_mat_env = matrix(c(400, 40, 0,
                         40, 200, 40,
                         0, 40, 200), nrow = 3, byrow = T)
small_mat_env

# variance of the error sigma^2
s_2_residuals = 25


# big genetic covariance  
c_gen = as(kronecker(A, small_mat_gen),"dgCMatrix")

# big env covariance
c_env = as(kronecker(diag(length(unique(df$id))), small_mat_env),"dgCMatrix")

# general covariance
c = as(bdiag(c_gen,c_env),"dgCMatrix")

# length of the vector u
ul = n_of_terms_gen*length(unique(df$id)) + n_of_terms_nongen*length(unique(df$id)) 
n = dim(df)[1] 

# compute the fixed effect
g = lm(trait ~ -1+ df$z1+ df$z2+ df$z3, data = df)

length_error_1 = n
mu_u = rep(0, ul)
mu_e = rep(0, length_error_1)
c_e = as(s_2_residuals*diag(length_error_1), "dgCMatrix")

f = df$z1*g$coefficients[1] + df$z2*g$coefficients[2]+ df$z3*g$coefficients[3]

u = as.vector(rmvn(1, mu_u, c, ncores = 2))
e = as.vector(rmvn(1, mu_e, c_e, ncores = 2))
form = lFormula(df$trait ~ -1 + df$z1 +df$z2 + df$z3 +(-1 + df$z1 +df$z2 + df$z3 |df$id) + (-1 + df$z1 +df$z2+ df$z3 |df$id), 
                data = df)
df$trait = as.vector(f + t(form$reTrms$Zt)%*%u + e)

formula = df$trait ~ -1 + df$z1 + df$z2 + df$z3 +  ( -1 + df$z1 + df$z2 + df$z3 |df$id) + ( -1 + df$z1 + df$z2 + df$z3  |df$id2)

fit <- fit_me( formula,
               data = df,
               label = df$id, 
               pedigree = pedigree, 
               n = n_of_terms_gen,
               optimizer = 'nloptwrap')

extract_variance = function(fit_object)
{
  result = as.matrix(as.data.frame(VarCorr(fit_object))["vcov"])
  return(result)
}

# extract the estimate genetic and env. variance
# genetic
data_gen = extract_variance(fit)
vec_gen = c(data_gen[1],
            data_gen[4],
            data_gen[5],
            data_gen[4],
            data_gen[2],
            data_gen[6],
            data_gen[5],
            data_gen[6],
            data_gen[3])
est_mat_gen = matrix(vec_gen, nrow = 3)
save(est_mat_gen, file = "est_mat_gen.Rdata")

# env
data_env = extract_variance(fit)
vec_env = c(data_env[1+6],
            data_env[4+6],
            data_env[5+6],
            data_env[4+6],
            data_env[2+6],
            data_env[6+6],
            data_env[5+6],
            data_env[6+6],
            data_env[3+6])
real_mat_env = matrix(vec_env, nrow = 3)

# get the matrix A 
position = unique(df$id)[order(unique(df$id))]
A = getA(pedigree)[position,position]

# big genetic covariance  
c_gen = as(kronecker(A, est_mat_gen),"dgCMatrix")

# big env covariance
c_env = as(kronecker(diag(length(unique(df$id))), real_mat_env),"dgCMatrix")

# general covariance
c = as(bdiag(c_gen,c_env),"dgCMatrix")

# length of the vector u
ul = n_of_terms_gen*length(unique(df$id)) + n_of_terms_nongen*length(unique(df$id)) 
iterations = iterations
n = dim(df)[1]  

sds = matrix(0, nrow = 0.5*n_of_terms_gen*(n_of_terms_gen+1)+0.5*(n_of_terms_nongen)*(n_of_terms_nongen+1)+1, ncol = iterations)

g = lm(trait ~ -1+df$z1+df$z2+df$z3, data = df)

n = dim(df)[1]
length_error_1 = n
mu_u = rep(0, ul)
mu_e = rep(0, length_error_1)
c_e = as(s_2_residuals*diag(length_error_1), "dgCMatrix")

f = df$z1*g$coefficients[1] + df$z2*g$coefficients[2] + df$z3*g$coefficients[3]

df$id2 = df$id

# progress bar
total <- iterations
pb <- txtProgressBar(min = 0, max = total, style = 3)

system.time(
  for( i in 1:iterations)
  {
    # generate the random effect vector composed of alpha and gamma
    u = as.vector(rmvn(1, mu_u, c, ncores = 2))
    
    # the error vector
    e = as.vector(rmvn(1, mu_e, c_e, ncores = 2))
    
    form = lFormula(df$trait ~ -1 + df$z1 +df$z2+df$z3+(-1 + df$z1 + df$z2 +df$z3|df$id) + (-1 + df$z1 + df$z2 +df$z3|df$id), 
                    data = df)
    
    # 1st version: in which i have as many generated points as in the dataset
    df$trait = as.vector(f + t(form$reTrms$Zt)%*%u + e)
    
    # # store it to analyze it later
    # responses = cbind(responses, df$trait)
    
    # now we fit the model, the result should give out the same cov matrix as before
    fm1 <- fit_me( df$trait ~ -1 + df$z1 +df$z2+df$z3+(-1 + df$z1 + df$z2 +df$z3|df$id) + (-1 + df$z1 + df$z2 +df$z3|df$id),
                   data = df,
                   label = id, 
                   pedigree = pedigree, 
                   n = n_of_terms_gen,
                   optimizer = 'nloptwrap')# nloptwrap, Nelder_mead, bobyqa
    
    samp = as.matrix(as.data.frame(VarCorr(fm1))["vcov"])
    sds[,i] = samp
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
  }
)

sample_27 = sds 
save(sample_27, file = "sample_27.Rdata")

# Pointwise confidence bands for covariance function and eigenvectors
fine_grid = 25

# covariance function
# real covariance function
real_gen_cov_func = creator_gen_cov_function(small_mat_gen)
t1_gen <- t2_gen <- seq(-1, 1, length= fine_grid)
real_gen_values <- outer(t1_gen, t2_gen, real_gen_cov_func)

# estimated genetic covariance function
est_gen_cov_fun = creator_gen_cov_function(est_mat_gen)
est_gen_values <- outer(t1_gen, t2_gen, est_gen_cov_fun)

persp(t1_gen, 
      t2_gen, 
      est_gen_values, 
      theta = 30,
      phi = 15,
      shade = 0.45,
      col = "aquamarine",
      ticktype = "detailed",
      nticks = 4)

# compute the lower and upper bound for each point

# ar will be the 3D matrix with all the values of the functions in the grid
zeros = rep(0,fine_grid*fine_grid*iterations)
ar_values_cov <- array(zeros, c(fine_grid, fine_grid, iterations))
ar_differences_cov <- array(zeros, c(fine_grid, fine_grid, iterations))
ar_values_eig <- array(zeros, c(fine_grid, iterations))
ar_differences_eig <- array(zeros, c(fine_grid, iterations))

# compute the function
for ( i in 1:iterations)
{
  data_gen = sds[,i]
  vec_gen = c(data_gen[1],
              data_gen[4],
              data_gen[5],
              data_gen[4],
              data_gen[2],
              data_gen[6],
              data_gen[5],
              data_gen[6],
              data_gen[3])
  mat_gen = matrix(vec_gen, nrow = 3)
  
  gen_cov_func = creator_gen_cov_function(mat_gen)
  
  t1 <- t2 <- seq(-1, 1, length = fine_grid)
  ar_values_cov[,,i] =  outer(t1, t2, gen_cov_func)
  ar_differences_cov[,,i] = ar_values_cov[,,i] - est_gen_values
  
  # eigenfunctions
  est_eigenfunction = eigen(est_gen_values)$vectors[,1]
  
  ar_values_eig[, i] = eigen(ar_values_cov[,,i])$vectors[,1]
  ar_differences_eig[, i] = ar_values_eig[, i]-est_eigenfunction
}
contenuto_cov = rep(0,fine_grid*fine_grid*2)
contenuto_eig = rep(0,fine_grid*2)
ul_bound_cov <- array(contenuto_cov, c(fine_grid, fine_grid, 2))
ul_bound_eig <- array(contenuto_eig, c(fine_grid, 2))

for ( i in 1:fine_grid)
{
  for (j in 1:fine_grid)
  {
    upper_lower_cov = quantile(ar_differences_cov[i, j,], 
                               probs = c(0.025, 0.975))
    ul_bound_cov[i, j,] = c(est_gen_values[i, j]-upper_lower_cov[2],
                            est_gen_values[i, j]-upper_lower_cov[1])
  }
  upper_lower_eig = quantile(ar_differences_eig[i,], 
                            probs = c(0.025, 0.975))
  ul_bound_eig[i,] = c(est_eigenfunction[i] -upper_lower_eig[2],
                       est_eigenfunction[i] -upper_lower_eig[1])
}

# plot covariance function
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~real_gen_values)
fig <- fig %>% add_surface(z = ~ul_bound_cov[,,1], opacity = 0.98, colorscale = list(c(0, 0), c("tan", "blue")))
fig <- fig %>% add_surface(z = ~ul_bound_cov[,,2], opacity = 0.98, colorscale = list(c(0, 0), c("tan", "blue")))
fig

real_eigenfunction = eigen(real_gen_values)$vectors[,1]
plot(1:20, -real_eigenfunction, type = 'l', ylim = c(-0.5,.2))
points(1:20, est_eigenfunction, col = "green", type = 'l')
points(1:20, ul_bound_eig[,1], col = "red", type = 'l')
points(1:20, ul_bound_eig[,2], col = "red", type = 'l')

# simultaneous confidence interval
vec_max_diff_cov = rep(0, iterations)
vec_max_diff_eig = rep(0, iterations)

for(k in 1:iterations)
{
  vec_max_diff_cov[k] = max(abs(ar_differences_cov[,,k]))
  vec_max_diff_eig[k] = max(abs(ar_differences_eig[,k]))
}

quantile_cov = quantile(vec_max_diff_cov, probs = 0.975)
quantile_eig_sim = quantile(vec_max_diff_eig, probs = 0.975)

ul_bound_cov_sim = array(contenuto_cov, c(fine_grid, fine_grid, 2))
ul_bound_eig_sim = array(contenuto_eig, c(fine_grid, 2))
  
ul_bound_cov_sim[,,1] = est_gen_values-quantile_cov
ul_bound_cov_sim[,,2] = est_gen_values+quantile_cov

ul_bound_eig_sim[, 1] = est_eigenfunction-quantile_eig_sim
ul_bound_eig_sim[, 2] = est_eigenfunction+quantile_eig_sim

# plotting
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~real_gen_values)
fig <- fig %>% add_surface(z = ~ul_bound_cov_sim[,,1], opacity = 0.98, colorscale = list(c(0, 0), c("tan", "blue")))
fig <- fig %>% add_surface(z = ~ul_bound_cov_sim[,,2], opacity = 0.98, colorscale = list(c(0, 0), c("tan", "blue")))
fig

real_eigenfunction = eigen(real_gen_values)$vectors[,1]
plot(1:20,- real_eigenfunction, type = 'l', ylim = c(-0.5,.5))
points(1:20, est_eigenfunction, col = "green", type = 'l')
points(1:20, ul_bound_eig_sim[,1], col = "red", type = 'l')
points(1:20, ul_bound_eig_sim[,2], col = "red", type = 'l')