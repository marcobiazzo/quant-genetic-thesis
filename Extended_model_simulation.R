library(MASS)
library(ggplot2)
#library(ggpubr)

data(milk)

# Create a smaller milk dataframe to let the fitting be faster (200 samples)
milk <- within(milk, sdMilk <- milk / sd(milk))
milk <- milk[sample(nrow(milk)),]
m = data.frame(milk[1:300,c(10,1,2,5,3)]) # we eliminate the columns we don't need, at first we fit a model just for the random intercept
colnames(m)= c("sdmi","ide","l","d","h")
head(m)

# We standardize the dataset m
m <- within(m, d <- d / sd(d))
head(m)

# To try the extended model since we have identifiability issues we need to duplicate some rows
m = rbind(m,m)

attach(m)

n_of_terms_gen = 2
n_of_terms_nongen = 2

# Get the matrix A and then multiplying by the variance get the cov matrix of the random effect
pos = unique(ide)[order(unique(ide))]
A = getA(pedCowsR)[pos,pos]
dim(A)

real_model = TRUE

# Fix a value for the multiplication constant of the matrix A to be the covariance matrix
s_2_reffect = 2
if(real_model)
{
  small_func_mat_reffect = 0.1*diag(n_of_terms_gen)
  small_func_mat_reffect[1,1] = 1
  small_func_mat_reffect = s_2_reffect*small_func_mat_reffect
}else{
  small_func_mat_reffect = s_2_reffect*diag(n_of_terms_gen)
}

small_func_mat_reffect
c_gen = kronecker(A, small_func_mat_reffect)# covariance matrix of the random effect vector


# Environmental covariance matrix

env = TRUE
small_func_mat_env = NULL

if (env)
{
  s_2_env = 2
  if(real_model)
  {
    small_func_mat_env = 0.1*diag(n_of_terms_nongen)
    small_func_mat_env[1,1] = 1
    small_func_mat_env = s_2_env*small_func_mat_env
  } else {
    small_func_mat_env = s_2_env*diag(n_of_terms_nongen)
  }
  
  small_func_mat_env
  c_env = kronecker(diag(length(unique(h))), small_func_mat_env)
  
  # General covariance matrix of the big u
  c = bdiag(c_gen,c_env)
} else {
  c = c_gen
}
dim(c)

# Let's fix a value for the variance of the error sigma^2
s_2_residuals = 0.01

ul = n_of_terms_gen*length(unique(ide)) + n_of_terms_nongen*length(unique(h)) 

iterations = 10

n = dim(m)[1]

sds = matrix(0, nrow = 5, ncol = iterations)

g = lmer(sdmi ~ l + log(d) + (1|ide) + (1|h) ,data = m)

f = rep(1,n)*g@beta[1] + m$l*g@beta[2] + log(m$d)*g@beta[3]

total <- iterations
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)

system.time(
  for( i in 1:iterations)
  {
    u = mvrnorm(n = 1, rep(0,ul), c)
    e = mvrnorm(n = 1, rep(0, n), s_2_residuals*diag(n))
    
    form = lFormula(sdmi ~ l + log(d) + (1 + log(d)|ide) + (1 + log(d)|h), data = m)
    
    # We create a synthetic response
    m$sdmi = as.vector(f + t(form$reTrms$Zt)%*%u + e)
    fm1 <- fit_me( form,
                   data = m,
                   label = ide, 
                   pedigree = pedCowsR, 
                   n = n_of_terms_gen, 
                   optimizer = "Nelder_Mead")
    
    # Fitting pedigreemm to compare the distributions
    # g = pedigreemm(sdmi ~ l + log(d) + (1|ide) ,data = m, pedigree = list(ide = pedCowsR))
    
    samp = as.matrix(as.data.frame(VarCorr(fm1))["vcov"])
    # if(samp[3,1]<2.26e-06 || samp[4,1]<2.26e-06)
    # {
    #   ss = getME(fm1,"theta")
    #   fm1 <- lmer(sdmi ~ l + log(d) + (1 + log(d)||ide)+(1+log(d)||h), m, start = ss)
    # }
    
    sds[,i] = samp
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, i)
  }
)

par(mfrow=c(3,3))
hist(sds[1,], breaks = 40)
hist(sds[2,], breaks = 40)
hist(sds[3,], breaks = 40)
hist(sds[4,], breaks = 40)
hist(sds[5,], breaks = 40)
hist(sds[6,], breaks = 40)
hist(sds[7,], breaks = 40)

sds_500_Nelder_Mead_800rowsdataset = sds

save(sds_500_Nelder_Mead_800rowsdataset, file="sds_500_Nelder_Mead_800rowsdataset.Rdata")

result = data.frame(var_comp_1 = sds[1,], var_comp_2 = sds[2,], var_res = sds[3,])

p_var_comp_1<- ggplot(result, aes(x= var_comp_1)) +
  geom_histogram(color="black", fill="aquamarine3", binwidth=0.1)+
  geom_vline(aes(xintercept=2), color="black", linetype="dashed", size=1)+
  labs(title="Residuals", x ="Fitted Residual Variance", y = "Count")+
  theme(plot.title = element_text(hjust = 0.5))
p_var_comp_1

p_var_comp_2<- ggplot(result, aes(x= var_comp_2)) +
  geom_histogram(color="black", fill="aquamarine3", binwidth=0.1)+
  geom_vline(aes(xintercept=2), color="black", linetype="dashed", size=1)+
  labs(title="Residuals", x ="Fitted Residual Variance", y = "Count")+
  theme(plot.title = element_text(hjust = 0.5))
p_var_comp_2

p_var_res<- ggplot(result, aes(x= var_res)) +
  geom_histogram(color="black", fill="salmon", binwidth=0.0015)+
  geom_vline(aes(xintercept=0.1), color="black", linetype="dashed", size=1)+
  labs(title="Residuals", x ="Fitted Residual Variance", y = "Count")+
  theme(plot.title = element_text(hjust = 0.5))
p_var_res

ggarrange(p_var_comp_1, p_var_comp_2, p_var_res, nrow = 2, ncol = 2)

