# sources header relative to the path for this file
source("../headers/header.R")

# gibbs sampling function
gibbs_laplace_regression = function(N_draws,x,y,beta_0,tau_0){
  # data manipulation
 
  data = list(y=y,x=x)
  
  # checking for missing initial values and making them if necessary
  if(missing(beta_0)) beta_0 = c(median(data$y),rep(0,ncol(data$x)-1))
  if(missing(tau_0)) tau_0 = 1/mean(abs(data$y-data$x%*%beta_0))
  
  # initial state creation
  current_state = list(beta=beta_0,tau=tau_0)
  
  # output creation
  out = list(beta=matrix(NA,ncol(x),N_draws),tau=rep(NA,N_draws))
  rownames(out$beta) = colnames(x)
  
  # loop for sampling
  for(i in 1:N_draws){
    current_state = update_current_state(date,current_state)
    out$beta[,i] = current_state$beta
    out$tau[i] = current_state$tau
  }
  
  # returning
  return(out)
}
# test 100
# gibbs_laplace_regression(100,x,y)

#get 10000 draws
result = gibbs_laplace_regression(10000,x,y)
# result$beta
# result$tau
apply(result$beta,1,mean)# posterior means of the beta draws
summary(result$tau)[4] # posterior means of the tau draws

##########################################
#make univariate 95 intervals for each element of beta#
mean_ = apply(result$beta,1,mean)
sd = apply(result$beta,1,sd)

lower = mean_ - 1.96*sd/100
upper = mean_ + 1.96*sd/100
ConfidenceInterval = rbind(lower,upper)
ConfidenceInterval

# I think the column 1,2,5,6,7,8 are non-zero in this Laplace regression.
