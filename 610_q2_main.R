# sources header relative to the path for this file
source("../headers/header.R")

# gibbs sampling function
gibbs_laplace_regression = function(N_draws,x,y,beta_0,tau_0){
  # data manipulation
  if(is.null(colnames(x))){
    colnames(x) = paste("x",1:ncol(x),sep="_")
  }
  if(!all(abs(x[,1]-1)<.Machine$double.eps)){
    x = cbind(c(1),x)
    colnames(x)[1] = "intercept"
  }
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

