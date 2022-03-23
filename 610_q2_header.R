library(mvtnorm)
library(SuppDists)
load("final_q2_data.RData")
x = as.matrix(x)
y = as.matrix(y)

draw_tau = function(data,state){
  shape = 1 + length(data$y)
  rate = 1 + sum(abs(data$y-data$x%*%state$beta))
  return(rgamma(1,shape,rate))
}

draw_beta = function(data,state){
  n = length(y)
  m = ncol(x)# m = 10
  means_z = 1/(state$tau*abs(data$y-data$x%*%state$beta)) 
  shape_z = 1 
  z <<- rinvGauss(n,means_z,shape_z) 
  
  Z = diag(z)
  a_ = solve(state$tau^2*t(data$x)%*%Z%*%(data$x)+diag(m)) 
  mean_beta =a_%*%(state$tau^2*t(data$x)%*%Z%*%data$y) 
  return(rmvnorm(m,mean_beta,diag(10)))
  
  
}
update_current_state = function(data, state){
  state$tau = draw_tau(data,state) 
  draw.beta = draw_beta(data,state) 
  state$beta = draw.beta[,sample(1:10,1)] 
  return(state)
}
