draw_tau = function(data,state){
  shape = 1 + length(data$y)
  rate = 1 + sum(abs(data$y-data$x%*%state$beta))
  return(rgamma(1,shape,rate))
}

update_current_state = function(data, state){
  state$beta = draw_beta(data,state)
  state$tau = draw_tau(data,state)
  return(state)
}