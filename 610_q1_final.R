
set.seed(1234567890)
Rprof(interval = 1e-4, filter.callframe=TRUE)
# direct inversion
r_squared_fun_1 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  xty = crossprod(x,y)
  xtx = crossprod(x)
  beta = solve(xtx)%*%xty
  SST = crossprod(xty,beta)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# solve linear system
r_squared_fun_2 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  xty = crossprod(x,y)
  xtx = crossprod(x)
  beta = solve(xtx,xty)
  SST = crossprod(xty,beta)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
  
}

# invert using eigen system
r_squared_fun_3 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  xty = crossprod(x,y)
  xtx = crossprod(x)
  eigen_xtx = eigen(xtx,symmetric=TRUE)
  u = eigen_xtx$vectors
  d_inv = diag(1/eigen_xtx$values)
  beta = u %*% d_inv %*% t(u) %*% xty
  SST = crossprod(xty,beta)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# avoid inversion using svd
r_squared_fun_4 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  svd_x = svd(x)
  u = svd_x$u
  uty = crossprod(u,y)
  SST = crossprod(uty)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}

# ignore right singular vectors
r_squared_fun_5 = function(x,y){
  y = y - mean(y)
  x = sweep(x,2,colMeans(x),"-")
  SSE_intercept_only = sum(y^2)
  svd_x = svd(x,nv=0)
  u = svd_x$u
  uty = crossprod(u,y)
  SST = crossprod(uty)[1,1]
  r_squared = SST/SSE_intercept_only
  return(r_squared)
}
output = matrix(NA,5,20)

r_squared_fun_i = function(i,x,y){
  if (i == 1){
    r_squared_fun_1(x,y)
  }
  if (i ==2){
    r_squared_fun_2(x,y)
  }
  if (i ==3){
    r_squared_fun_3(x,y)
  }
  if (i ==4){
    r_squared_fun_4(x,y)
  }
  if (i ==5){
    r_squared_fun_5(x,y)
  }
}

totaltime = function(i,x,y){
  Rprof(filename='q.out',filter.callframes = TRUE,interval=0.01)
  r_squared_fun_i(i,x,y)
  Rprof(NULL)
  summaryRprof('q.out')$by.total[1,1]
}

for (j in 1:20){
  x = matrix(rnorm(1000*500),1000,500)
  y = rnorm(1000)
  for (i in 1:5){
    output[i,j] = totaltime(i,x,y)
  }
}
output

#function 2 is the best among all
#because SVD function is time consuming.
profiling_matrix <- output
f1 = summary(profiling_matrix[1,])
f2 = summary(profiling_matrix[2,])
f3 = summary(profiling_matrix[3,])
f4 = summary(profiling_matrix[4,])
f5 = summary(profiling_matrix[5,])
f = rbind(f1,f2,f3,f4,f5)
summary_matrix = as.matrix(f)
summary_matrix


