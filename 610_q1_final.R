
set.seed(1234567890)
x1<-rnorm(500000)
x <- matrix(x1, nrow=1000, ncol=500, byrow=TRUE)
x

y<-rnorm(1000)
length(y)
y




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
Rprof(NULL)
Rprof_summ1= summaryRprof()
Rprof_summ1
tt1<-Rprof_summ1$by.total$total.time

Rprof(interval = 1e-4, filter.callframe=TRUE)
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
Rprof(NULL)
Rprof_summ2 = summaryRprof()
Rprof_summ2
tt2<-Rprof_summ2$by.total$total.time

Rprof(interval = 1e-4, filter.callframe=TRUE)
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
Rprof(NULL)
Rprof_summ3 = summaryRprof()
Rprof_summ3
tt3<-Rprof_summ3$by.total$total.time

Rprof(interval = 1e-4, filter.callframe=TRUE)
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
Rprof(NULL)
Rprof_summ4 = summaryRprof()
Rprof_summ4
Rprof_summ4 = summaryRprof()
Rprof_summ4
tt4<-Rprof_summ4$by.total$total.time


Rprof(interval = 1e-4, filter.callframe=TRUE)
# avoid inversion using svd
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
Rprof(NULL)
Rprof_summ5 = summaryRprof()
Rprof_summ5
tt5<-Rprof_summ5$by.total$total.time

f<-cbind(tt1,tt2,tt3,tt4,tt5)

Profiling_matrix <- as.data.frame(t(f))
Profiling_matrix

rownames(Profiling_matrix)<-c("r_squared_fun_1","r_squared_fun_2","r_squared_fun_3","r_squared_fun_4","r_squared_fun_5")
print(Profiling_matrix)


#20 columns cannot be created as there are only 11 features when the profiling function is run r squared function.
# therefore we just obtain 5 by 11 matrix.
# the row names are given as the function names and the column indicates the total time spent on each function.

summary(Profiling_matrix)

final_matrix<- Profiling_matrix[,-c(1:4,11)]
final_matrix

rowSums(final_matrix)

#the r_squared_fun_2 takes the least time to compute when compared with the other function.
#beacuse it runs in very less time.
#also t(x)%*%x was calculated in function1 itself , so the remaining part was just very simple and small.
#therefore the computation took less time.


