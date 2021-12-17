

eigen_fun = function(X,tol=.Machine$double.eps^0.75){
  if(!isSymmetric(X)) stop("x is not symmetric")
  d = ncol(X)
  U = diag(1,d,d,FALSE)
  while(TRUE){
    did_we_update = FALSE
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        # test whether abs(X[i,j]) is too large
        if(abs(X[i,j])>tol){
          # if so, get the Givens Rotation matrix R
          # and update X and U
          R=diag(1,d,d,FALSE)
          if(X[i,i]==X[j,j]){
            sin=sin(pi/4)
            cos=cos(pi/4)
            R[i,i]=sin
            R[i,j]=cos
            R[j,j]=sin
            R[j,i]=-1*cos
            }else{
            theta=1/2*atan(2*X[j,i]/(X[i,i]-X[j,j]))
            sin=sin(theta)
            cos=cos(theta)
            R[i,i]=sin
            R[i,j]=cos
            R[j,j]=sin
            R[j,i]=-1*cos
          }
            
           X=R%*%X%*%t(R)
           U=U%*%R
            
          
          
          # track that in this sweep through we have done at least one update
          did_we_update = TRUE
        }
        
      }
    }
    if(!did_we_update) break
  }
  values=diag(X)
  sorted_values = sort(values,index.return=TRUE,decreasing = True)
  values = sorted_values$X
  vectors = U[,sorted_values$iX]
  return(list(values=values, vectors=vectors))
}
d=2
z = matrix(rnorm(d^d),d,d)
X = z + t(z)
set.seed(1234567890)
eigen_fun(X)

d=3
z = matrix(rnorm(d^d),d,d)
X = z + t(z)
set.seed(1234567890)
eigen_fun(X)

d=4
z = matrix(rnorm(d^d),d,d)
X = z + t(z)
set.seed(1234567890)
eigen_fun(X)

d=5
z = matrix(rnorm(d^d),d,d)
X = z + t(z)
set.seed(1234567890)
eigen_fun(X)

