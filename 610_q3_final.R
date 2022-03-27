

eigen_fun <- function(X,tol=.Machine$double.eps^0.75){
  d = ncol(X)
  U = diag(1,d,d,FALSE)
  
  Rlocal <- function(A){
    d = ncol(X)
    A <- C
  B <- diag(1,d,d,FALSE)
    for (i in 1:(d-1)) {
      X <- matrix(rep(0,d),nrow = d)
      X[i:d] <- R2[i:d,i]
      V <- X
      V[i] <- X[i]+norm(X,"F")*sign(X[i])
      xlocal <- V/norm(V,"F")
      A <- A-2*(xlocal%*%t(xlocal)%*%A)
     B <- B-2*(xlocal%*%t(xlocal)%*%B)
    }
    return(list("B" = t(B), "C" = C))
  }
  
  while(TRUE){
    did_we_update = FALSE
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        # test whether abs(X[i,j]) is too large
        if(abs(X[i,j])>tol){
          # if so, get the Givens Rotation matrix R
          # and update X and U
          for (k in 1:1000){
          update <- Rlocal(U)
          Q = update$B
          R = update$A
          U = X %*% Q
          }
          # track that in this sweep through we have done at least one update
          did_we_update = FALSE
        }
        
      }
    }
    if(!did_we_update) break
  }
  
  values = diag(R)
  sorted_values = sort(values,index.return=TRUE,decreasing=TRUE)
  values = sorted_values$x
  vectors = Q[,sorted_values$ix]
  return(list(values=values, vectors= -vectors))
  
}
set.seed(1234567890)
d =3
z = matrix(rnorm(d^d),d,d)
X = z + t(z)
# X
eigen_fun(X)
eigen(X)

test_that("this is d =3 value test",{
  expect_equal(eigen_fun(X)$values,eigen(X)$values)
}
)


set.seed(1234567890)
d =2
z = matrix(rnorm(d^d),d,d)
X = z + t(z)
# X
eigen_fun(X)
eigen(X)

test_that("this is d =2 value test",{
  expect_equal(eigen_fun(X)$values,eigen(X)$values)
}
)

set.seed(1234567890)
d =4
z = matrix(rnorm(d^d),d,d)
X = z + t(z)
# X
eigen_fun(X)
eigen(X)

test_that("this is d =4 value test",{
  expect_equal(eigen_fun(X)$values,eigen(X)$values)
}
)

set.seed(1234567890)
d =5
z = matrix(rnorm(d^d),d,d)
X = z + t(z)
# X
eigen_fun(X)
eigen(X)

test_that("this is d =5 value test",{
  expect_equal(eigen_fun(X)$values,eigen(X)$values)
}
)
