library(matlib)
n = 100
d = 200
M = 500
theta_0 = rep(1/sqrt(d),d)
I <- diag(d)

arr <- 1:d
diff <- abs(outer(arr, arr, "-"))

q = 0.2
Sigma_matrix <- 0.1 ^ diff

Err_1 = rep(0, M)
Err_2 = rep(0, M)
Err_3 = rep(0, M)
Err_4 = rep(0, M)
Err = rep(0, M)
starttime = Sys.time()
for (i in 1:M){
  set.seed(2000+i)
  if (i%%100 == 0){
    print(i)
    print(Sys.time() - starttime)
    starttime = Sys.time()
  }
  Z_1 = matrix(rnorm( n * d ), n, d) %*% Sigma_matrix
  Z_2 = matrix(rnorm( n * d ), n, d) %*% Sigma_matrix
  
  X_1 = Z_1 %*% theta_0 + rnorm(n)
  X_2 = Z_2 %*% theta_0 + rnorm(n)
  Y = Z_2%*% theta_0 + rnorm(n)

  theta_hat = t(Z_1) %*% solve(Z_1 %*% t(Z_1)) %*% X_1 * d/n
  err_1 = X_2 - Z_2 %*% theta_hat
  err_2 = sum(err_1 * Y)
  Err_1[i] = sum((theta_0 - theta_hat)^2)
  Err_2[i] = sum((theta_0 - theta_hat)^2)
  Err_3[i] = sum(err_1^2)
  Err_4[i] = err_2
}
summary(Err_1)
summary(Err_2)
summary(Err_3)
summary(Err_4)
