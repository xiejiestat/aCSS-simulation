library(glmnet)
library(gglasso)
library(abess)
library(MASS)
library(ncvreg)

set.seed(123)
n = 50 
d = 100 
tol = 1e-6
M = 200 # copies number
N = 5000 # replicate times

theta_0 = c(rep(5, 5), rep(0, d - 5)) 
supp = 5
xi_0 = c(rep(1, 5), rep(0, d-5))
group = rep(1:(d/5), each = 5)
sigma_all = seq(from = 1, to = 8, by = 1)
nsigma =length(sigma_all)
signal = seq(from = 0, to = 1, by = 0.2)
nsignal = length(signal)

SSOSP = function(theta, Z){
  index = which(abs(theta) > tol) 
  submatrix <- (t(Z) %*% Z)[index, index, drop = FALSE]
  if (nrow(submatrix) != ncol(submatrix)) {
    stop("The submatrix is not square. Please provide a valid column index vector.")
  }
  is_positive_definite <- all(eigen(submatrix)$values > 0)
  return(is_positive_definite)
}

estimate_lasso = function(X_ran, Z){
  lasso_model = cv.glmnet(Z, X_ran, alpha = 1, intercept = FALSE)
  best_lambda_lasso = lasso_model$lambda.min
  lasso_best_model = glmnet(Z, X_ran, alpha = 1, lambda = best_lambda_lasso, intercept = FALSE)
  lasso_coef = coef(lasso_best_model)
  theta = as.vector(lasso_coef)[-1]
  
  theta
}

estimate_group = function(X_ran, Z){
  group_lasso_model = cv.gglasso(Z, X_ran, group = group, intercept = FALSE)
  best_lambda_group = group_lasso_model$lambda.min
  group_lasso_best_model = gglasso(Z, X_ran, group = group, lambda = best_lambda_group, intercept = FALSE)
  group_lasso_coef = coef(group_lasso_best_model)
  theta = as.vector(group_lasso_coef)[-1]
  
  theta
}

estimate_subset = function(X_ran, Z, supp){
  abess_fit = abess(Z, X_ran,fit.intercept = FALSE)
  print(abess_fit)
  abess_coef= coef(abess_fit, support.size = supp)
  theta = as.vector(abess_coef)[-1]
  
  theta
}

estimate_ridgeless = function(X_ran, Z) {
  ZtZ_inv = ginv(t(Z) %*% Z)  # Moore-Penrose pseudo-inverse
  theta = as.vector(ZtZ_inv %*% t(Z) %*% X_ran)
  
  theta
}

estimate_SCAD = function(X_ran, Z){
  scad_fit <- ncvreg(Z, X_ran, penalty = "SCAD")
  cv_scad <- cv.ncvreg(Z, X_ran, penalty = "SCAD")
  best_lambda_scad <- cv_scad$lambda.min
  scad_coef <- coef(scad_fit, lambda = best_lambda_scad)
  theta = as.vector(scad_coef)[-1]
  
  theta
}

estimate_MCP = function(X_ran,Z){
  mcp_fit <- ncvreg(Z, X_ran, penalty = "MCP")
  cv_mcp <- cv.ncvreg(Z, X_ran, penalty = "MCP")
  best_lambda_mcp <- cv_mcp$lambda.min
  mcp_coef <- coef(mcp_fit, lambda = best_lambda_mcp)
  theta = as.vector(mcp_coef)[-1]
  
  theta
}

lasso_test = function(X ,Y, Z){
  penalty_factor = rep(0, d+1)
  penalty_factor[2:d+1] = 1
  X_full = cbind(X, Z)
  lasso_model = cv.glmnet(X_full, Y, alpha = 1, intercept = FALSE, penalty_factor = penalty_factor)
  best_lambda_lasso = lasso_model$lambda.min
  lasso_best_model = glmnet(X_full, Y, alpha = 1, lambda = best_lambda_lasso, penalty_factor = penalty_factor, intercept = FALSE)
  abs(coef(lasso_best_model)[2])
}

pval_oracle_compute = function(X, Y, Z, theta, M){
  pval = rep (0,3)
  var = 1
  mu = Z%*%theta
  
  # marginal covariance
  delta_1 = (t(Y)%*%(mu-X))
  var_1 = sqrt(var*sum(Y^2))
  pval[1] = pnorm(delta_1/var_1)
  
  # lasso based dCRT
  xi_hat = estimate_lasso(Y, Z)
  Y_hat = Y-Z%*%xi_hat
  delta_2 = (t(Y_hat)%*%(mu-X))
  var_2 = sqrt(var*sum(Y_hat^2))
  pval[2] = pnorm(delta_2/var_2)
  
  # lasso estimate
  Xcopies = list()
  for (m in 1:M){
    Xcopies[[m]] = c(mu) + c(sqrt(var)*rnorm(n))
  }
  lasso_ = function(x){lasso_test(x,Y,Z)}
  pval[3] = 
    (1+sum(unlist(lapply(Xcopies,lasso_))>=lasso_(X)))/(1+M)
  
  pval
}

pval_compute = function(X, Y, Z, W, theta, sigma, M){ #compute p-value
  pval = rep (0,3)
  var = sigma^2/(1+sigma^2)
  mu = var*(Z%*%theta) + (X+sigma*W)/(1+sigma^2)
  
  # marginal covariance
  delta_1 = (t(Y)%*%(mu-X))
  var_1 = sqrt(var*sum(Y^2))
  pval[1] = pnorm(delta_1/var_1)
  
  # lasso based dCRT
  xi_hat = estimate_lasso(Y, Z)
  Y_hat = Y-Z%*%xi_hat
  delta_2 = (t(Y_hat)%*%(mu-X))
  var_2 = sqrt(var*sum(Y_hat^2))
  pval[2] = pnorm(delta_2/var_2)
  
  # lasso estimate
  Xcopies = list()
  for (m in 1:M){
    Xcopies[[m]] = c(mu) + c(sqrt(var)*rnorm(n))
  }
  lasso_ = function(x){lasso_test(x,Y,Z)}
  pval[3] = 
    (1+sum(unlist(lapply(Xcopies,lasso_))>=lasso_(X)))/(1+M)
  
  pval
}

Pval_oracle = vector("list", N * nsignal)
dim(Pval_oracle) = c(N, nsignal)

Pval_lasso = Pval_group = Pval_ridgeless = Pval_subset = Pval_SCAD = Pval_MCP = vector("list", N * nsignal * nsigma) 
dim(Pval_lasso) = dim(Pval_group) = dim(Pval_ridgeless) = dim(Pval_subset) = dim(Pval_SCAD) = dim(Pval_MCP) = c(N, nsignal, nsigma)

for (i in 1:N){
  print(i)
  for (k in 1:nsigma){
    #simulate data
    Z = matrix(rnorm(n * d), n, d) 
    X = Z %*% theta_0 + rnorm(n)
    
    W = rnorm(n)
    sigma = sigma_all[k]/sqrt(d)
    X_ran = X + sigma*W
    
    # estimation
    theta_hat_lasso = estimate_lasso(X_ran, Z)
    theta_hat_group = estimate_group(X_ran, Z)
    theta_hat_ridgeless = estimate_ridgeless(X_ran, Z)
    theta_hat_subset = estimate_subset(X_ran, Z, supp)
    theta_hat_SCAD = estimate_SCAD(X_ran, Z)
    theta_hat_MCP = estimate_MCP(X_ran, Z)
    theta_oracle = theta_0
    
    # different magnitude
    for (j in 1:nsignal){
      Y = signal[j]*X + c(Z %*% xi_0) + rnorm(n)
      
      # p-value
      Pval_oracle[[i,j]] = pval_oracle_compute(X, Y, Z,theta_oracle, M)
      Pval_lasso[[i,j,k]] = pval_compute(X, Y, Z, W, theta_hat_lasso, sigma, M)
      Pval_group[[i,j,k]] = pval_compute(X, Y, Z, W, theta_hat_group, sigma, M)
      Pval_ridgeless[[i,j,k]] = pval_compute(X, Y, Z, W, theta_hat_ridgeless, sigma, M)
      Pval_subset[[i,j,k]] = pval_compute(X, Y, Z, W, theta_hat_subset, sigma, M)
      Pval_SCAD[[i,j,k]] = pval_compute(X, Y, Z, W, theta_hat_SCAD, sigma, M)
      Pval_MCP[[i,j,k]] = pval_compute(X, Y, Z, W, theta_hat_MCP, sigma, M)
    }
  }
}


#power plot
alpha = 0.1

power_oracle  = vector("list", nsignal * 3)
dim(power_oracle) = c(nsignal, 3)

power_lasso = power_group = power_ridgeless = power_subset = power_SCAD = power_MCP = vector("list", nsignal * nsigma * 3) 
dim(power_lasso) = dim(power_group) = dim(power_ridgeless) = dim(power_subset) = dim(power_SCAD) = dim(power_MCP) = c(nsignal, nsigma, 3)

for (j in 1:nsignal){
  for(k in 1:nsigma){
    for (t in 1:3){
      # initial value
      power_oracle[[j, t]] = 0
      power_lasso[[j, k, t]] = 0
      power_group[[j, k, t]] = 0
      power_ridgeless[[j, k, t]] = 0
      power_subset[[j, k, t]] = 0
      power_SCAD[[j, k, t]] = 0
      power_MCP[[j, k, t]] = 0
      
      # power summation
      for (i in 1:N){
        power_oracle[[j, t]] = power_oracle[[j, t]] + as.numeric(Pval_oracle[[i,j]][t] < alpha)
        power_lasso[[j, k, t]] = power_lasso[[j, k, t]] + as.numeric(Pval_lasso[[i, j, k]][t] < alpha) 
        power_group[[j, k, t]] = power_group[[j, k, t]] + as.numeric(Pval_group[[i, j, k]][t] < alpha) 
        power_ridgeless[[j, k, t]] = power_ridgeless[[j, k, t]] + as.numeric(Pval_ridgeless[[i, j, k]][t] < alpha) 
        power_subset[[j, k, t]] = power_subset[[j, k, t]] + as.numeric(Pval_subset[[i, j, k]][t] < alpha) 
        power_SCAD[[j, k, t]] = power_SCAD[[j, k, t]] + as.numeric(Pval_SCAD[[i, j, k]][t] < alpha) 
        power_MCP[[j, k, t]] = power_MCP[[j, k, t]] + as.numeric(Pval_MCP[[i, j, k]][t] < alpha) 
      }
      
      # power average
      power_oracle[[j, t]] = power_oracle[[j, t]]/N
      power_lasso[[j, k, t]] = power_lasso[[j, k, t]]/N
      power_group[[j, k, t]] = power_group[[j, k, t]]/N
      power_ridgeless[[j, k, t]] = power_ridgeless[[j, k, t]]/N
      power_subset[[j, k, t]] = power_subset[[j, k, t]]/N
      power_SCAD[[j, k, t]] = power_SCAD[[j, k, t]]/N
      power_MCP[[j, k, t]] = power_MCP[[j, k, t]]/N
    }
  }
}

sink("basic_1.txt")  
cat("n:", n)
cat("d:", d)
cat("M:", M)
cat("N:", N)
cat("sigma_values:", sigma_all)
cat("signal_values:", signal)
cat("order: oracle, lasso, group, ridgeless, subset, SCAD, MCP")
sink()

sink("output_1.txt")
matrix_MC = cbind(as.matrix(power_oracle[, 1],size = c(nsignal ,1)), as.matrix(power_lasso[, , 1], size = c(nsignal, nsigma)), as.matrix(power_group[, , 1], size = c(nsignal, nsigma)), as.matrix(power_ridgeless[, , 1], size = c(nsignal, nsigma)), as.matrix(power_subset[, , 1], size = c(nsignal, nsigma)), as.matrix(power_SCAD[, , 1], size = c(nsignal, nsigma)), as.matrix(power_MCP[, , 1], size = c(nsignal, nsigma)))

matrix_dCRT = cbind(as.matrix(power_oracle[, 2],size = c(nsignal ,1)), as.matrix(power_lasso[, , 2], size = c(nsignal, nsigma)), as.matrix(power_group[, , 2], size = c(nsignal, nsigma)), as.matrix(power_ridgeless[, , 2], size = c(nsignal, nsigma)), as.matrix(power_subset[, , 2], size = c(nsignal, nsigma)), as.matrix(power_SCAD[, , 2], size = c(nsignal, nsigma)), as.matrix(power_MCP[, , 2], size = c(nsignal, nsigma)))


matrix_est = cbind(as.matrix(power_oracle[, 3],size = c(nsignal ,1)), as.matrix(power_lasso[, , 3], size = c(nsignal, nsigma)), as.matrix(power_group[, , 3], size = c(nsignal, nsigma)), as.matrix(power_ridgeless[, , 3], size = c(nsignal, nsigma)), as.matrix(power_subset[, , 3], size = c(nsignal, nsigma)), as.matrix(power_SCAD[, , 3], size = c(nsignal, nsigma)), as.matrix(power_MCP[, , 3], size = c(nsignal, nsigma)))
print(rbind(matrix_MC, matrix_dCRT, matrix_est))
sink()
