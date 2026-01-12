library(glmnet)
library(gglasso)
library(sparsegl)
library(MASS)
library(ncvreg)
library(ggplot2)
source("lasso_inference.r")

n = 50 
d = 200
tol = 1e-6
N = 500

s_theta = 5
s_xi = 5
theta_0 = c(rep(4, s_theta), rep(0, d-s_theta)) * sqrt(d)
xi_0 = c(rep(0.2, s_xi),rep(0,d-s_xi) ) * sqrt(d)
# group = rep(1:(d/10), each = 10)

sigma_all = seq(from = 0.4, to = 0.6, by = 0.1)
nsigma =length(sigma_all)
signal = seq(from = 0, to = 1, by = 0.2)
nsignal = length(signal)

arr <- 1:d
diff <- abs(outer(arr, arr, "-"))
q = 0.2
Sigma_matrix <- q ^ diff

estimate_lasso = function(X_ran, Z){
  lasso_model = cv.glmnet(Z, X_ran, alpha = 1, intercept = FALSE)
  best_lambda_lasso = lasso_model$lambda.min
  lasso_best_model = glmnet(Z, X_ran, alpha = 1, lambda = best_lambda_lasso, intercept = FALSE)
  lasso_coef = coef(lasso_best_model)
  mu_hat = Z %*% as.vector(lasso_coef)[-1]
  
  return (list(theta = as.vector(lasso_coef)[-1], mu = mu_hat))
}

estimate_group = function(X_ran, Z, groups){
  cv_fit <- cv.sparsegl(Z, X_ran, groups, intercept = FALSE)
  theta = as.vector(coef(cv_fit))[-1]
  mu_hat = Z %*% theta
  
  return (list(theta = theta, mu = mu_hat))
}

estimate_subset = function(X_ran, Z, supp){
  abess_fit = abess(Z, X_ran,fit.intercept = FALSE)
  abess_coef= coef(abess_fit, support.size = supp)
  mu_hat = Z %*% as.vector(abess_coef)[-1]
  
  return (list(theta = as.vector(abess_coef)[-1], mu = mu_hat))
}

# estimate_ridgeless = function(X_ran, Z) {
#   ZtZ_inv = ginv(t(Z) %*% Z)  # Moore-Penrose pseudo-inverse
#   theta = as.vector(ZtZ_inv %*% t(Z) %*% X_ran)
#   
#   theta
# }

estimate_SCAD = function(X_ran, Z){
  scad_fit <- ncvreg(Z, X_ran, penalty = "SCAD")
  cv_scad <- cv.ncvreg(Z, X_ran, penalty = "SCAD")
  best_lambda_scad <- cv_scad$lambda.min
  scad_coef <- coef(scad_fit, lambda = best_lambda_scad)
  mu_hat = scad_coef[1] + Z %*% as.vector(scad_coef)[-1]
  
  return (list(theta = as.vector(scad_coef)[-1], mu = mu_hat))
}

estimate_MCP = function(X_ran, Z){
  mcp_fit <- ncvreg(Z, X_ran, penalty = "MCP", intercept= FALSE)
  cv_mcp <- cv.ncvreg(Z, X_ran, penalty = "MCP")
  best_lambda_mcp <- cv_mcp$lambda.min
  mcp_coef <- coef(mcp_fit, lambda = best_lambda_mcp)
  mu_hat = mcp_coef[1] + Z %*% as.vector(mcp_coef)[-1]
  
  return (list(theta = as.vector(mcp_coef)[-1], mu = mu_hat))
}


Test_glmnet = function(X, d, Z, Y){
  # Set the penalty factor for the covariates
  penalty_factor = rep(0, (d+1))
  penalty_factor[2:(d+1)] = 1
  X_cbind = cbind(X, Z)
  fit = glmnet(X_cbind, Y, alpha = 0.7, lambda = 0.2, penalty.factor = penalty_factor, intercept = FALSE)
  abs(coef(fit)[2])
}

pval_lasso_oracle_compute = function(d, X, Y, Z, mu_0, N = 100){
  T_0 = Test_glmnet(X, d, Z, Y)
  
  # oracle method
  Xcopies = list()
  for(m in 1:N){
    Xcopies[[m]] = mu_0 + rnorm(n)
  }
  
  Tcopies = rep(0, N)
  for (m in 1:N){
    Tcopies[m] = Test_glmnet(Xcopies[[m]], d, Z, Y)
  }
  pval =
    (1+sum(Tcopies>= T_0))/(1+N)
  
  pval
}

pval_lasso_compute = function(d, X, Y, Z, mu_hat, X_ran, sigma, N = 100){
  T_0 = Test_glmnet(X, d, Z, Y)
  
  var = sigma^2/(1+sigma^2)
  mu = var * (mu_hat) + (1-var) * X_ran
  Xcopies = list()
  for(m in 1:N){
    Xcopies[[m]] = mu + sqrt(var)*(rnorm(n))
  }
  Tcopies = rep(0, N)
  for (m in 1:N){
    Tcopies[m] = Test_glmnet(Xcopies[[m]], d, Z, Y)
  }
  pval =
    (1+sum(Tcopies>= T_0))/(1+N)
  
  
  pval
}

pval_oracle_compute = function(X, Y, mu_0){
  Tvalue = t(Y) %*% (mu_0-X)/ sqrt(sum(Y^2))
  pval =  pnorm(Tvalue)
  
  pval
}
# 
pval_compute = function(X, Y, mu_hat, X_ran, sigma){
  var = sigma^2/(1+sigma^2)
  mu = var * (mu_hat) + (1-var) * X_ran
  Tvalue = t(Y) %*% (mu-X)/ sqrt(sum(Y^2) * var)
  pval =  pnorm(Tvalue)
  pval
}


pval_debiased = function(X, Y, m = NULL, A = NULL, lambda = NULL,  intercept = FALSE, resol = 1.3, maxiter= 50, threshold = 1e-2, verbose = TRUE){
  p <- ncol(X)
  n <- nrow(X)
  col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X))
  X <- X %*% diag(col.norm)
  htheta <- Lasso (X,Y,lambda=lambda,intercept=intercept);
  
  sigma.hat <- (1/n)*(t(X)%*%X)
  mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)))
  if(is.null(m)){
    m <- InverseLinftyOneRow(sigma.hat, i = 1, mu=mu, maxiter=maxiter, threshold=threshold)$optsol;
    A <- t(m) %*% sigma.hat %*% (m);
  }
  
  unbiased.Lasso <- htheta[1] + (m %*% t(X) %*%(Y - X %*% htheta))/n
  p.val <- 1-pnorm(sqrt(n)*abs(unbiased.Lasso)/sqrt(A))
  
  return (list(p.val = p.val, m = m , A = A))
}

Pval_oracle = Pval_des = rep(0, N * nsignal)
dim(Pval_oracle) = dim(Pval_des) = c(N, nsignal)

Error_lasso = Error_SCAD = Error_MCP = Error_subset = rep(0, N * nsigma) 
dim(Error_lasso) = dim(Error_SCAD) = dim(Error_MCP) = dim(Error_subset) =  c(N, nsigma)

Pval_lasso = Pval_group = Pval_subset = Pval_SCAD = Pval_MCP = Pval_sub_oracle = rep(0, N * nsignal * nsigma) 

dim(Pval_lasso) = dim(Pval_group) = dim(Pval_subset) = dim(Pval_SCAD) = dim(Pval_MCP) = dim (Pval_sub_oracle) = c(N, nsignal, nsigma)

starttime = Sys.time()
for (i in 1:N){
  set.seed(1000 + 2 * i)
  if (i%%10 == 0){
    print(i)
    print(Sys.time() - starttime)
    starttime = Sys.time()
  }
  indice_ran = sample(1:d)
  theta_ran = theta_0[indice_ran]
  xi_ran = xi_0[indice_ran]
  
  Z = matrix(rnorm(n * d), n, d) / sqrt(d)
  mu_0 = Z %*% theta_ran
  X = mu_0 + rnorm(n)
  W = rnorm(n)
  epsilon = rnorm(n)
  
  X_full=cbind(X,Z)
  # mat = lasso.proj(X_full, 1:n, return.Z = TRUE)$Z
  
  Res = pval_debiased(X_full, 1:n)
  m = Res$m
  A = Res$A
  for (t in 1:nsignal){
    Y =  signal[t]* X + c(Z %*% xi_ran) + epsilon
    Yhat = Y - estimate_lasso(Y, Z)$mu 
    Pval_oracle[i,t] = pval_oracle_compute(X, Yhat, mu_0)
    # Pval_des[i,t] = lasso.proj(X_full, Y, Z = mat)$pval[1]
    # Pval_oracle[i,t] = pval_lasso_oracle_compute(d, X, Y, Z, mu_0)
    
    Pval_des[i,t] = pval_debiased(X_full, Y, m = m, A = A)$p.val
  }
  for (k in 1:nsigma){
    sigma = sigma_all[k]
    X_ran = X + sigma * W
    
    # estimation
    mu_hat_lasso = estimate_lasso(X_ran, Z)
    mu_hat_subset = estimate_subset(X_ran, Z, supp = s_theta)
    mu_hat_SCAD = estimate_SCAD(X_ran, Z)
    mu_hat_MCP = estimate_MCP(X_ran, Z)
    Error_lasso[i,k] = sum((mu_hat_lasso$mu-mu_0)^2)
    Error_SCAD[i,k] = sum((mu_hat_SCAD$mu-mu_0)^2)
    Error_MCP[i,k] = sum((mu_hat_MCP$mu-mu_0)^2)
    Error_subset[i,k] = sum((mu_hat_subset$mu-mu_0)^2)
    # different magnitude
    for (j in 1:nsignal){
      Y =  signal[j]* X + c(Z %*% xi_0) + epsilon
      Yhat = Y - estimate_lasso(Y, Z)$mu 
      Pval_lasso[i,j,k] = pval_compute(X, Yhat, mu_hat_lasso$mu, X_ran, sigma)
      Pval_SCAD[i,j,k] = pval_compute(X, Yhat, mu_hat_SCAD$mu, X_ran, sigma)
      Pval_MCP[i,j,k] = pval_compute(X, Yhat, mu_hat_MCP$mu, X_ran, sigma)
      Pval_subset[i,j,k] = pval_compute(X, Yhat, mu_hat_subset$mu, X_ran, sigma)
      Pval_sub_oracle[i,j,k] = pval_compute(X, Yhat, mu_0, X_ran, sigma)
      # Pval_lasso[i,j,k] = pval_lasso_compute(d, X, Y, Z, mu_hat_lasso$mu, X_ran, sigma)
      # Pval_SCAD[i,j,k] = pval_lasso_compute(d, X, Y, Z, mu_hat_SCAD$mu, X_ran, sigma)
      # Pval_MCP[i,j,k] = pval_lasso_compute(d, X, Y, Z, mu_hat_MCP$mu, X_ran, sigma)
      # Pval_subset[i,j,k] = pval_lasso_compute(d, X, Y, Z, mu_hat_subset$mu, X_ran, sigma)
      # Pval_sub_oracle[i,j,k] = pval_lasso_compute(d, X, Y, Z, mu_0, X_ran, sigma)
    }
  }
}

Mean_error_lasso = apply(Error_lasso, c(2), function(x) sum(x)/N)
Mean_error_SCAD = apply(Error_SCAD, c(2), function(x) sum(x)/N)
Mean_error_MCP = apply(Error_MCP, c(2), function(x) sum(x)/N)
Mean_error_subset = apply(Error_subset, c(2), function(x) sum(x)/N)

print("error_lasso")
print(Mean_error_lasso)
print("error_SCAD")
print(Mean_error_SCAD)
print("error_MCP")
print(Mean_error_MCP)
print("error_subset")
print(Mean_error_subset)

alpha = 0.05

power_oracle = apply(Pval_oracle, c(2), function(x) sum(as.numeric(x<alpha))/N)
power_des = apply(Pval_des, c(2), function(x) sum(as.numeric(x<alpha))/N)
power_lasso = apply(Pval_lasso, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
# power_ridgeless = apply(Pval_ridgeless, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
power_subset = apply(Pval_subset, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
power_group = apply(Pval_group, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
power_SCAD = apply(Pval_SCAD, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
power_MCP = apply(Pval_MCP, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)
power_sub_oracle = apply(Pval_sub_oracle, c(2,3), function(x) sum(as.numeric(x<alpha))/ N)

print("power_oracle")
print(power_oracle)
print("power_des")
print(power_des)
print("power_lasso")
print(power_lasso)
print("power_SCAD")
print(power_SCAD)
print("power_MCP")
print(power_MCP)
print("power_subset")
print(power_subset)

for (i in 1:length(sigma_all)){
  plot_data <- data.frame(
    signal = rep(signal, times = 6),
    power = c(unlist(power_oracle), unlist(power_des), unlist(power_lasso[,i]),  unlist(power_SCAD[,i]), unlist(power_MCP[,i]),  unlist(power_subset[,i])),
    method =rep(c("oracle", "debiased-lasso", "lasso",  "SCAD", "MCP", "subset"), each = 6)) 
  plot_data$method <- factor(plot_data$method, levels = c("oracle", "debiased-lasso", "lasso", "SCAD", "MCP", "subset"))
  plot_data$delta <- sqrt(plot_data$power * (1 - plot_data$power) / N) 
  p <- ggplot(plot_data, aes(x = signal, y = power, linetype = method, color = method)) +
    geom_line(aes( color = method), size = 0.5) +
    geom_ribbon(aes(ymin = power - delta, ymax = power +  delta, fill = method), alpha = 0.3, color = NA) +
    scale_color_manual(values=c("#D8C033","#B9A2A3","#225F6D","#B0946F",
                                "#D8C033","#B9A2A3","#225F6D")) +
    labs(x = expression(beta[0]), y = "power") +
    xlim(0, 1) +
    ylim(0, 1) +
    geom_hline(yintercept = 0.1, linetype = "dotted", color = "red", size = 0.5) +
    theme_bw() +
    theme(
      legend.text = element_text(size = 14),  # Font size for legend text
      legend.title = element_text(size = 16),  # Font size for legend title
      axis.title.x = element_text(size = 16),  # Font size for x-axis label
      axis.title.y = element_text(size = 16),  # Font size for y-axis label
      axis.text.x = element_text(size = 14),  # Font size for x-axis tick labels
      axis.text.y = element_text(size = 14)  # Font size for y-axis tick labels
    )
  print(p) 
  myfile=paste("myplot",sigma_all[i],".pdf",sep="")
  ggsave( p, filename =myfile, height = 5, width = 7, units = "in")
}

data <- data.frame(
  sigma = rep(sigma_all, 6),
  type1_error = c(
    rep(power_oracle[1], nsigma), rep(power_des[1], nsigma),
    power_lasso[1,],  power_SCAD[1, ], power_MCP[1, ], power_subset[1, ]
  ), 
  method = rep(c("oracle", "debiased-lasso", "lasso","SCAD", "MCP", "subset"), each = length(sigma_all))
)
data$method = factor(data$method, levels = c("oracle", "debiased-lasso", "lasso","SCAD", "MCP", "subset"))

data$delta <- sqrt(data$type1_error * (1 - data$type1_error) / N) 
color_set <- c("#D8C033","#B9A2A3","#225F6D","#B0946F",
               "#D8C033","#B9A2A3","#225F6D")
plot <- ggplot(data, aes(x = sigma, y = type1_error, color = method)) +
  geom_line(aes(linetype = method), size = 0.5) +
  geom_ribbon(aes(ymin = type1_error - delta, ymax = type1_error +  delta, fill = method), alpha = 0.3, color = NA) +
  labs(x = expression(sigma), y = "Type I Error") + ylim(0, 1) +
  theme_minimal() + 
  scale_color_manual(values = color_set) +
  scale_fill_manual(values = color_set) +
  geom_hline(yintercept = 0.1, linetype = "dotted", color = "red", size = 0.5) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 14),  # Font size for legend text
    legend.title = element_text(size = 16),  # Font size for legend title
    axis.title.x = element_text(size = 16),  # Font size for x-axis label
    axis.title.y = element_text(size = 16),  # Font size for y-axis label
    axis.text.x = element_text(size = 14),  # Font size for x-axis tick labels
    axis.text.y = element_text(size = 14)  # Font size for y-axis tick labels
  )
print(plot)
ggsave('type1_sparse_d100_147.pdf', plot, width = 7, height = 5, units = "in")

save(Pval_oracle, file = "pval_oracle.RData")
save(Pval_des, file = "pval_des.RData")
save(Pval_lasso, file = "pval_lasso.RData")
save(Pval_SCAD, file = "pval_SCAD.RData")
save(Pval_MCP, file = "pval_MCP.RData")
save(Pval_subset, file = "pval_subset.RData")
save(Pval_sub_oracle, file = "pval_sub_oracle.RData")
