library(lbfgs)
library(glmnet)
library(gglasso)
generate_parameters_gl = function(d = 100, sigma = sqrt(10), n = 50, l1 = 2, l2 = 0.01, theta0 = c(rep(0.5, 5), rep(0,100 - 5)), isnorm = 1){
  parameters = list()
  parameters$d = d 
  parameters$sigma = sigma
  parameters$l1 = l1
  parameters$l2 = l2 
  parameters$tol = 1e-6   
  parameters$signal = seq(from=0,to=1,by=0.2) 
  parameters$isnormalize = ifelse(isnorm, sqrt(parameters$d), 1) 
  parameters$theta0 = theta0*parameters$isnormalize
  
  # example specific parameters
  parameters$example = list()
  parameters$example$n = n 
  parameters$example$coefyZ = c(rep(0.1,5), rep(0,parameters$d-5))*parameters$isnormalize  
  parameters
}

generate_experiment_gl = function(isignal,parameters){
  experiment = list()
  Z = matrix(rnorm(parameters$example$n*parameters$d),parameters$example$n,parameters$d)/parameters$isnormalize
  X = c(Z%*%parameters$theta0) + rnorm(parameters$example$n)
  
  Y_mlr = parameters$signal[isignal]*X + c(Z%*%parameters$example$coefyZ) + rnorm(parameters$example$n)
  experiment$X = X 
  
  experiment$example = list()
  experiment$example$Z = Z
  
  proj_acss = solve(tcrossprod(Z, Z)*parameters$d/parameters$sigma**2 + diag(parameters$example$n))
  proj_acss_chol = t(chol(proj_acss))
  
  experiment$example$Proj_acss = proj_acss
  experiment$example$Proj_acss_chol = proj_acss_chol
  
  experiment$example$Y_mlr = Y_mlr 
  
  experiment
}


generate_X_null_gl = function(theta, parameters, experiment){
  X = c(experiment$example$Z%*%theta) + rnorm(parameters$example$n)
  X
}



generate_copies_acss = function(M,parameters,experiment, thetah){  
  mu =  experiment$example$Z%*%thetah + parameters$l2*experiment$example$Proj_acss%*%experiment$example$Z%*%thetah*parameters$d/parameters$sigma**2
  Xcopies = list()
  for(m in 1:M){
    Xcopies[[m]] = c(mu) + c(experiment$example$Proj_acss_chol%*%rnorm(parameters$example$n))
  }
  Xcopies
  
}







Tfun_glmnet = function(x,parameters,experiment){  
  # Set the penalty factor for the covariates
  penalty_factor = rep(0, parameters$d+1)
  penalty_factor[2:parameters$d+1] = 1
  X = cbind(x, experiment$example$Z)
  fit = glmnet(X, experiment$example$Y_mlr, alpha = 0.7, lambda = 0.2, penalty.factor = penalty_factor, intercept = FALSE) 
  abs(coef(fit)[2])
}

projection = function(theta,parameters,experiment){
  theta_sign = sign(theta)
  zero_index = which(theta==0)
  nonzero_index = which(abs(theta) > parameters$tol) 
  if(length(nonzero_index)>0){
    A0 =  matrix(0, ncol = parameters$d, nrow  = length(zero_index)+1)
    A0[,nonzero_index] = matrix(rep(theta_sign[nonzero_index], length(zero_index)+1), byrow = 1, nrow =  length(zero_index)+1) 
    if(length(zero_index)>0){
      A0[1: length(zero_index), zero_index] = diag(length(zero_index))
    }
  }
  else{
    A0 = diag(parameters$d)
  }
  #pro_act = t(A0)%*%solve(A0%*%t(A0))%*%A0
  pro_act = crossprod(A0,solve(A0%*%t(A0)))%*%A0
  pro_act =  pro_act * (abs(pro_act) >  parameters$tol)
  pro_inact = diag(parameters$d) - pro_act  
  
  pro = list() 
  pro$pro_act = pro_act
  pro$pro_inact = pro_inact
  
  pro
}

gradient = function(x,theta,parameters,experiment){
  c(crossprod(experiment$example$Z, experiment$example$Z%*%theta-x))  + parameters$l2*theta
} 


check_SSOSP_con = function(x,w,theta,parameters,experiment){
  # Checks the minimum eigenvalue of the Hessian is positive 
  # and that the gradient is within a small tolerance of zero
  pro = projection(theta,parameters,experiment)   
  sum((pro$pro_inact%*%(gradient(x,theta,parameters,experiment)+parameters$sigma*w))^2)<=parameters$tol
}






thetahat_con_gl = function(x,w,parameters,experiment){ 
  Lw = function(theta){sum((x - experiment$example$Z%*%theta)**2)/2 + parameters$sigma*sum(w*theta) + sum(theta**2)/2*parameters$l2}
  gradLw = function(theta){c(crossprod(experiment$example$Z, experiment$example$Z%*%theta-x))+parameters$sigma*w + parameters$l2*theta} 
  lasso_fit = glmnet(experiment$example$Z, x)
  theta_init = coef(lasso_fit, s = parameters$l1/parameters$d)[2:(parameters$d+1)]  
  thetahat = lbfgs::lbfgs(call_eval=Lw,call_grad = gradLw,vars = theta_init, invisible = 1, orthantwise_c = parameters$l1)$par 
  thetahat = thetahat*(abs(thetahat) > parameters$tol)
  thetahat
}




generate_copies_acsscon = function(M,parameters,experiment, thetah, W){
  g = gradient(experiment$X,thetah,parameters,experiment)+parameters$sigma*W
  mu =  experiment$example$Z%*%thetah + experiment$example$Proj_acss%*%experiment$example$Z%*%(parameters$l2*thetah - g)*parameters$d/parameters$sigma**2
  
  Xcopies = list()
  for(m in 1:M){
    Xcopies[[m]] = c(mu) + c(experiment$example$Proj_acss_chol%*%rnorm(parameters$example$n))
  }
  Xcopies
}




run_one_trial_gl = function(seed, M, d, sigma, n, l1,l2, theta0, isnorm){ 
  set.seed(seed) 
  parameters = generate_parameters_gl(d, sigma, n, l1,l2, theta0, isnorm)
  nsignal = length(parameters$signal)
  pval_aCSS = pval_oracle = pval_aCSS_con = rep(0,nsignal) 
  for(isignal in 1:nsignal){ 
    experiment = generate_experiment_gl(isignal,parameters) 
    Tfun_ = function(x){Tfun_glmnet(x,parameters,experiment)}
    
    # oracle method
    Xcopies = list()
    for(m in 1:M){
      Xcopies[[m]] = generate_X_null_gl(parameters$theta0, parameters, experiment)
    }
    pval_oracle[isignal] = 
      (1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+M)
    
    
    
    # aCSS method
    W = rnorm(parameters$d)/sqrt(parameters$d) 
    
    
    thetah_acss = c(solve(crossprod(experiment$example$Z, experiment$example$Z)+ parameters$l2*diag(parameters$d))%*%(crossprod(experiment$example$Z,experiment$X) - parameters$sigma*W))
    
    thetah_acss_con = thetahat_con_gl(experiment$X,W,parameters,experiment)
    
    
    Xt_acss = list()
    Xt_acss_con = list()
    
    Xt_acss = generate_copies_acss(M,parameters,experiment, thetah_acss)
    pval_aCSS[isignal] = 
      (1+sum(unlist(lapply(Xt_acss,Tfun_))>=Tfun_(experiment$X)))/(1+M)
    
    if(check_SSOSP_con(experiment$X,W,thetah_acss_con,parameters,experiment)){ 
      Xt_acss_con = generate_copies_acsscon(M,parameters,experiment, thetah_acss_con,W)
      pval_aCSS_con[isignal] = 
        (1+sum(unlist(lapply(Xt_acss_con,Tfun_))>=Tfun_(experiment$X)))/(1+M)
    }else{ 
      pval_aCSS_con[isignal] = 1
    }
    
    
    
  }
  
  pvals = cbind(pval_oracle, pval_aCSS,pval_aCSS_con)
  colnames(pvals) = c('oracle', 'aCSS','aCSS_con')
  rownames(pvals) = parameters$signal
  
  pvals
}



