BVAR_draw = function(data, lags = lags, const = T, trend = F, trend_qua = F, ex = NULL, ex_lag = NULL, xlags = 0, chol = T, draws = draws) {
  # data = df[["FR"]][,4:8] ; lags =1; draws = 10;const= T;trend =T; trend_qua = F; ex = NULL; ex_lag = NULL; xlags = 0 
  # Starts the analysis
  # Run a simple VAR
  var = VAR_estimate(data , lags = lags, const = const, trend = trend, trend_qua = trend_qua,ex = ex, ex_lag = ex_lag, xlags = xlags)
  
  # Take the argument from the VAR
  n = var$data_final$number_of_variables
  #const = var$const 
  
  # S = crossprod(var$res)/var$data_final$df  Maybe
  S = crossprod(var$res) # SSE
  # Draw Sigma from IW ~ (S,v), where v = T-(K*p+1)-K-1 :
  Sigma = MCMCpack::riwish(v = var$data_final$dof, S)    
  
  # Compute the Kronecker product Q x inv(X'X) and vectorize the matrix
  # pi_hat in order to obtain vec(pi_hat):
  beta_variance = Sigma %x% solve(crossprod(as.matrix(var$data_final$x_lhs)))
  
  # vec(b)
  s = dim(var$bols)
  vec_beta_hat = matrix(var$bols, s[1]*s[2],1)
  
  # Draw beta from a multivariate normal distribution with mean vec(b) and variance Sigma x inv(X'X):
  beta =  t(mvtnorm::rmvnorm(mean = vec_beta_hat,sigma = beta_variance,n = 1))
  
  # reshape vec(b) such that Y=X*b+e, i.e. b is (n*lag+nexo)x(n).
  beta = matrix(beta, nrow = (n*lags + var$nexo), ncol = n)
  
  # Create the companion form representation matrix A:
  CM = rbind(t(beta[(1+var$nexo):nrow(beta),]), cbind(diag(n*lags-n), matrix(0,n*lags-n,n))) #(K*p)x(K*p) matrix
  
  #check stability
  eig = c(eigen(CM)[[1]])
  stability = max(Mod(eig)) < 1
  
  # Store residuals and fitted values:
  yfit = as.matrix(var$data_final$x_lhs) %*% beta
  res = as.matrix(var$data_final$y_lhs - yfit)
  

  bvar = list(beta, CM, Sigma,res,yfit, stability)
  names(bvar)= c("beta","CM", "Sigma","res", "y_fit", "stability")
  return(bvar)  
  
}


