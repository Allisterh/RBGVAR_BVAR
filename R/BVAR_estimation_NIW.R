BVAR_estimation_NIW = function(data, lags = 4, reps = 1000, burn = 500, const = T, trend = F, trend_qua = F, star = NULL, ex = NULL, ex_lag = NULL, xlags = 0,dummy = NULL,
                               
                               lamda = 1, tau = 0, epsilon = 0.1 , epsilonxl=0.5, epsilonx = 0.5, TVV = FALSE, feedback = F) {
  

  lamdaP    = lamda;
  tauP      = tau;
  epsilonP  = epsilon; # When high -> constant prior tends to 0.
  epsilonXP = epsilonx # When high -> constant prior tends to 0.
  epsilonXLP= epsilonxl
  
  DTA     =   data
  L       =   lags
  XL      =   xlags
  STAR    =   star
  EX      =   ex
  EXLAG   =   ex_lag
  TREND   =   trend
  TRENDQ  =   trend_qua
  CONST   =   const
  DUMMY   =   dummy
  FDBACK  =   feedback
  MAXTRIES=   10^6
  
  vardata = vardata_form(data      =  DTA,
                         lags      =  L ,
                         const     =  CONST, 
                         trend     =  TREND, 
                         trend_qua =  TRENDQ, 
                         star      =  STAR,
                         xlags     =  XL, 
                         dummy     =  DUMMY,
                         ex        =  EX,
                         feedback_only = FDBACK
  )
  
  
  # Take the arguments from the
  n          <-  vardata$number_of_endogenous
  nexo       <-  vardata$number_of_exogenous
  L          <-  vardata$number_of_lags
  t_lag      <-  vardata$number_of_obs_lags
  draws      <-  reps - burn
  dum        <-  dummy_obs_prior(data = vardata, lamda = lamdaP ,tau = tauP , epsilon = epsilonP, epsilonx = epsilonXP, epsilonxl =epsilonXLP )
  
  # conditional mean of the VAR coefficients
  X0 = dum$X_dummy; Y0 = dum$Y_dummy
  
  mstar= matrix(ols(x = X0,y = Y0))  #ols on the appended data
  
  ixx= solve(crossprod(X0))  # inv(X0'X0) to be used later in the Gibbs sampling algorithm
  
  sigma= diag(n); #starting value for sigma
  
  
  CM = array(0, list(n * L , n * L, draws))
  Sigma = array(0, list(n, n, draws))
  res = array(0, list(t_lag, n, draws))
  fit = array(0, list(t_lag, n, draws))
  beta_s = array(0, list((n*L + nexo),n, draws), dimnames=list( colnames(vardata$x_lhs) , vardata$names_of_endog_variables , NULL ))
  
  stab = 0
  i=1
  while (i <= reps){ # while 1

    if ((i%%20) == 0 ){
      print(paste0("This is my output of iteration ", i, "over",reps,"."))
    }
    
    if (i <= burn){ # if 1
        vstar = sigma %x% ixx
        beta =  mstar + t(matrix(rnorm(n*(n*L+nexo)), nrow = 1) %*% chol(vstar))
        #beta =  t(mvtnorm::rmvnorm(mean = mstar, sigma = vstar ,n = 1))
        
        # draw covariance
        e = Y0 - X0 %*% matrix(beta,nrow = n*L+nexo, ncol = n)
        # 
        # if(i == 1){
        #   b_w = tvReg::bwCov(as.matrix(e))
        # }  
        # 
        # time varying volatility
        # if (TVV){
        #   
        #   Sigma.hat = tvReg::tvCov( x = as.matrix(e),est = "ll")
        #   #sigma_tv = matrix(0, ncol = n, nrow = t_lag)
        #   # for (ii in 1:t_lag){
        #   #   sigma_tv[ii, ] = diag(Sigma.hat[,,ii])
        #   #
        #   # }
        #   #plot(y = sigma_tv[,5], x = 1:t_lag)
        #   
        #   ##The mean over time of all estimates
        #   scale = apply(Sigma.hat, 1:2, median)
        #   
        # } else {
          
          scale = crossprod(e)
          
        #}
        
        sigma = MCMCpack::riwish( (t_lag + nrow(Y0)) , scale)  #flag    
        i = i + 1
        
    } 
    
    if (i > burn){ # if 2
      
      mxtries = 0
      
      while ( mxtries <= MAXTRIES ){ #while 2
          
            vstar = sigma %x% ixx
            beta =  mstar + t(matrix(rnorm(n*(n*L+nexo)), nrow = 1) %*% chol(vstar))
            
            # draw covariance
            e = Y0 - X0 %*% matrix(beta,nrow = n*L+nexo, ncol = n)
            # time varying volatility
            # if (TVV){
            #   
            #   Sigma.hat = tvReg::tvCov( x = as.matrix(e),est = "lc", bw = b_w)
            #   #sigma_tv = matrix(0, ncol = n, nrow = t_lag)
            #   # for (ii in 1:t_lag){
            #   #   sigma_tv[ii, ] = diag(Sigma.hat[,,ii])
            #   #
            #   # }
            #   #plot(y = sigma_tv[,5], x = 1:t_lag)
            #   
            #   ##The mean over time of all estimates
            #   scale = apply(Sigma.hat, 1:2, mean)
            #   
            # } else {
              
              scale = crossprod(e)
              
            #}
            
            
            sigma = MCMCpack::riwish( (t_lag + nrow(Y0)) , scale)  #flag
            
            beta_pos =  matrix(beta,nrow = n*L+nexo, ncol = n)
            CM_pot   =  rbind(t(beta_pos[(1+nexo):nrow(beta_s), ]), cbind(diag(n*L-n), matrix(0,n*L-n,n)))
            eig = eigen(CM_pot)[[2]] 
            if (max(Mod(eig)) < 1) { #if3
              mxtries <- MAXTRIES+1
              stab = stab + 1
              beta_s[, , i-burn ] = beta_pos
              CM[, , i-burn ]  =  CM_pot
              Sigma[, , i-burn] = sigma
              res[, , i-burn ] =  e[1:t_lag,]
              i = i + 1
             } else {
              mxtries = mxtries + 1   
             }# if3
            
      } #while 2
        
    #  i = i + 1
     
     } #if 2
     
} # while 1
  
  #beta_s[, , i-burn ] = matrix(beta,nrow = n*L+nexo, ncol = n)
  #CM[, , i-burn ]  =   rbind(t(beta_s[(1+nexo):nrow(beta_s), , i-burn]), cbind(diag(n*L-n), matrix(0,n*L-n,n)))
  arguments = list("tau" = tau, "lamda" = lamda, "epsilon" =  epsilon)
  bvar = list(beta_s, Sigma, CM, res, fit, stab, vardata, draws, arguments)
  names(bvar) = c("posterior","Sigma", "CM", "res", "yfit","stability", "vardata", "draws","priors")
  return(bvar)
  
  
}


# lamdaP    = 0.5;
# tauP      = 0;
# epsilonP  = 0.001; # When high -> constant prior tends to 0.
# epsilonXP = 0.5 # When high -> constant prior tends to 0.
# epsilonXLP= 0.5
# 
# DTA     =   df1$dfr$ES[,5:10]
# L       =   1
# XL      =   1
# STAR    =   df1$dfr$ES[,11:14]
# EX      =   NULL
# EXLAG   =   NULL
# TREND   =   T
# TRENDQ  =   T
# CONST   =   T
# DUMMY   =   NULL
# FDBACK  =   T
# MAXTRIES=   20000
# 
# reps = 10000
# burn = 2000

