Stack_attack = function(big_df, star) {
  # big_df = df1; ith_draw = 1; i = "AGR"; l=1
  # Arguments taken from built structures NOT varying across countries - Auxiliary data structures
  draws = big_df$var_mod[[1]]$draws
  # Names for the position of country-variables in the global vector
  NAMES = data.frame("Names" = names(GV[, -1])) %>% separate(Names, into = c("Country", "Variable"), remove = F, sep = "___" )
  # Number of lags of the endogenous and exogenous variables
  nxlags = big_df$var_mod[[1]]$vardata$number_of_xlags
  nlags = big_df$var_mod[[1]]$vardata$number_of_lags
  
  # Auxiliary-In each draw collects the G matrix 
  Gith = list()
  
  # final outputs structure
  #Global Covariance
  Smat_global = array(0, list(nvar_global, nvar_global, draws), dimnames = list(names(star$GV[, -1]), names(star$GV[, -1]) , NULL))
  
  # G inverse matrix - Coefficients of contemporaneous coefficients
  Ginv = array(0, list(nvar_global, nvar_global, draws), dimnames=list(names(star$GV[,-1]), names(star$GV[,-1]) , NULL  ))
  
  # Global Companion Matrix
  #CM_global = array(0,list((nvar_global * nlags), nvar_global * nlags, draws))
  Hmat = array(0,list(nvar_global, nvar_global,nlags, draws), dimnames = list(names(star$GV[, -1]), names(star$GV[, -1]) ,NULL, NULL))
  
  # Hmat[,,,1]
  for (ith_draw in 1:draws) {
    
    for (i in big_df$iso2c) {
      
      # Gobal covariance matrix - Block Diagonal
      position =  which(NAMES$Country == i)
      Smat_global[position, position, ith_draw] =  big_df$var_mod[[i]]$Sigma[, , ith_draw]
      
      
      #crossprod(big_df$residuals[[i]][, , ith_draw])
      
      #   G  Matrix   # 
      #names of laged coefficients
      names_lamda = big_df$var_mod[[i]]$vardata$names_of_star_variables
      # Number of equatioins per country
      n = big_df$var_mod[[i]]$vardata$number_of_endogenous
      # Selection matrix for each country
      W =  rbind(star$S[[i]], star$W[[i]][names_lamda, ])
      
      # Bind the Identity matrix In with Lamda0 
      I = diag(n)
      # Lamda0 Matrix
      if (length(names_lamda) == 1) {
        LAMDA = big_df$var_mod[[i]]$posterior[grep(paste0("_s_", 0), rownames(big_df$var_mod[[i]]$posterior[,,1])), ,ith_draw]
      } else {
        LAMDA = t(big_df$var_mod[[i]]$posterior[grep(paste0("_s_", 0), rownames(big_df$var_mod[[i]]$posterior[,,1])), ,ith_draw])
      }
      
      if(length(LAMDA)==0){
        LAMDA = matrix(0, nrow = nrow(LAMDA), ncol = (nrow(W) - ncol(I)) ) 
      } 
      # big_df$var_mod[[i]]$posterior[,,1]
      # A0 matrix 
      A0 = as.matrix(cbind(I , -LAMDA))
      
      # G matrix for draw ith
      Gith[[i]] =  A0 %*% W
      #View(t(Gith[[i]] %*% t(GV[,-1])))
      # # # # # # # # # # # # # # # # # # # # # #
      
      #    H Matrix    #
      # Collects the lag coeffs for each country    
      
      names_coefficients = rownames(big_df$var_mod[[i]]$posterior[, , 1])
      names_phi = big_df$var_mod[[i]]$vardata$names_of_endogenous
      nexo = big_df$var_mod[[i]]$vardata$number_of_exogenous
      #l =1
      for (l in 1:nlags){
        
        if (l <= nxlags) {
          Lamda = big_df$var_mod[[i]]$posterior[grep(paste0("_s_", l), names_coefficients), , ith_draw]
          Phi = big_df$var_mod[[i]]$posterior[nexo + grep(paste0("_", l), names_coefficients[-c(1:nexo)]), , ith_draw]
          
          if (is.array(Lamda)){
            B1 = as.matrix(cbind(t(Phi) , t(Lamda)))
            Hlith = B1 %*% W
          } else  {
            B1 = as.matrix(cbind(t(Phi) ,Lamda))
            Hlith = B1 %*% W
          }
          
          rownames(Hlith) = paste(i, rownames(Hlith), sep = "___")
          #position = grep(pattern = i, rownames(Hmat[,,1,1]))
          
        } else {
          
          Phi = big_df$var_mod[[i]]$posterior[nexo + grep(paste0("_", l), names_coefficients[-c(1:nexo)]), , ith_draw]
          Lamda = matrix(0, ncol = nrow(Phi), nrow = length(names_lamda))
          
          rownames(Lamda) = paste(names_lamda, l, sep = "_")
          colnames(Lamda) = rownames(Phi)
          
          if (is.array(Lamda)){
            B1 = as.matrix(cbind(t(Phi) , t(Lamda)))
            Hlith = B1 %*% W
          } else  {
            B1 = as.matrix(cbind(t(Phi) ,Lamda))
            Hlith = B1 %*% W
          }
          
          rownames(Hlith) = paste(i, rownames(Hlith), sep = "___")
          #position = grep(pattern = i, rownames(Hmat[,,1,1]))
          
        } # end if l<=nxlags
        
        Hmat[position, ,l,ith_draw] = Hlith
        
      } # lags loop
      
      
      # # # # # # # # # # # # # # # # # # # # # #
      
      
    } # cross section loop
    
    Ginv[,,ith_draw] = solve(do.call(rbind, Gith))
    
  } # draws loop
  
  ret = list(Ginv, Smat_global, Hmat)
  names(ret) = c("Ginv", "Smat_Global", "Hmat")
  return(ret)
  
} # stack_attack



