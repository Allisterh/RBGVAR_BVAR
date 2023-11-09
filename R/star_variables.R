#   This is an auxillary Fuction that builds the star variables for every country
#   Data should be submitted in stacked cross-sections format. 
#   There must be a uniqe identifier for time and crosssectional unit
#   The first column should be the Time ID and the second the Cross-section ID
#   Weight matrix should be submitted in a matrix form with country colnames and rownames 
#   being the same with the cross-sectional identifier.
#   Data shoulld not contain any NAs. Function cannot handle NAs.
#   The function will return a list of 3 elements.
#   The first two element will be the W and S(election) matrices for the Solution of the GVAR 
#   the third will be a list of the star variables for every cross-section.
#   data = gv
star_variables = function(data, w_mat, lags = 1, universal_variables = NULL){
  
   #universal_variables = c("lr","r", "UNRATE", "FE", "Tax", "EQ")
   #data = gv;w_mat
  names(data)[1:2] = c("timeID", "crossID")
  
  GV = data %>% 
    pivot_longer(cols = -c(timeID, crossID), names_to = "Variable") %>% 
    mutate(crossID = factor(crossID, rownames(w_mat))) %>% 
    arrange(crossID, timeID) %>% 
    pivot_wider(names_from = c(crossID,Variable), names_sep = "___") %>% 
    select_if(~ !any(is.na(.)))
   # select_if(~sum(!is.na(.)) > 0)
  
  #  This is he lagged gobal vector
  #GV_lag =  VAR_data_prep(GV[,-1], lags = lags)
  #GV = cbind(GV[-lags,1], GV_lag$LHS)
  
  nglobal = ncol(GV[,-1])
  names_nglobal = names(GV[,-1])
  # Processing the names pf the global vector
  NAMES = data.frame( "Names" = names(GV[,-1]) ) %>% 
    separate(Names, into = c("Country", "Variable"), remove = F, sep = "___")
  
  W = list()
  S = list()
  Star_variables = list()
  uniq_variables = list()
  nrowsW = c()
  nrowsS = c()
  nvariables = list()
  k = 0
  #i = "AGR"
  for (i in unique(NAMES$Country)){
    # Selection matrix S
    nvariables[[i]] = length(unique(NAMES$Variable[NAMES$Country == i]))
    S[[i]] = matrix(0, ncol = nglobal, nrow = nvariables[[i]])
    S[[i]][,(1+k):(nvariables[[i]]+k)] = diag(nvariables[[i]])
    rownames(S[[i]]) = NAMES$Variable[NAMES$Country == i]
    colnames(S[[i]]) = names_nglobal
    k = k + nvariables[[i]]
    # matrix W
    # Unique variables after I drop the variables of the country I am interested in
    uniq_variables[[i]] = unique(NAMES$Variable[!NAMES$Country == i])    
    nrowsW[i] = length(uniq_variables[[i]])
    W[[i]] = matrix(0, nrow = nrowsW[[i]], ncol = nglobal )
    rownames(W[[i]]) = uniq_variables[[i]]
    for (j in uniq_variables[[i]] ) {
      #j = "y"
      position = which(NAMES$Variable == j)
      country = NAMES$Country[which(NAMES$Variable == j)]
      W[[i]][j, position ] = w_mat[i, country] 
      
      # Universal variables
      if (!is.null(universal_variables)) {
        if (j %in% universal_variables) {
          W[[i]][j, position ] = 1  
        }
      }
      
    }
    rownames(W[[i]]) = paste0(rownames(W[[i]]), "_s")
    # Star variables
    Date = GV[, 1]
    Star_variables[[i]] = t(W[[i]] %*% t(GV[,-1]))
    #namesstar = str_c(i, colnames(Star_variables[[i]]), "s", sep = "_")
    Star_variables[[i]] = data.frame(cbind(Date,Star_variables[[i]]))
    #names(Star_variables[[i]]) = c("timeID", namesstar)
  }
  
  return(list("W" = W,"S" = S, "Star_variables"= Star_variables, "GV" = GV))
  
}
