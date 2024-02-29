#Likelihood function modified to include b0 and b1 (changes are noted by: ######### NEW #########, at the beginning and end)
scale_fit_fun <- function(gamma1, df, pred.list, extrap = F, fit.test =F){ #gamma1= vector of parameters to optimize
  
  full_gamma[names(gamma1)]=gamma1

  ######### NEW #########
  # Extract parameters
  b0 <- full_gamma["b0"]
  b1 <- full_gamma["time"]
  b2=full_gamma["time2"]   	
  t3=df$t2 - b2 #include a temporal threshold
  t3[t3<0]=0
  
#NEEDED FOR THRESHOLD MODEL
#  if(! ("time2" %in% pred.list))
#  {
#  	b2=b1
#  }else{
#	  
#	  if(b2 < b1)
#	  	return(NA)
#  }
#  if(b0<0 || b0 >1)
#  {	
#  	return(NA)
#  }
  ######### NEW #########
#  if(b1 != 0)	
  if("time" %in% pred.list)
  {
	#logistic version
	df$suitability2= df$suitability / (1 + exp(-(b0 + b1*t3))) #works best, but not clear why it is better than the threshold since it is basically the same thing
	
	  pos=which(is.na(df$t2))
	  df$suitability2[pos]=df$suitability[pos]
#	  if(! ("time2"  %in% pred.list))	  
#	  {
#	  	pos=which(t3>0)
#	  	df$suitability2[pos]=df$suitability[pos]
#	  }



#threshold segmented version 
#	df$suitability2=df$suitability
#	pos=which(!is.na(df$t2) & df$t2 < b1)
#	df$suitability2[pos]=df$suitability[pos]*b0
#	pos=which(!is.na(df$t2) & df$t2 >= b1 & df$t2 <b2)
#	slope=(1-b0)/(b2-b1) #rise/run
#	if(length(pos)>0)
#		df$suitability2[pos] = (b0+slope*t3[pos])*df$suitability[pos]



  }else{
	  df$suitability2=df$suitability
  }  

  ln <- pred.list #replace pred.list by ln (list.name) so that we don't overwrite it
  
    if(length(gamma1) == 1){
      names(gamma1) <- "sd"
    }
   
    scale_fun <- function(dff, fit.test){
      pos <- which(dff$type == "invacost")
      x <- dff[pos,] #invacost data
      y <- dff[-pos,] #extrapolation data
      if(nrow(y) == 1){
        y$cost_pred <- NA 
        if(fit.test ==T){ #return the full dataset to calculate metrics of fit
          return(y)
        } else {
		return(NA)
	}  
	
      } else {
        
        for(i in 1:nrow(y)){
          
          ######### NEW #########
          # For y dataset
          
          pred1 <- x$Cost_estimate_per_year_2017_USD_exchange_rate[-i] * (y$GDP[i] / x$GDP[-i])^full_gamma["GDP"] * 
            (y$Pop_size[i] / x$Pop_size[-i])^full_gamma["Pop_size"] * 
            (y$Agriculture[i] / x$Agriculture[-i])^full_gamma["Agriculture"] * 
            (y$suitability2[i]  / x$suitability2[-i])^full_gamma["suitability"]
          y$cost_pred[i] <- mean(pred1)
      }
        
        if(fit.test ==T){
          return(y)
          
        } else {
          ll <- log(dnorm(y$Cost_estimate_per_year_2017_USD_exchange_rate - y$cost_pred, sd= full_gamma["sd"]))
#	print(dff$comb2[1])
#	print(head(ll))
#          ll <- sum(log(dnorm(y$cost_mil - y$cost_pred, sd= full_gamma["sd"])))
          return(-sum(ll)) 
        }
        
    }
   }
    df.by <- by(df, df$comb2, scale_fun, fit.test= fit.test) #apply the function "scale_fun" to each "comb2" column. comb2 is the combination of species/sector/type of cost
    
    
    
    if(fit.test == T){
      df.full <- do.call(rbind, df.by) 
      #df.full<- rbind(df.by)
      cat('*')
      return(df.full)
      
    } else {
      measure.fit <- sum(df.by, na.rm = T) #sum of the log-likelihood
#	print(measure.fit)
      return(measure.fit)
    }
    
}

