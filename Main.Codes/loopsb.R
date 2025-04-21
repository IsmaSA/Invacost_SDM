   for(m in 1:length(pred.l)){ #loop over each model
  
    print(m)  
    model_id <- m
    pred.list <- unlist(pred.l[[m]]) #to obtain the correct predictors
    # Add b0 and b1 to the initial parameters  NEW (!!!!!!!!!!!)
    print( pred.list)

    if("time" %in% pred.list)
    {
    	par1 <- c(rep(0, length(pred.list)), sd(eu.df1$Cost_estimate_per_year_2017_USD_exchange_rate, na.rm=T), b0 = .1)#,s=0)#10) ######### NEW ######### (bo and b1)
	    names(par1) <- c(pred.list, "sd", "b0")#,"s")   ######### NEW ######### (bo and b1)
    }else{
    	par1 <- c(rep(0, length(pred.list)), sd(eu.df1$Cost_estimate_per_year_2017_USD_exchange_rate, na.rm=T)) ######### NEW ######### (bo and b1)
	    names(par1) <- c(pred.list, "sd")   ######### NEW ######### (bo and b1)
    
    }
        
    #choose the best of the previous parameter - using the one with the most similar variables, and the best likelihood
	len=length(pred.list)
    if(len>1)
    {
    	mlk= sapply(1:(m-1),function(m1){
		#take only the simpler models where all terms will be subset of more complex model 
		x=pred.l[[m1]]
	    	if(length(which(x %in% pred.list))==length(x))
	    	{
	    		return(lst[[m1]]$Likelihood)
	    		
	    	}else{
			return(NA)    	
		}
    	})
    	tmpp=params[[which.min(mlk)]]
    	par1[names(tmpp)]=tmpp
#    	if("time2" %in% pred.list) #ONLY NEEDED FOR THRESHOLD MODEL
#    	{
#    		if(par1["time2"] < par1["time"])
#    			par1["time2"]=par1["time"]
#    	}
    }
     op.gamma=rep_optim(par1,df.fit.full2, pred.list)  # include pred.list or give error
    
      params[[m]]=op.gamma$par
#    }
	
    gamma1=op.gamma$par
#    names(gamma1) <- names(par1) #c(pred.list, "sd", "b0", "b1") ######### NEW ######### (bo and b1)
#    gamma1
    ######### NEW #########
    #print(paste("Convergence code:", op.gamma$convergence))
    ######### NEW #########
    
    
    ##Metrics of fit --> calculate AIC for the given set of parameters gamma
    
    # log-likelihood:, 
    #fit.t1 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, extrap = F, fit.test=T) 
    fit.t1 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, bound = NULL, ret.dat =T) 
    
    # fit.t11 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, extrap = F, fit.test=F) 
    fit.t11 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, bound = NULL, ret.dat =F) 
    likelihoods <- fit.t11
    v.AIC=2*likelihoods + 2 * length(gamma1)  
    
    
    lst[[m]] <- data.frame(
      "AIC" = v.AIC    ,
      "Likelihood" = likelihoods ,
      "Type_cost" = i, 
      "Model" = m                ,
      'predictors' = paste(pred.list, collapse = ";"),
      "Full.Partial_model" = "Partial"
    )
  }
results_df <- do.call(rbind,lst)
saveRDS(results_df,paste0("result_",oname,nm_outlier,i,".rds"))
params1=lapply(params, function(x){
full_gamma2[names(x)]=x
return(full_gamma2)})

params2[[i]]=do.call(rbind,params1)
saveRDS(params2[[i]], paste0("params_",oname,nm_outlier,i,".rds"))
lst1[[i]]=results_df

