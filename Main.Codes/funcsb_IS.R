#Likelihood function modified to include b0 and b1 (changes are noted by: ######### NEW #########, at the beginning and end)
scale_fit_fun <- function(gamma1, df, pred.list, bound = NULL, ret.dat =F){ #gamma1= vector of parameters to optimize
  
  full_gamma[names(gamma1)]=gamma1

  ######### NEW #########
  # Extract parameters
  b0 <- full_gamma["b0"]
  b1 <- full_gamma["time"]
  b2=full_gamma["time2"]   	
  t3=df$t2 - b2 #include a temporal threshold
  t3[t3<0]=0
  if(full_gamma["sd"]<=0)
  	return(NA)
  ######### NEW #########
  if("time" %in% pred.list)
  {
	#logistic version
	df$suitability2= df$suitability / (1 + exp(-(b0 + b1*t3))) #works best, but not clear why it is better than the threshold since it is basically the same thing
	
	  pos=which(is.na(df$t2))
	  if(length(pos)>0)
	  	df$suitability2[pos]=df$suitability[pos]

  }else{
	  df$suitability2=df$suitability
  }  

  ln <- pred.list #replace pred.list by ln (list.name) so that we don't overwrite it
  
    if(length(gamma1) == 1){
      names(gamma1) <- "sd"
    }
   
    scale_fun <- function(dff, ret.dat)
    {
      pos <- which(dff$type == "invacost")
      x <- dff[pos,] #invacost data
      y <- dff[-pos,] #extrapolation data
     if(nrow(y) == 0) #this should never happen
     	return(NULL)	

     	if(nrow(x) == 1 && nrow(y) ==1) #if the same thing, don't calculate 
	{
		if(x$Index == y$Index)
		{

			if(ret.dat ==T){ #return the full dataset to calculate metrics of fit
			  y$cost_pred=NA
			  return(y)
			} else {
				return(NA) #don't contribute to the fit
			}  
		}	
  	}
	for(i in 1:nrow(y))
	{
	  
	  ######### NEW #########
	  # For y dataset
	  pos=which(x$Index == y$Index[i])
	  	if(length(pos)==0)
			pos=nrow(x)+1
		if(is.null(bound))
		{
		  pred1 <- x$Cost_estimate_per_year_2017_USD_exchange_rate[-pos] * (y$GDP[i] / x$GDP[-pos])^full_gamma["GDP"] * 
		    (y$Pop_size[i] / x$Pop_size[-pos])^full_gamma["Pop_size"] * 
		    (y$Agriculture[i] / x$Agriculture[-pos])^full_gamma["Agriculture"] * 
		    (y$suitability2[i]  / x$suitability2[-pos])^full_gamma["suitability"]
		}else{
		  pred1 <- x$Cost_estimate_per_year_2017_USD_exchange_rate[-pos] * bound_ratio(y$GDP[i] / x$GDP[-pos],bound$GDP[1],bound$GDP[2])^full_gamma["GDP"] * 
		    bound_ratio(y$Pop_size[i] / x$Pop_size[-pos],bound$Pop_size[1],bound$Pop_size[2])^full_gamma["Pop_size"] * 
		    bound_ratio(y$Agriculture[i] / x$Agriculture[-pos],bound$Agriculture[1],bound$Agriculture[2])^full_gamma["Agriculture"] * 
		    bound_ratio(y$suitability2[i]  / x$suitability2[-pos],bound$suitability[1],bound$suitability[2])^full_gamma["suitability"]
		
		}
	  	y$cost_pred[i] <- mean(pred1)
    	}
	
	if(ret.dat ==T){
	  return(y)
	  
	} else {
	  	ll <- log(dnorm(y$Cost_estimate_per_year_2017_USD_exchange_rate - y$cost_pred, sd= full_gamma["sd"]))
	  return(-sum(ll)) 
	}
        
    	
   }
    df.by <- by(df, df$comb2, scale_fun, ret.dat= ret.dat) #apply the function "scale_fun" to each "comb2" column. comb2 is the combination of species/sector/type of cost
    
    
    
    if(ret.dat == T){
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


ratio_bound <- function(dff, pred.list){
  
  ln <- pred.list 
  pos=grep("time", ln)
  if (length(pos)>0) {
    ln <- ln[-pos]
  }
  
  bound=as.data.frame(matrix(NA, ncol=length(ln),nrow=2))
  names(bound)=ln
  for(j in 1:length(ln)){
  	
	x=dff[,ln[j]] #get the column
	x <- as.matrix(x) 
	x1 =(1/x) %*% t(x) #get all  the ratios (using matrix multiplication)
        bound[,j] <- range(x1)
    }
    
  return(bound)
}
restrict<-function(pred,fit,cols)
{
	for(i in cols)
	{
		m=max(fit[,i])
		pos=which(pred[,i] > max)
		if(length(pos)>0)
			pred[pos,i]=max
		
		m=min(fit[,i])
		pos=which(pred[,i] < min)
		if(length(pos)>0)
			pred[pos,i]=min
		
	
	}
}

rep_optim<-function(par1, df, pred.list, bound=NULL)
{
	mem_val=0; mem_val2=0

      while(mem_val==0 | mem_val > mem_val2+.1) #is it the first one, or did we find an improvement last time?
      {
      	op.gamma <- optim(par1, scale_fit_fun, df= df, pred.list= pred.list, bound= bound, ret.dat=F) ##optimisation function
#run single func call to debug
#      	test <- scale_fit_fun(par1, df= df.fit.full2, pred.list= pred.list, extrap= F, ret.dat=F) 

       par1 <- op.gamma$par #extract gamma1
	mem_val=mem_val2
	mem_val2=op.gamma$value
	print(mem_val2)
      }
      return(op.gamma)
}

get_best<-function()
{
	l=list()
	pos=which.min(results_df$AIC)
	AIC=results_df$AIC[pos]
#	pos=18
	print(pos)
	if(is.matrix(params1))
	{
		params=params1[pos,]
	}else{
	
		params=params1[[pos]]
	}
	
	print(params)
	l$params1=params

	params=params[!is.na(params)]
	l$params=params
	l$pred.list=pred.l[[pos]]
	return(l)
}

val_interp<-function()
{
	b<-get_best()

	df=scale_fit_fun(b$params,df.fit.full2, b$pred.list,NULL,T)
	
	print("interpolation")
	if(use_exp==T)
	{
		ss=summary(lm(df$Cost_estimate_per_year_2017_USD_exchange_rate~df$cost_pred))
	}else{
		ss=summary(lm(exp(df$cost_pred)~exp(df$Cost_estimate_per_year_2017_USD_exchange_rate)))
	}
		print(ss)

	if(validate==T)	
	  {
		#re_use df dataframe, since it is the right size. However, change cost_pred to NA, to ensure no erroneous usage of fitted values 
		df$cost_pred=NA 

		for(i in 1:nrow(eu.df1)) #iterate through each datapoint, and fit each row. Return the extrapolated datapoint. 
		{
			print(paste("val point",i))
			pos2=which(df.fit.full2$Index == eu.df1$Index[i]) #remove to fit, put back in to validate
			if(do_bound==T)
			{
				rb= ratio_bound(eu.df1[-i,], pos.pred)	

			}else{
				rb=NULL
			}

			op.gamma=rep_optim(b$params, df.fit.full2[-pos2,],b$pred.list)

		      #use the new best fitting parameters, and get predicted value for the one datapoint
		
     	      		df[i,]=scale_fit_fun(op.gamma$par, rbind(eu.df1,df.fit.full2[pos2[2],]),b$pred.list,rb,T)

		}
		if(use_exp==T)
		{
			ss_val=summary(lm(df$Cost_estimate_per_year_2017_USD_exchange_rate~df$cost_pred))
		}else{
			ss_val=summary(lm(exp(df$Cost_estimate_per_year_2017_USD_exchange_rate)~exp(df$cost_pred)))
		
		}
		r2val=ss_val$r.squared
		
	
		} else {
		r2val=NA
		
	}
	l=list()
	l[[1]]=df

	l[[2]]=c(ss$r.squared,r2val,b$params1[1:6])
	l[[3]]=AIC
	
	### Bootstrapping for interpolation model 
	#  if ((tp == "Damage" & all(c("Pop_size", "Agriculture", "suitability", "time", "time2") %in% pred.list)) |
	#   (tp == "Management" & all(c("GDP", "Agriculture") %in% pred.list))) {
	if (uncertainty) {
		print("Bootstrapping interpolation") #basically, we want to generate an SE for each parameter value here
	  boot_results <<- replicate(1000, {
	    sample_indices <- sample(1:nrow(df), replace = TRUE)
	    sample_indices=c(sample_indices,nrow(df)+sample_indices) #because of the structure of the interpolation algorithm
	    df_sample <- df.fit.full2[sample_indices, ]
	    df_sample$Index=1:nrow(df)
	    # op.gamma <- optim(par1, scale_fit_fun, df= df_sample, pred.list= pred.list, extrap= F, ret.dat=F) ##optimisation function
	    # return(op.gamma$par)
	    return(rep_optim(b$params,df_sample,b$pred.list)$par)
	  } )
	  
	  boot_se <- apply(boot_results,1,sd)
	  l[[4]] <- boot_se
	  
	    saveRDS(boot_results, paste0("inter_boot_results_", tp,  ".rds")) 
	    saveRDS(boot_se, paste0("inter_boot_se_", tp,".rds")) 
	    
	    # Include predictions here
	    pred_interp = boot_results
	    for(k in 1:ncol(pred_interp)) {
	      pred_interp1 = pred_interp[,k ] 
	      rb= ratio_bound(eu.df1, names(pred_interp1)[!names(pred_interp1) %in% c("sd","b0")]  )
	      extrap <- scale_fit_fun(gamma1 = pred_interp1, as.data.frame(df.fit.full1), names(pred_interp1),rb, T)
	      save_name =  file.path(folder_a, paste0("interp_",k,"_",tp,".rds")) 
	      saveRDS(extrap,  save_name)
	    }
	      
	      
	}
	# }
	return(l)
	}

compare_algorithms <-function(interp_pred=NULL, interp_AIC=NULL)
{
	
	#COMPARISON AGAINST OTHER ALGORITHMS
	#create a fitting dataset, removing the validation set, for a fair comparison. This will also only look at non-singlets to make analyses comparable 

	fit_df=eu.df1
	fit_df$t2[is.na(fit_df$t2)]=max(fit_df$t2,na.rm=T)
	#Using GAMS first
	#print("GAM")
	#K=5
	#g=gam(Cost_estimate_per_year_2017_USD_exchange_rate ~ s(GDP,k=K)+s(Pop_size,k=K)+s(Agriculture,k=K)+s(suitability,k=K)+s(t2,k=K),
	#data=fit_df, select = TRUE, method="GCV.Cp")
#	g=gam(Cost_estimate_per_year_2017_USD_exchange_rate ~ s(t2,k=K)+s(suitability,k=K)+s(GDP,k=K)+s(Agriculture,k=K)+s(Pop_size,k=K),
#	data=fit_df, method="REML")
	#sg=summary(g)
#	print(sg) #fitted analyses
	#fit_df$cost_pred=NA
	#for(i in 1:nrow(fit_df)) #iterate through each datapoint, and fit each row. Return the extrapolated datapoint. 
	#{
#		print(paste("val point GAM",i))
	#	g=gam(Cost_estimate_per_year_2017_USD_exchange_rate ~ s(GDP,k=K)+s(Pop_size,k=K)+s(Agriculture,k=K)+s(suitability,k=K)+s(t2,k=K)+comb2,
	#		data=fit_df[-i,], select = TRUE, method="GCV.Cp")
	#		fit_df$cost_pred[i]=predict.gam(g,newdata=fit_df[i,])
	#}	
	#gam_r2val=summary(lm(fit_df$Cost_estimate_per_year_2017_USD_exchange_rate~fit_df$cost_pred))$r.squared
	#print(gam_r2val)	
	#tmp_results[2,5:11]=c(sg$r.sq,gam_r2val, sg$s.table[,3]^.5)


#	print("RF")

	#	fit_df$cost_pred=NA
	#	for(i in 1:nrow(fit_df)) #iterate through each datapoint, and fit each row. Return the extrapolated datapoint. 
	#	{
#		print(paste("val point RF",i))
	#		g=randomForest(Cost_estimate_per_year_2017_USD_exchange_rate ~ GDP+Pop_size+Agriculture+suitability+t2+comb2,
	#		data=fit_df[-i,], proximity=F, ntree= 800, importance = TRUE)
	#		fit_df$cost_pred[i]=predict(g,newdata=fit_df[i,])
	#}	
	#rf_r2val=summary(lm(fit_df$Cost_estimate_per_year_2017_USD_exchange_rate~fit_df$cost_pred))$r.squared
	#print(rf_r2val)	



	#mixed model 
	print("mixed")
	factors=c("GDP","Pop_size","Agriculture","suitability","t2")
	y="Cost_estimate_per_year_2017_USD_exchange_rate ~ (1| comb2)"
	
	f=formula(paste(c(y,factors),collapse="+"))
#	m=lmer(Cost_estimate_per_year_2017_USD_exchange_rate ~ GDP+Pop_size+Agriculture +suitability+t2+ (1|comb2) ,data=fit_df, REML=T)
	m=lmer(f ,data=fit_df, REML=T)
	sm=summary(m)
	glmm_AIC=sm$AICtab
	factors1=factors[which(sm$coefficients[-1,"Pr(>|t|)"]<0.05)]
	f1=formula(paste(c(y,factors1),collapse="+"))
	m=lmer(f1 ,data=fit_df, REML=T)
	print(r.squaredGLMM(m)) #from MuMIn package

	if(validate==T)
	{
		fit_df$cost_pred=NA
		for(i in 1:nrow(fit_df)) #iterate through each datapoint, and fit each row. Return the extrapolated datapoint. 
		{
	#		print(paste("val point GLMM",i))
	#		m_val=lmer(Cost_estimate_per_year_2017_USD_exchange_rate ~ GDP+Pop_size+Agriculture +suitability+t2+ (1|comb2) ,data=fit_df[-i,], REML=T)
			m_val=lmer(f1 ,data=fit_df[-i,], REML=T)
			fit_df$cost_pred[i]=predict(m_val,newdata=fit_df[i,],allow.new.levels=T)
		}	

		if(use_exp==T)
		{
			glmm_r2val=summary(lm(fit_df$Cost_estimate_per_year_2017_USD_exchange_rate~fit_df$cost_pred))$r.squared
		}else{
			glmm_r2val=summary(lm(exp(fit_df$Cost_estimate_per_year_2017_USD_exchange_rate)~exp(fit_df$cost_pred)))$r.squared
		
		}
		print(glmm_r2val)
	}else{
		glmm_r2val=NA
		fit_df$cost_pred=predict(m,newdata=fit_df,allow.new.levels=T)
	}
	tmp_results[2,5:11]=c(r.squaredGLMM(m)[2],glmm_r2val,sm$coefficients[-1,"Pr(>|t|)"])
	
	#		if(!is.null(interp_pred))   
	#	{
	#	
	#		minAIC=min(c(glmm_AIC , interp_AIC ) ) 
	#		glmm_AIC=exp(-.5 *(glmm_AIC-minAIC))
	#		interp_AIC=exp(-.5 *(interp_AIC-minAIC))
	#		ratio=tmp_results[2,5]/(tmp_results[1,5]+tmp_results[2,5])
		#	ratio=glmm_AIC/(glmm_AIC+interp_AIC)     # ceck if works with this 
	#		ensem_cost=(fit_df$cost_pred*ratio+interp_pred*(1-ratio))
	#		tmp_results[3,6]=summary(lm(fit_df$Cost_estimate_per_year_2017_USD_exchange_rate~ensem_cost))$r.squared
	#}
	
	# Calculate SE for GLMM model 
	#if ((tp == "Damage" & all(c("GDP","Pop_size", "Agriculture", "suitability", "time" ) %in% pred.list)) |
	#   (tp == "Management" & all(c("GDP","Pop_size", "Agriculture", "suitability", "time" ) %in% pred.list)) ) {

	#*****

	#****NEED TO GET SE FOR EACH FIXED EFFECT PARAMETER, NOT r.squared
	cat("Measuring SE for glmm")	  
	
	if (uncertainty) {
	  # jsut to get the data
	  pred<- readRDS(paste0("params_main_rm_outl_", tp,  ".rds")) 
	  
	  # In fact for both type of costs the best model is n30 
	  if(tp =="Damage") {
	    #  pre=  brian[brian$type == tp & brian$model == "glmm", ]
	    pred=  pred[30,]
	    factors = paste0(names(pred)[!names(pred) %in% c("sd","b0", "time", "time2")], collapse = "+")
	    f1=formula(paste(c(y,factors),collapse="+"))
	  } else{
	    # pre=  brian[brian$type == tp & brian$model == "glmm", ]
	    pred=  pred[30,]
	    factors = paste0(names(pred)[!names(pred) %in% c("sd","b0", "time", "time2")], collapse = "+")
	    f1=formula(paste(c(y,factors),collapse="+"))
	  }
	  
	  rb= ratio_bound(eu.df1, names(pred)[!names(pred) %in% c("sd","b0")]  )
	  extrap <- scale_fit_fun(gamma1 = pred, as.data.frame(df.fit.full1), names(pred),rb, T)
	  
	  bootstrap_glmm_params <- function(data, formula, extrap, R = 1000, folder_b) {
	    h <- 1  
	    outputs <- vector("list", R)
	    
	    for (i in 1:R) {
	      sample_indices <- sample(1:nrow(data), replace = TRUE)
	      data_sample <- data[sample_indices, ]
	      model_sample <- lmer(formula, data = data_sample, REML = TRUE)
	      
	      pred_glmm <- extrap
	      pred_glmm$t2[is.na(pred_glmm$t2)] <- max(data$t2, na.rm = TRUE)
	      pred_glmm$cost_pred_glmm <- predict(model_sample, newdata = as.data.frame(pred_glmm), allow.new.levels = TRUE)
	      
	      save_name2 = file.path(folder_b, paste0("glmm_", h, "_", tp, ".rds")) # h == i 
	      saveRDS(pred_glmm, save_name2)
	      h <- h + 1  
	      
	      outputs[[i]] <- as.numeric(fixef(model_sample))
	    }
	    
	    param_samples <- do.call(cbind, outputs)
	    
	    param_ses <- apply(param_samples, 1, sd)
	    param_means <- apply(param_samples, 1, mean)
	    
	    # Refit model on the original data for name extraction
	    m2 <- lmer(formula, data = data, REML = TRUE)
	    param_names <- names(fixef(m2))
	    names(param_ses) <- param_names
	    names(param_means) <- param_names
	    
	    list(param_samples = param_samples, param_ses = param_ses, param_means = param_means)
	  }
	  
	  glmm_results <- bootstrap_glmm_params(fit_df, f1, extrap, folder_b = folder_b)
	  

	  # original vs bootstrapped:
	  m=lmer(f1 ,data=fit_df, REML=T) #original model 
	  original <- fixef(m)
    fixed_eff<- data.frame(rbind(as.data.frame(original)),as.data.frame(glmm_results$param_means), as.data.frame(glmm_results$param_ses))
    fixed_eff$percent_error <- (fixed_eff$glmm_results.param_ses / fixed_eff$glmm_results.param_means) * 100
    colnames(fixed_eff) <- c("original_model","mean_boots","se_boots","percent_error")
	  
    saveRDS(glmm_results$param_samples, paste0("glmm_boot_sample_", tp,"_",model_id, ".rds"))
	  saveRDS(fixed_eff, paste0("glmm_se_", tp,"_",model_id, ".rds"))
	}
	#	}
	
	return(tmp_results)	
}


# Functions for extrapolation:

do_df_extrap <- function(df.sp, df.fit){ #this is the main one to check
  #Keep only species/sector/type of cost that are NOT present in the invacost database
  df.sp$comb2 <-  paste(df.sp$Species, df.sp$Type_of_cost_merged) #df.sp is df_extrapolation1
  comb_extrap <- paste(df.sp$Official_country, df.sp$Species,  df.sp$Type_of_cost_merged)
  comb_fit <- paste(df.fit$Official_country, df.fit$Species, df.fit$Type_of_cost_merged)
  
  same <- which(comb_extrap %in% comb_fit)
  df.extrap <- df.sp[- same,] #remove all rows of data that are present in df for country/species/type combo
  #keep only ones that can be extrapolated. Technically, glmm could be applied to everything, but don't worry about that right now.
  #  same2 <- which(df.extrap$comb2 %in% df.fit$comb2) # Note: This remove most of the sp, I only got 11-12
  #  df.extrap <-   df.extrap[same2,] #remove all rows of data that are present in df for country/species/type combo
  
  df.fit$type <- "invacost"
  df.extrap$type <- "extrapol" #to keep track of what is what
  
  df.extrap$Cost_estimate_per_year_2017_USD_exchange_rate <- NA
  df.extrap$cost_mil <- NA
  df.extrap$Cost_ID <- paste("extrap", "_", 1:nrow(df.extrap), sep="")
  
  #reorder columns to use rbind
  df.extrap <- df.extrap[names(df.fit)]
  
  ##We want to keep only the extrapolation part of the dataframe, because the "invacost" df will be the bootstrap
  return(df.extrap)
} 


bound_ratio <- function(v, min.r, max.r){
  v[v < min.r] <- min.r
  v[v > max.r] <- max.r
  return(v)
} 


