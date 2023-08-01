

### Code for Brian: 

# Ismael Soto 
# 13/07/2023

# Load functions: 
do_df_extrap <- function(df.sp, df.fit){
  #Keep only species/sector/type of cost that are NOT present in the invacost database
  df.sp$comb2 <-  paste(df.sp$Species, df.sp$Type_of_cost_merged)
  df.sp$comb <- paste(df.sp$Official_country, df.sp$Species,  df.sp$Type_of_cost_merged)
  df.fit$comb <- paste(df.fit$Official_country, df.fit$Species, df.fit$Type_of_cost_merged)
  
  same <- which(df.sp$comb %in% df.fit$comb)
  df.extrap <- df.sp[-same,] #remove all rows of data that are present in df
  
  df.fit$type <- "invacost" #to keep track of waht is what
  df.extrap$type <- "extrapol" #to keep track of what is what
  
  df.extrap$Cost_estimate_per_year_2017_USD_exchange_rate <- NA
  df.extrap$cost_mil <- NA
  df.extrap$Cost_ID <- paste("extrap", "_", 1:nrow(df.extrap), sep="")
  
  #reorder columns to use rbind
  df.extrap <- df.extrap[names(df.fit)]
  
  ##We want to keep only the extrapolation part of the dataframe, because the "invacost" df will be the bootstrap
  return(df.extrap)
} 
get_ratio <- function(dff, pred.list){
  
  ln <- pred.list #replace pred.list by ln (list.name) so that we don't overwrite it
  
  pos <- which(dff$type == "invacost")
  x <- dff[pos,] #invacost data
  y <- dff[-pos,] #extrapolation data
  
  
  mxr=matrix(NA, nrow=nrow(y), ncol=length(ln)*2)
  
  for(i in 1:nrow(y)){
    if(is.null(ln)){
      return(mxr)
      
    } else {
      
      for(j in 1:length(ln)){
        mxr[i, j] <- min(y[i, ln[j]]/x[,ln[j]])
        mxr[i, j+length(ln)] <- max(y[i, ln[j]]/x[,ln[j]])
      }
    }
  }
  mxr2 <- as.data.frame(mxr)
  colnames(mxr2) <- c(paste(rep("min", length(pred.list)), pred.list, sep="_"), paste(rep("max", length(pred.list)), pred.list, sep="_"))
  return(mxr2)
} 
bound_ratio <- function(v, min.r, max.r){
  v[v < min.r] <- min.r
  v[v > max.r] <- max.r
  return(v)
} 
scale_fit_fun <- function(gamma1, df, pred.list, extrap = F, fit.test =F){ #gamma1= vector of parameters to optimize
  
  ln <- pred.list #replace pred.list by ln (list.name) so that we don't overwrite it
  
  if(extrap == F){ #default (extrap = F) is when we want to estimate the gamma parameters (with optim)
    if(length(gamma1) == 1){
      names(gamma1) <- "sd"
    }
    full_gamma[names(gamma1)]=gamma1
    
    scale_fun <- function(dff, fit.test){
      pos <- which(dff$type == "invacost")
      x <- dff[pos,] #invacost data
      y <- dff[-pos,] #extrapolation data
      
      
      if(nrow(y) == 1){
        y$cost_pred <- NA 
        
        
        if(fit.test ==T){ #return the full dataset to calculate metrics of fit
          return(y)
        } else {
          ll <- sum(log(dnorm(y$cost_mil - y$cost_pred, sd= full_gamma["sd"]))) #calculate log-likelihood sd= gamma1[length(gamma1)]
          return(-ll) #return -ll because optim finds the minimum
        }
        
      } else {
        
        for(i in 1:nrow(y)){
          
          pred1 <- x$cost_mil[-i] * (y$GDP[i] / x$GDP[-i])^full_gamma["GDP"] * 
            (y$Pop_size[i] / x$Pop_size[-i])^full_gamma["Pop_size"] * 
            (y$Agriculture[i] / x$Agriculture[-i])^full_gamma["Agriculture"] * 
            (y$suitability[i] / x$suitability[-i])^full_gamma["suitability"]
          y$cost_pred[i] <- mean(pred1)
        }
        
        if(fit.test ==T){
          return(y)
          
        } else {
          ll <- sum(log(dnorm(y$cost_mil - y$cost_pred, sd= full_gamma["sd"])))
          return(-ll) 
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
      return(measure.fit)
    }
    
  } else { #When extrap = T we want to obtain an estimation of costs
    
    scale_fun <- function(dff){
      
      ln <- pred.list
      
      pos <- which(dff$type == "invacost")
      x <- dff[pos,] #invacost data
      y <- dff[-pos,] #extrapolation data
      
      
      v <- list() #v will contains the ratios
      for(i in 1:nrow(y)){
        if(is.null(ln)){
          pred1 <- x$cost_mil
          y$cost_pred[i] <- mean(pred1)
        } else {
          
          
          for(j in 1:length(ln)){
            v[[j]] <- bound_ratio(y[i, ln[j]]/x[,ln[j]], min.r= min(ratio.fit[,paste("min", pred.list[j], sep="_")]), max.r= max(ratio.fit[,paste("max", pred.list[j], sep="_")])) #use pred.list instead of ln because the name for sector is "sector" not the specific one
          }
          
          vm <- matrix(unlist(v), ncol=length(ln)) #create a matrix where each column is a predictor (so that we can multiply across rows)
          colnames(vm) <- ln
          
          prod.vm <-apply(vm, MARGIN= 1, function(x){ prod(x^gamma1[-length(gamma1)])}) #multiply each row and apply the gamma1 exp
          
          pred1 <- x$cost_mil * prod.vm
          y$cost_pred[i] <- mean(pred1)
          
        }
      }
      
      return(y) #Here we want to keep the entire dataframe to compare/analyse the results
    }
    
    df.by <- by(df, df$comb2, scale_fun)
    df.full <- do.call(rbind, df.by) 
    
    return(df.full)
  }
} ## seems that now is working
aic.fit <- function(x, gamma1){
  LL <- sum(log(dnorm(x$cost_mil - x$cost_pred, sd= gamma1[length(gamma1)])), na.rm=T)
  aic <- -2*(LL) + 2*(length(gamma1)-1) #-1 because one of the parameters is the standard deviation
  return(aic)
} 


# Load some packages:
suppressMessages({
  library(dplyr, quiet = TRUE, warn.conflicts = FALSE)
  library(reshape, quiet = TRUE, warn.conflicts = FALSE)
  library(ggplot2)
  library(doBy)
  library(tidyr)
  library(dplyr)
  library(tidyr)  
  library(stringr)
  library(invacost)
  library(biomod2)
  library(ggplot2)
  library(rgbif)
  library(gridExtra)
  library(raster)
  library(fasterize)
  library(sf)
  library(rasterVis)
  library(sp)
  library(mgcv)
  library(caret)
  library(pROC)
  library(xlsx)
  library(readxl)
  library(writexl)
  library(rgdal)
  Sys.setenv(LANGUAGE = "en")
})



# InvaCost data:
invacost33= read_xlsx('invacost33.xlsx')
eu.df.full <- invacost33
eu.df.full$comb2 <- paste(eu.df.full$species, eu.df.full$Type_of_cost_merged) #create a column with the combination of species name, sector impacted, and type of costs (damage, management, mixed)
eu.df.full <- eu.df.full[order(eu.df.full$comb2),] #reorder so everything matches later

# Transformation (!)
eu.df.full$cost_mil <- log1p(eu.df.full$cost_mil)


# Extrapolation data
df_extrapolation1= read_xlsx('df_extrapolation1.xlsx')
dfsp <- df_extrapolation1 #data for extrapolation
dfsp = dfsp[ ,-3]




##Give all the combinations of predictors (in order to tell scale_fit_fun what it should run)
pos.pred <- c("GDP", "Pop_size", "Agriculture", "suitability")
pred.l <- list(m1= c(pos.pred[1]), m2= c(pos.pred[2]), m3= c(pos.pred[3]), m4= c(pos.pred[4]), m5= pos.pred[c(1,2)], m6= pos.pred[c(1,4)], m7= pos.pred[c(1,3)], 
               m8= pos.pred[c(4,2)], m9= pos.pred[c(4,3)], m10= pos.pred[c(2,3)], m11= pos.pred[c(1,3,4)], m12= pos.pred[c(1,2,3)], m13= pos.pred[c(1,2,4)],
               m14= pos.pred[c(2,3,4)], m15= pos.pred) #create all combinations of predictors
full_gamma <-c(0,0,0,0,NA) #named "full" gamma vector to make the optimization function (scale_fit_fun) faster
names(full_gamma)=c(pos.pred, "sd")




#Change names to the databases
eu.df <- eu.df.full
colnames(eu.df)[1] = 'Species'
colnames(eu.df)[2] = 'Official_country'
colnames(eu.df)[4]= 'Cost_estimate_per_year_2017_USD_exchange_rate'

colnames(dfsp)[1] = 'Species'
colnames(dfsp)[2] = 'Official_country'

cost.types = unique(invacost33$Type_of_cost_merged)
v.AIC.costtypes <- numeric(length(cost.types))
names(v.AIC.costtypes) = cost.types

likelihoods <- numeric()

results <- data.frame(
  model = character(), 
  aic_difference = numeric(), 
  likelihood_difference = numeric(), 
  Best_model = character(),
  stringsAsFactors = FALSE
)


pred.list <- pred.l[['m7']]  #Example

for(m in 1:length(pred.l)){ #loop over each model
  pred.list <- pred.l[[m]] #to obtain the correct predictors
  
  par1 <- c(rep(0.1, length(pred.list)), sd(eu.df$cost_mil, na.rm=T)) #starting values for optimization function
  
  names(par1) <- c(pred.list, "sd")
  par1
  
  df.extrap <- do_df_extrap(df.sp= dfsp, df.fit= eu.df) #this contains ONLY the extrapolation data, no invacost data
  
  df.fit1 <- eu.df
  eu.df$type <- "invacost"
  df.fit1$type <- "extrapol"
  df.fit.full <- rbind(eu.df, df.fit1)
  
  df.fit.full <- as.data.frame(df.fit.full)
  
  
  
  for (i in cost.types) {
    df.fit.full2 <- df.fit.full[df.fit.full$Type_of_cost_merged == i,]
    
    if(length(par1) == 1){ #cannot use optim with a single precictor
      par1 <- c(0,1000) #reassign par1 beacause optim need to know the min and max to look for the optimized value
      op.gamma <- optimize(scale_fit_fun, par1, maximum = F, df= df.fit.full2, pred.list= pred.list, extrap= F, fit.test=F) #optimisation function
      gamma1 <- op.gamma$minimum #extract gamma1
      
    } else { #when more than a single parameter to fit
      op.gamma <- optim(par1, scale_fit_fun, df= df.fit.full2, pred.list= pred.list, extrap= F, fit.test=F) ##optimisation function
      gamma1 <- op.gamma$par #extract gamma1
    }
    
    names(gamma1) <- c(pred.list, "sd")
    gamma1
    
    
    
    ##Metrics of fit --> calculate AIC for the given set of parameters gamma
    fit.t1 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, extrap = F, fit.test=T) #extract df to calculate measure of fit
    
    v.AIC <- aic.fit(fit.t1, gamma1= gamma1)
    v.AIC.costtypes[i] <- v.AIC #store the AIC
    
    # log-likelihood:
    fit.t11 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, extrap = F, fit.test=F) 
    likelihoods[i] <- fit.t11
    
    
    ##Estimation of projected costs
    ratio.fit <- get_ratio(df.fit.full2, pred.list= pred.list) #obtain all the ratios
    df.extrap <- df.extrap[names(eu.df)] #make sure the columns names match
    
    
    # Adapt the eu.df dataframe based on the type of cost:
    if (i =="Damage") {
      eu.df1<-  eu.df %>% filter(Type_of_cost_merged=="Damage")
      df.full.extrap <- rbind(eu.df1, df.extrap) 
    } else if (i =="Management") {
      eu.df1<-  eu.df %>% filter(Type_of_cost_merged=="Management")
      df.full.extrap <- rbind(eu.df1, df.extrap) 
    } else if (i =="Mixed") {
      eu.df1<-  eu.df %>% filter(Type_of_cost_merged=="Mixed")
      df.full.extrap <- rbind(eu.df1, df.extrap) 
    }
    
    
    
    
    #We can't estimate the cost for a certain comb2 that is NOT present in the invacost df
    #Need to remove comb2 that are present in extrapol but absent in invacost and vice versa
    absent <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]), unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"])) #those are in the invacost dataset but NOT in the extrapol
    absent2 <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"]), unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]))
    df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent),]
    df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent2),]
    
    
    df.full.extrap = as.data.frame(df.full.extrap)
    
    
    
    #Run the function scale_fit_fun to obtain predictde costs
    extrap <- scale_fit_fun(gamma1 = gamma1, df.full.extrap, pred.list=pred.list, extrap = T, fit.test= F)
    
    # Transform back the cost (back-transforming):
    extrap$original_scale = expm1(extrap$cost_pred)
    
   # write.csv2(extrap, file=paste0("Ismael_SDM", "_",m,"_",i, ".csv"))
    sum_pred_cost <- sum(extrap$cost_pred, na.rm= T)
    
    gamma.df2 <- c(gamma1, v.AIC, sum_pred_cost, i)
    names(gamma.df2) <- c(pred.list, "sd", "v.AIC", "sum_pred_cost", "type_of_cost")
    
   # write.csv2(gamma.df2, paste0("./", "Ismael_SDM_pred", "_", names(pred.l[m]),"_",i, ".csv"))
    cat( pred.list,'Done  ===', "\n")
  }
  
  
  # sum of the likelihood for the model m across type of costs
  total_likelihood <- sum(likelihoods)
  total_likelihood 
  
  
  
  
  #### Run the overall model   #### 
  
  # Fit the full model after the model loop
  par1_full <- c(rep(0.1, length(pred.list)), sd(eu.df$cost_mil, na.rm=T)) #starting values for optimization function
  names(par1_full) <- c(pred.list, "sd")
  par1_full
  
  
  
  if(length(par1) == 1){ #cannot use optim with a single precictor
    par1 <- c(0,1000) #reassign par1 beacause optim need to know the min and max to look for the optimized value
    op.gamma_full <- optimize(scale_fit_fun, par1_full, maximum = F, df= df.fit.full, pred.list= pred.list, extrap= F, fit.test=F) #optimisation function
    gamma1_full <- op.gamma_full$minimum #extract gamma1
    
  } else { #when more than a single parameter to fit
    op.gamma_full <- optim(par1_full, scale_fit_fun, df= df.fit.full, pred.list= pred.list, extrap= F, fit.test=F) ##optimisation function
    gamma1_full <- op.gamma_full$par #extract gamma1
  }
  names(gamma1_full) <- c(pred.list, "sd")
  gamma1_full  
  
  
  
  fit.t1_full <- scale_fit_fun(gamma1 = gamma1_full, df.fit.full, pred.list= pred.list, extrap = F, fit.test=T) #extract df to calculate measure of fit
  
  v.AIC_full <- aic.fit(fit.t1_full, gamma1= gamma1_full)
  v.AIC_full
  
  
  
  # log-likelihood:
  fit.t11_full <- scale_fit_fun(gamma1 = gamma1_full, df.fit.full, pred.list= pred.list, extrap = F, fit.test=F) 
  likelihood_full <- fit.t11_full
  
  
  
  ratio.fit <- get_ratio(df.fit.full, pred.list= pred.list) #obtain all the ratios
  df.extrap <- df.extrap[names(eu.df)] #make sure the columns names match
  df.full.extrap <- rbind(eu.df, df.extrap)
  
  
  #We can't estimate the cost for a certain comb2 that is NOT present in the invacost df
  #Need to remove comb2 that are present in extrapol but absent in invacost and vice versa
  absent <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]), unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"])) #those are in the invacost dataset but NOT in the extrapol
  absent2 <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"]), unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]))
  df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent),]
  df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent2),]
  
  
  df.full.extrap = as.data.frame(df.full.extrap)
  
  #Run the function scale_fit_fun to obtain predictde costs
  extrap <- scale_fit_fun(gamma1 = gamma1, df.full.extrap, pred.list=pred.list, extrap = T, fit.test= F)
  
  
  # Transform back the cost (back-transforming):
  extrap$original_scale = expm1(extrap$cost_pred)
  
  
  #write.csv2(extrap, file=paste0("Ismael_SDM", "_",m,"_", ".csv"))
  sum_pred_cost <- sum(extrap$cost_pred, na.rm= T)
  
  
  gamma.df3 <- c(gamma1_full, v.AIC, sum_pred_cost)
  names(gamma.df3) <- c(pred.list, "sd", "v.AIC", "sum_pred_cost")
  gamma.df3
  
  #write.csv2(gamma.df3, paste0("./", "Ismael_SDM_pred", "_", names(pred.l[m]),"_", ".csv"))
  cat( pred.list,'Overall  ===', "\n")
  
  
  ## Check AIC partial models (Damage, Management & Mixed costs) vs. full model
  
  aic_difference <- sum(v.AIC.costtypes) - v.AIC_full
  aic_difference 
  
  
  likelihood_difference <- total_likelihood - likelihood_full
  likelihood_difference
  
  results <- rbind(results, data.frame(
    model = m,  
    aic_difference = aic_difference,
    likelihood_difference = likelihood_difference,
    best_model = if (likelihood_difference > 0) "Separate type of costs" else "full model",
    stringsAsFactors = FALSE
  ))
  
}

####################################################################################

####################################################################################

####################################################################################


cost.types = unique(invacost33$Type_of_cost_merged)
v.AIC.costtypes <- numeric(length(cost.types))
names(v.AIC.costtypes) = cost.types

results_df <- data.frame()




for(m in 1:length(pred.l)){ #loop over each model
  pred.list <- pred.l[[m]] #to obtain the correct predictors
  
  par1 <- c(rep(0.1, length(pred.list)), sd(eu.df$cost_mil, na.rm=T)) #starting values for optimization function
  
  names(par1) <- c(pred.list, "sd")
  par1
  
  df.extrap <- do_df_extrap(df.sp= dfsp, df.fit= eu.df) #this contains ONLY the extrapolation data, no invacost data
  
  df.fit1 <- eu.df
  eu.df$type <- "invacost"
  df.fit1$type <- "extrapol"
  df.fit.full <- rbind(eu.df, df.fit1)
  
  df.fit.full <- as.data.frame(df.fit.full)
  
  
  
  for (i in cost.types) {
    df.fit.full2 <- df.fit.full[df.fit.full$Type_of_cost_merged == i,]
    
    if(length(par1) == 1){ #cannot use optim with a single precictor
      par1 <- c(0,1000) #reassign par1 beacause optim need to know the min and max to look for the optimized value
      op.gamma <- optimize(scale_fit_fun, par1, maximum = F, df= df.fit.full2, pred.list= pred.list, extrap= F, fit.test=F) #optimisation function
      gamma1 <- op.gamma$minimum #extract gamma1
      
    } else { #when more than a single parameter to fit
      op.gamma <- optim(par1, scale_fit_fun, df= df.fit.full2, pred.list= pred.list, extrap= F, fit.test=F) ##optimisation function
      gamma1 <- op.gamma$par #extract gamma1
    }
    
    names(gamma1) <- c(pred.list, "sd")
    gamma1
    
    
    
    ##Metrics of fit --> calculate AIC for the given set of parameters gamma
    fit.t1 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, extrap = F, fit.test=T) #extract df to calculate measure of fit
    
    v.AIC <- aic.fit(fit.t1, gamma1= gamma1)
    v.AIC.costtypes[i] <- v.AIC #store the AIC
    
    # log-likelihood:
    fit.t11 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, extrap = F, fit.test=F) 
    likelihoods[i] <- fit.t11
    
    
    ##Estimation of projected costs
    ratio.fit <- get_ratio(df.fit.full2, pred.list= pred.list) #obtain all the ratios
    df.extrap <- df.extrap[names(eu.df)] #make sure the columns names match
    
    
    # Adapt the eu.df dataframe based on the type of cost:
    if (i =="Damage") {
      eu.df1<-  eu.df %>% filter(Type_of_cost_merged=="Damage")
      df.full.extrap <- rbind(eu.df1, df.extrap) 
    } else if (i =="Management") {
      eu.df1<-  eu.df %>% filter(Type_of_cost_merged=="Management")
      df.full.extrap <- rbind(eu.df1, df.extrap) 
    } else if (i =="Mixed") {
      eu.df1<-  eu.df %>% filter(Type_of_cost_merged=="Mixed")
      df.full.extrap <- rbind(eu.df1, df.extrap) 
    }
    
    
    
    
    #We can't estimate the cost for a certain comb2 that is NOT present in the invacost df
    #Need to remove comb2 that are present in extrapol but absent in invacost and vice versa
    absent <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]), unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"])) #those are in the invacost dataset but NOT in the extrapol
    absent2 <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"]), unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]))
    df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent),]
    df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent2),]
    
    
    df.full.extrap = as.data.frame(df.full.extrap)
    
    
    
    #Run the function scale_fit_fun to obtain predictde costs
    extrap <- scale_fit_fun(gamma1 = gamma1, df.full.extrap, pred.list=pred.list, extrap = T, fit.test= F)
    
    # Transform back the cost (back-transforming):
    extrap$original_scale = expm1(extrap$cost_pred)
    
    # write.csv2(extrap, file=paste0("Ismael_SDM", "_",m,"_",i, ".csv"))
    sum_pred_cost <- sum(extrap$cost_pred, na.rm= T)
    
    gamma.df2 <- c(gamma1, v.AIC, sum_pred_cost, i)
    names(gamma.df2) <- c(pred.list, "sd", "v.AIC", "sum_pred_cost", "type_of_cost")
    
    # write.csv2(gamma.df2, paste0("./", "Ismael_SDM_pred", "_", names(pred.l[m]),"_",i, ".csv"))
    new_results_df <- data.frame(
      "AIC" = sum(v.AIC.costtypes),
      "Likelihood" = sum(likelihoods),
      "Model" = m,
      "Full/Partial_model" = "Partial"
    )
    
    cat( pred.list,'Done  ===', "\n")
  
    }
  results_df <- rbind(results_df, new_results_df)
}




##########################      FULL MODEL     ##########################

##########################      FULL MODEL     ##########################


results_df1 <- data.frame()


for(m in 1:length(pred.l)){ #loop over each model
  pred.list <- pred.l[[m]] #to obtain the correct predictors
  
  par1 <- c(rep(0.1, length(pred.list)), sd(eu.df$cost_mil, na.rm=T)) #starting values for optimization function
  
  names(par1) <- c(pred.list, "sd")
  par1
  
  df.extrap <- do_df_extrap(df.sp= dfsp, df.fit= eu.df) #this contains ONLY the extrapolation data, no invacost data
  
  df.fit1 <- eu.df
  eu.df$type <- "invacost"
  df.fit1$type <- "extrapol"
  df.fit.full <- rbind(eu.df, df.fit1)
  
  df.fit.full <- as.data.frame(df.fit.full)
  
  
    
    if(length(par1) == 1){ #cannot use optim with a single precictor
      par1 <- c(0,1000) #reassign par1 beacause optim need to know the min and max to look for the optimized value
      op.gamma <- optimize(scale_fit_fun, par1, maximum = F, df= df.fit.full, pred.list= pred.list, extrap= F, fit.test=F) #optimisation function
      gamma1 <- op.gamma$minimum #extract gamma1
      
    } else { #when more than a single parameter to fit
      op.gamma <- optim(par1, scale_fit_fun, df= df.fit.full, pred.list= pred.list, extrap= F, fit.test=F) ##optimisation function
      gamma1 <- op.gamma$par #extract gamma1
    }
    
    names(gamma1) <- c(pred.list, "sd")
    gamma1
    
    
    
    ##Metrics of fit --> calculate AIC for the given set of parameters gamma
    fit.t1 <- scale_fit_fun(gamma1 = gamma1, df.fit.full, pred.list= pred.list, extrap = F, fit.test=T) #extract df to calculate measure of fit
    
    v.AIC <- aic.fit(fit.t1, gamma1= gamma1)
    v.AIC.costtypes[i] <- v.AIC #store the AIC
    
    # log-likelihood:
    fit.t11 <- scale_fit_fun(gamma1 = gamma1, df.fit.full, pred.list= pred.list, extrap = F, fit.test=F) 
    likelihoods[i] <- fit.t11
    
    
    ##Estimation of projected costs
    ratio.fit <- get_ratio(df.fit.full, pred.list= pred.list) #obtain all the ratios
    df.extrap <- df.extrap[names(eu.df)] #make sure the columns names match
    
    df.full.extrap <- rbind(eu.df, df.extrap) 
    
    
    
    #We can't estimate the cost for a certain comb2 that is NOT present in the invacost df
    #Need to remove comb2 that are present in extrapol but absent in invacost and vice versa
    absent <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]), unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"])) #those are in the invacost dataset but NOT in the extrapol
    absent2 <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"]), unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]))
    df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent),]
    df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent2),]
    
    
    df.full.extrap = as.data.frame(df.full.extrap)
    
    
    
    #Run the function scale_fit_fun to obtain predictde costs
    extrap <- scale_fit_fun(gamma1 = gamma1, df.full.extrap, pred.list=pred.list, extrap = T, fit.test= F)
    
    # Transform back the cost (back-transforming):
    extrap$original_scale = expm1(extrap$cost_pred)
    
    # write.csv2(extrap, file=paste0("Ismael_SDM", "_",m,"_",i, ".csv"))
    sum_pred_cost <- sum(extrap$cost_pred, na.rm= T)
    
    gamma.df2 <- c(gamma1, v.AIC, sum_pred_cost, i)
    names(gamma.df2) <- c(pred.list, "sd", "v.AIC", "sum_pred_cost", "type_of_cost")
    
    # write.csv2(gamma.df2, paste0("./", "Ismael_SDM_pred", "_", names(pred.l[m]),"_",i, ".csv"))
    new_results_df1 <- data.frame(
      "AIC" = sum(v.AIC.costtypes),
      "Likelihood" = sum(likelihoods),
      "Model" = m,
      "Full/Partial_model" = "Full"
    )
    
    results_df1 <- rbind(results_df1, new_results_df1)
    cat( pred.list,'Done  ===', "\n")
  }
  


z <- bind_rows(results_df, results_df1)



####################################################################################

####################################################################################

####################################################################################

cost.types = unique(invacost33$Type_of_cost_merged)
v.AIC.costtypes <- numeric(length(cost.types))
names(v.AIC.costtypes) = cost.types

results_df <- data.frame()



for(m in 1:length(pred.l)){ #loop over each model
  pred.list <- pred.l[[m]] #to obtain the correct predictors
  
  
  df.extrap <- do_df_extrap(df.sp= dfsp, df.fit= eu.df) #this contains ONLY the extrapolation data, no invacost data
  
  df.fit1 <- eu.df
  eu.df$type <- "invacost"
  df.fit1$type <- "extrapol"
  df.fit.full <- rbind(eu.df, df.fit1)
  
  df.fit.full <- as.data.frame(df.fit.full)
  
  for (i in cost.types) {
    
    if (i =="Damage") {
      eu.df1<-  eu.df %>% filter(Type_of_cost_merged=="Damage")
     
    } else if (i =="Management") {
      eu.df1<-  eu.df %>% filter(Type_of_cost_merged=="Management")
   
    } else if (i =="Mixed") {
      eu.df1<-  eu.df %>% filter(Type_of_cost_merged=="Mixed")

    }
    
    
    par1 <- c(rep(0.1, length(pred.list)), sd(eu.df1$cost_mil, na.rm=T)) #starting values for optimization function
    
    names(par1) <- c(pred.list, "sd")
    par1
    
    
    
    df.fit.full2 <- df.fit.full[df.fit.full$Type_of_cost_merged == i,]
    
    if(length(par1) == 1){ #cannot use optim with a single precictor
      par1 <- c(0,1000) #reassign par1 beacause optim need to know the min and max to look for the optimized value
      op.gamma <- optimize(scale_fit_fun, par1, maximum = F, df= df.fit.full2, pred.list= pred.list, extrap= F, fit.test=F) #optimisation function
      gamma1 <- op.gamma$minimum #extract gamma1
      
    } else { #when more than a single parameter to fit
      op.gamma <- optim(par1, scale_fit_fun, df= df.fit.full2, pred.list= pred.list, extrap= F, fit.test=F) ##optimisation function
      gamma1 <- op.gamma$par #extract gamma1
    }
    
    names(gamma1) <- c(pred.list, "sd")
    gamma1
    
    
    
    ##Metrics of fit --> calculate AIC for the given set of parameters gamma
    fit.t1 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, extrap = F, fit.test=T) #extract df to calculate measure of fit
    
    v.AIC <- aic.fit(fit.t1, gamma1= gamma1)
    #v.AIC.costtypes[i] <- v.AIC #store the AIC
    
    # log-likelihood:
    fit.t11 <- scale_fit_fun(gamma1 = gamma1, df.fit.full2, pred.list= pred.list, extrap = F, fit.test=F) 
    likelihoods <- fit.t11
    
    
    ##Estimation of projected costs
    ratio.fit <- get_ratio(df.fit.full2, pred.list= pred.list) #obtain all the ratios
    df.extrap <- df.extrap[names(eu.df)] #make sure the columns names match
    
    
    # Adapt the eu.df dataframe based on the type of cost:
    if (i =="Damage") {
      df.full.extrap <- rbind(eu.df1, df.extrap) 
    } else if (i =="Management") {
      df.full.extrap <- rbind(eu.df1, df.extrap) 
    } else if (i =="Mixed") {
      df.full.extrap <- rbind(eu.df1, df.extrap) 
    }
    
    
    
    
    #We can't estimate the cost for a certain comb2 that is NOT present in the invacost df
    #Need to remove comb2 that are present in extrapol but absent in invacost and vice versa
    absent <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]), unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"])) #those are in the invacost dataset but NOT in the extrapol
    absent2 <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"]), unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]))
    df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent),]
    df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent2),]
    
    
    df.full.extrap = as.data.frame(df.full.extrap)
    
    
    
    #Run the function scale_fit_fun to obtain predictde costs
    extrap <- scale_fit_fun(gamma1 = gamma1, df.full.extrap, pred.list=pred.list, extrap = T, fit.test= F)
    
    # Transform back the cost (back-transforming):
    extrap$original_scale = expm1(extrap$cost_pred)
    
    # write.csv2(extrap, file=paste0("Ismael_SDM", "_",m,"_",i, ".csv"))
    sum_pred_cost <- sum(extrap$cost_pred, na.rm= T)
    
    gamma.df2 <- c(gamma1, v.AIC, sum_pred_cost, i)
    names(gamma.df2) <- c(pred.list, "sd", "v.AIC", "sum_pred_cost", "type_of_cost")
    
    # write.csv2(gamma.df2, paste0("./", "Ismael_SDM_pred", "_", names(pred.l[m]),"_",i, ".csv"))
    new_results_df <- data.frame(
      "AIC" = v.AIC    ,
      "Likelihood" = likelihoods ,
      "Type_cost" = i, 
      "Model" = m                ,
      "Full/Partial_model" = "Partial"
    )
    results_df <- rbind(results_df, new_results_df)
    cat( pred.list,'Done  ===', "\n")
    
  }
  
}











