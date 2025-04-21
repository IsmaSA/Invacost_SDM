library(raster)
library(sp)


#z = readRDS("glmm_1_Damage.rds")
#z1 = readRDS("inter_boot_se_Management.rds")

#z1 = readRDS("./interpolation_results/interp_1_Damage.rds")
#z2 = readRDS("./glmm_results/glmm_1_Damage.rds")
#nrow(z1)
#nrow(z2)
#b <- readRDS("op_main.rds")
#re = readRDS("op_main.rds")
#results_df=readRDS("result_main_rm_outl_Damage.rds")
#params1=readRDS("params_main_rm_outl_Damage.rds")
library(mgcv)
library(lme4)
library(randomForest)
library(MuMIn)
library(afex)
#rm_time_na=T
rm_outlier=T
use_exp=T
validate=T
#tp="Management"
tp="Damage"
uncertainty=T
do_bound=F
use_exist=F

folder_a <- "interpolation_results"  
folder_b <- "glmm_results"  

op_nm="op_main.rds"
model_id=""
tmp_results=as.data.frame(matrix(NA,nrow=3,ncol=12)) #- get direction and significance.
names(tmp_results)=c("type", "rm_outlier", "rm_time_na", "model", "r2", "r2val","gdp", "pop_size", "agri", "suitability", "time", "time2")
tmp_results[,4]=c("interp","glmm","ensemble")
mem_results=NULL
for(tp in c("Management","Damage"))
#for(tp in c("Damage"))
  #for(tp in c("Management"))
{
print(tp)
for(rm_outlier in c(T))
{
print(paste("remove outlier",rm_outlier))

for(rm_time_na in c(F))
{
print(paste("remove missing time",rm_time_na))

if(tp=="Management")
{
#	outlier_bound=1500
	outlier_bound=600
}else{
	outlier_bound=5000
}
#oname=""
oname="main_"
source("funcsb_IS.R")

if(rm_outlier==T){
  nm_outlier="rm_outl_"
}else{
  nm_outlier=""
}

#Load datasets: 
eu.df <- read_xlsx("eu.df2.xlsx", sheet = "Sheet1")

#eu.df <- read.csv("eu.df2b.csv")
#eu.df$Type_of_cost_merged = "Damage"
eu.df$Impacted_year = as.numeric(eu.df$Impacted_year)
eu.df$t2= eu.df$Impacted_year - eu.df$t
eu.df$t2[eu.df$t2<0]=0
eu.df$GDP=eu.df$GDP/10^9
eu.df$Pop_size=eu.df$Pop_size/10^6
#eu.df$suitability=eu.df$gbif_cells+1
#eu.df$GDP=log(eu.df$GDP)
#eu.df$Pop_size=log(eu.df$Pop_size)
if(rm_time_na==T)
	eu.df=eu.df[! is.na(eu.df$t2),]
eu.df$Cost_estimate_per_year_2017_USD_exchange_rate=exp(eu.df$Cost_estimate_per_year_2017_USD_exchange_rate)/1000000
if(rm_outlier==T){
  eu.df=eu.df[ eu.df$Cost_estimate_per_year_2017_USD_exchange_rate<=outlier_bound,]
  #eu.df=eu.df[ eu.df$cost_bil<=5,] #me
}
if(use_exp==F)
	eu.df$Cost_estimate_per_year_2017_USD_exchange_rate=log(eu.df$Cost_estimate_per_year_2017_USD_exchange_rate)

#for the purpose of evaluation/fitting,  remove all of the singlets, since cannot fit on these. Makes it run faster 
tmp_nms=table(eu.df$comb2)
tmp_nms=tmp_nms[tmp_nms>1] 
single_df=eu.df[!eu.df$comb2 %in% names(tmp_nms),]
single_df$Index=0
single_df$type="single"
eu.df=eu.df[eu.df$comb2 %in% names(tmp_nms),] 
eu.df$Index=1:nrow(eu.df)

##Give all the combinations of predictors (in order to tell scale_fit_fun what it should run)
pos.pred <- c("GDP", "Pop_size", "Agriculture", "suitability","time","time2")

pred.l <- list(m1= c(pos.pred[1]), m2= c(pos.pred[2]), m3= c(pos.pred[3]), m4= c(pos.pred[4]), m5= pos.pred[c(1,2)], m6= pos.pred[c(1,4)], 
               m7= pos.pred[c(1,3)],m8= pos.pred[c(4,2)], m9= pos.pred[c(4,3)], m10= pos.pred[c(2,3)], m11= pos.pred[c(1,3,4)], m12= pos.pred[c(1,2,3)], m13= pos.pred[c(1,2,4)],
               m14= pos.pred[c(2,3,4)], m15= pos.pred[1:4]) #create all combinations of predictors

tm_list=which(unlist(lapply(pred.l,function(x){return("suitability" %in% x)}))) #try suitability both with and without time

for(i in tm_list){
  pred.l[[paste0(i,"t")]]=c(pred.l[[i]],"time")
  pred.l[[paste0(i,"t2")]]=c(pred.l[[i]],"time","time2")
}


full_gamma <-c(0,0,0,0,0,0,NA,.10) #named "full" gamma vector to make the optimization function (scale_fit_fun) faster
names(full_gamma)=c(pos.pred, "sd","b0")
full_gamma2=full_gamma 
full_gamma2[] <-NA



#df_extrapolation1= read_xlsx('df_extrapolation1_version4.0.xlsx') %>% as.data.frame()
#dfsp <- df_extrapolation1 
#dfsp = dfsp[ ,-c(1,4,10,11,13,15,16)]# Cost
#dfsp<- dfsp %>% filter(!Type_of_cost_merged=="Mixed")
#dfsp$t2 = 2024-dfsp$t

#colnames(dfsp)[1] = 'Species'
#colnames(dfsp)[2] = 'Official_country'
#colnames(dfsp)[3]= 'Cost_estimate_per_year_2017_USD_exchange_rate'
#colnames(dfsp)[7] <- "suitability"

#dfsp$GDP <- as.numeric(dfsp$GDP)
#dfsp$Pop_size <- as.numeric(dfsp$Pop_size)
#dfsp$Agriculture <- as.numeric(dfsp$Agriculture)
#dfsp$GDP=dfsp$GDP/10^9
#dfsp$Pop_size=dfsp$Pop_size/10^6

#dfsp$comb2 <- paste(dfsp$Species, dfsp$Type_of_cost_merged)  # ensure this column
#eu.df$comb2 <- paste(eu.df$Species, eu.df$Type_of_cost_merged)

# ensure type of cost
#eu.df2 <- eu.df[eu.df$Type_of_cost_merged==tp, -c(3,6,11,12,16:19)]
#eu.df2$type="invacost"
#dfsp1 <- dfsp[dfsp$Type_of_cost_merged==tp, ] #many things have suitability zero and time NA. Why is any of this needed?

#df.extrap <- do_df_extrap(df.sp= dfsp1, df.fit= eu.df2[, names(dfsp1)])  
#df.extrap$type="extrapol"
#df.fit.full1 <- rbind(eu.df2[, names(df.extrap)], df.extrap) 
#df.fit.full1$Index = 1:nrow(df.fit.full1)

# For extrapolation invavost data:
eu_extra <- read_xlsx("eu.df2.xlsx", sheet = "Sheet1")
eu_extra$Impacted_year = as.numeric(eu_extra$Impacted_year)
eu_extra$t2= eu_extra$Impacted_year - eu_extra$t
eu_extra$t2[eu_extra$t2<0]=0
eu_extra$GDP=eu_extra$GDP/10^9
eu_extra$Pop_size=eu_extra$Pop_size/10^6
eu_extra$Cost_estimate_per_year_2017_USD_exchange_rate=exp(eu_extra$Cost_estimate_per_year_2017_USD_exchange_rate)/1000000
eu_extra=eu_extra[ eu_extra$Cost_estimate_per_year_2017_USD_exchange_rate<=outlier_bound,]
eu_extra$Index=1:nrow(eu_extra)


df_extrapolation1= read_xlsx('df_extrapolation1_version4.0.xlsx') %>% as.data.frame()
dfsp <- df_extrapolation1 
dfsp = dfsp[ ,-c(1,4,10,11,13,15,16)]# Cost
dfsp<- dfsp %>% filter(!Type_of_cost_merged=="Mixed")
dfsp$t2 = 2024-dfsp$t

colnames(dfsp)[1] = 'Species'
colnames(dfsp)[2] = 'Official_country'
colnames(dfsp)[3]= 'Cost_estimate_per_year_2017_USD_exchange_rate'
colnames(dfsp)[7] <- "suitability"

dfsp$GDP <- as.numeric(dfsp$GDP)
dfsp$Pop_size <- as.numeric(dfsp$Pop_size)
dfsp$Agriculture <- as.numeric(dfsp$Agriculture)
dfsp$GDP=dfsp$GDP/10^9
dfsp$Pop_size=dfsp$Pop_size/10^6
dfsp$Index = 1:nrow(dfsp)

dfsp$comb2 <- paste(dfsp$Species, dfsp$Type_of_cost_merged)  # ensure this column
eu.df$comb2 <- paste(eu.df$Species, eu.df$Type_of_cost_merged)
eu_extra$comb2 <- paste(eu_extra$Species, eu_extra$Type_of_cost_merged)

# ensure type of cost
eu_extra1 <- eu_extra[eu_extra$Type_of_cost_merged==tp, -c(3,6,11,12,16:19)]
eu_extra1$type="invacost"
dfsp1 <- dfsp[dfsp$Type_of_cost_merged==tp, ] #many things have suitability zero and time NA. Why is any of this needed?

#df.extrap <- do_df_extrap(df.sp= dfsp1, df.fit= eu.df1)  
df.extrap <- do_df_extrap(df.sp= dfsp1, df.fit= eu_extra1[, names(dfsp1)])  
eu_extra1 <- eu_extra1[, -c(5)] # to remove
eu_extra1 <- eu_extra1[ -10]
df.full.extrap <- rbind(eu_extra1, df.extrap[,names(eu_extra1)])  # to remove

absent <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]), unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"])) #those are in the invacost dataset but NOT in the extrapol
absent2 <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"]), unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]))
df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent),]
df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent2),]
df.fit.full1 <- df.full.extrap
df.fit.full1$Index = 1:nrow(df.fit.full1)
print(unique(df.fit.full1$Species))




results_df <- data.frame()

#eu.df$Cost_estimate_per_year_2017_USD_exchange_rate = eu.df$cost_bil * 1000
#This code include the information for b0 and b1 
df.fit1 <- eu.df

eu.df$type <- "invacost"
df.fit1$type <- "extrapol"
df.fit.full <- rbind(eu.df, df.fit1)
df.fit.full <- as.data.frame(df.fit.full)
lst1=list()
params2=list()
i = tp #{
  
  lst<-list()
  params<-list()

eu.df1<-  eu.df[eu.df$Type_of_cost_merged==i,]
df.fit.full2 <- df.fit.full[df.fit.full$Type_of_cost_merged == i,]

if(use_exist==F) 
{
	  source("loopsb.R")
}else{
	results_df=readRDS(paste0("params_",oname,nm_outlier,i,".rds"))
	params1=readRDS(paste0("params_",oname,nm_outlier,i,".rds"))
	
}
  
#}
#check r^2 as well
print(paste("type",tp))
print(paste("remove outlier",rm_outlier))
print(paste("remove unknown times",rm_time_na))
tmp_results[,1]=tp
tmp_results[,2]=rm_outlier
tmp_results[,3]=rm_time_na
tmp_interp=val_interp()
tmp_results[1,5:12]=tmp_interp[[2]]
tmp_results=compare_algorithms(tmp_interp[[1]]$cost_pred,tmp_interp[[3]])
mem_results=rbind(mem_results,tmp_results)

}}}
#saveRDS(mem_results,op_nm)
## Damage
#dama <- readRDS("result_test2_rm_outl_Damage.rds")
#dama<- dama %>% filter(Model==29)

#dama1 <- readRDS("params_test2_rm_outl_Damage.rds")
#dama1 <- dama1[29,]

### Managemente
#mana <- readRDS("result_test2_rm_outl_Management.rds")
#mana<- mana %>% filter(Model==18)

#mana1 <- readRDS("params_test2_rm_outl_Management.rds")
#mana1<- mana1[18,]
