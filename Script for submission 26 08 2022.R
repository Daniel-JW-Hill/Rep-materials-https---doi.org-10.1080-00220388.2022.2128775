#Technology heterogeneity and poverty traps: A latent class approach to technology gap drivers of chronic poverty
#
#This script accompanies the paper Technology heterogeneity and poverty traps: A latent class approach to technology gap drivers of chronic poverty
#Published in the Journal of Development Studies

#Prepared by Daniel Hill
#University of New England (Australia)
#August 2022
#dhill41@une.edu.au
#hill88756@gmail.com


rm(list=ls())

library(readxl)     #reads data
library(plm)        #For panel data regressions
library(fmsb)       #for radar diagrams
library(radarchart) #for radar diagrams
library(stargazer)  #to print tables

#############################################
#  RETRIEVE DATA
#############################################

#set Working Directory
# Set relevant working directory here with data

load()

############################################
#  Define required functions
#############################################

#Generate translog function 
gen_translog=function(x,timevar=NULL,tech_type=c("full","hicks-neutral","interactions only"),time_inc_type=c("time interact","time-dummy","linear","quad"),logged=TRUE,x_names=NULL){
  #Generates a translog x matrix from a matrix on (possible logged) input levels
  #Allows for different types of technology as specified in 'tech-type'
  #Allows for inclusion of a time variable and specification of its effect on technology via 'time_inc_type'
  
  #logged argument refers to whether the x matrix supplied is already in logs (TRUE) or not (FALSE)
  
  #TECH TYPE:
  #'full' is the full set of second order effects assumed in a translog including interactions and quadratic terms
  #'hicks-neutral' generates a matrix of first order and quadratic effects only - no interaction effects.
  #'interactions only' generates a matrix of first order and interaction effects only - no quadratic effects.
  
  #TIME TYPE:
  #all of these generate additional columns appended to the generated matrix according to:
  #'time interact' time interactions on all input variables included - explicitly allows for changing usage of inputs over time (changing technology)
  #'time-dummy' specifies a dummy for each time period included in the timevar vector. Implies non-parametric shifts in the PPC over time.
  #'linear' imposes a linear trend on the PPC over time
  #'quad' imposes a quadratic trend on the PPC over time
  
  #The optional x_names argument
  
  #Get lnx if not already. return warning if unlogged x matrix has non-positive inputs
  if(logged==FALSE){
    if(sum(x<=0)>0){return(warning("input matrix has non-positive inputs which are not allowed in a translog due to the strict essentiality condition"))}
    lnx=log(x)
  } else {
    lnx=x
  }
  
  quadeffects=0.5*lnx^2
  for(i in 1:(ncol(lnx)-1)){
    for(j in (i+1):(ncol(lnx))){
      if(i==1 && j==2){
        inteffects=0.5*x[,i]*x[,(j)]
      } else {
        inteffects=cbind(inteffects,0.5*x[,i]*x[,(j)])
      }
    }
  }
  
  if(tech_type=="full"){
    lnx=cbind(lnx,quadeffects,inteffects)
  } else if(tech_type=="hicks-neutral"){
    lnx=cbind(lnx,quadeffects)
  } else {
    lnx=cbind(lnx,inteffects)
  }
  
  #get time effects if timevar included
  if(is.null(timevar)==FALSE){
    if(time_inc_type=="time interact"){
      tx=timevar*lnx
    } else if(time_inc_type=="time-dummy"){
      tx=c()
      for(i in 2:length(levels(as.factor(timevar)))){
        tx[,i]=1*timevar==levels(as.factor(timevar))[i]
      }
    } else if(time_inc_type=="linear"){
      tx=timevar
    } else {
      tx=cbind(timevar,timevar^2)
    }
    lnx=cbind(tx,lnx)
  }
  
  if(is.null(x_names)==FALSE){
    xnames=colnames(x)
    quadnames=paste(xnames,"2",sep="")
    j=ncol(x)
    intnames=c()
    for(i in 1:(j-1)){
      if(i==1){
        intnames=paste(xnames[1],xnames[-1],sep="x")
      } else {
        intnames=c(intnames,paste(xnames[i],xnames[-1:-i],sep="x"))
      }
    }
    if(is.null(timevar)==FALSE){
      if(time_inc_type=="time-dummy"){
        xnames=c(paste("T",seq(1,max(timevar),1),sep=""),xnames)
      }
      if(time_inc_type=="linear"){
        xnames=c("T",xnames)
      }
      if(time_inc_type=="quad"){
        xnames=c("T","T2",xnames)
      }
    }
    if(tech_type=="full"){
      colnames(lnx)=c(xnames,quadnames,intnames)
    } else {
      if(tech_type=="interactions only"){
        colnames(lnx)=c(xnames,intnames)
      } else {
        colnames(lnx)=c(xnames,quadnames)
      }
    }
  } 
  return(lnx)
}


#use inverse hyperbolic sin to allow non-negative values and approximate a log.
invhypsin=function(x){
  out=log(x+(x^2+1)^(1/2))
  return(out)}


#############################################
#  Initial analysis
#############################################
#Prices are already corrected for CPI using 2009. this is the CPI series.
cpi=cbind(seq(2001,2014,1),
         c(1.5046,
           1.4550, 1.3798, 1.2958, 1.2414,
           1.1776, 1.1241, 1.0376, 1.0000,
           0.9128, 0.8380, 0.7600, 0.6838,
           0.6410)
         )

YEARS = levels(as.factor(data$YEAR))
#correct for cpi - only monetary variables
for(year in 1:length(YEARS)){
    cpi_factor=cpi[cpi[,1]==YEARS[year],2]
    data[data$YEAR==as.numeric(YEARS[year]),7:15]=data[data$YEAR==as.numeric(YEARS[year]),7:15]*cpi_factor
    }

#get requirement index and do asset regression
TIME = data$YEAR - min(data$YEAR) + 1
REQUIREMENT = data$HHSIZE*500*9.7790 #This is Naschold's (2012) requirements index as 500Rp in 1975 prices
LIVELIHOODS = (data$INCOME_TOTAL)/REQUIREMENT
#8.4% of observations across the entire sample are 'poor' on this measure.
#8.4% in extreme poverty in 2001
#5% in 2014. 
#Many falling into extreme poverty in 2005-06
# Badiani, R., Dercon, S. and Krishnan, P., 2007. Changes in Living Standards in Villages in India 1975-2004: Revisiting the ICRISAT village level studies. Chronic Poverty Research Centre Working Paper, (85).


#Poverty line calculations
# $1.90 USD a day (PPP 2011) is the accepted poverty line
#PPP conversion factor (rupee to USD) in 2009 (World Bank data). 
ppp_factor_2009 = 13.36286066 
ppp_factor_2011 = 15.5495491
Inc_USD_PPP = data$INCOME_TOTAL/(ppp_factor_2009*365) #USD per day in 2009 prices. 
Inc_USD_PPP = Inc_USD_PPP*(ppp_factor_2011/ppp_factor_2009) /data$HHSIZE #change to 2011 prices and scale for HH size

# 17# of observations are below the poverty line. 
# 27% were under the poverty line in 2001
# 9% in 2014 

#Summary of income per person per day (USD ppp)
summary(Inc_USD_PPP)

assets = data.frame(data[,1:2],
                    TIME=TIME,
                    INCOME_PC=(data$INCOME_TOTAL)/data$HHSIZE/10000,
                    LIVELIHOOD=LIVELIHOODS,
                    LIVESTOCK_PC=(data$SML_LIVESTOCK+data$LRG_LIVESTOCK)/data$HHSIZE/10000,
                    SML_LIVESTOCK_PC=data$SML_LIVESTOCK/data$HHSIZE/10000,
                    LRG_LIVESTOCK_PC=data$LRG_LIVESTOCK/data$HHSIZE/10000,
                    BUILDINGS_PC=data$BUILDINGS/data$HHSIZE/10000,
                    DRYLAND_PC = data$DRYLAND/data$HHSIZE,
                    IRRIGATED_PC = data$IRRIGATED/data$HHSIZE,
                    MOBILE_CAPITAL_PC=(data$DURABLES+data$FARM_EQUIP+data$STOCK)/data$HHSIZE/10000,
                    SAVINGS_PC=data$SAVINGS/data$HHSIZE/10000,
                    EDUC_PC = data$EDUC_TOT/data$HHSIZE ,
                    HHSIZE = data$HHSIZE
                    )

#asset data summary
ass_summ=matrix(NA,ncol=7,nrow=ncol(assets[,-1:-3]))
for(r in 1:nrow(ass_summ)){
  ass_summ[r,1] = min((assets[,-1:-3][,r]),na.rm=TRUE)
  ass_summ[r,2] = quantile((assets[,-1:-3][,r]),0.25,na.rm=TRUE)
  ass_summ[r,3] = median((assets[,-1:-3][,r]),na.rm=TRUE)
  ass_summ[r,4] = quantile((assets[,-1:-3][,r]),0.75,na.rm=TRUE)
  ass_summ[r,5] = max((assets[,-1:-3][,r]),na.rm=TRUE)
  ass_summ[r,6] = mean((assets[,-1:-3][,r]),na.rm=TRUE)
  ass_summ[r,7] = sqrt(var((assets[,-1:-3][,r]),na.rm=TRUE))
  }

rownames(ass_summ) = colnames(assets[,-1:-3])
colnames(ass_summ) = c("min","Q25","median","Q75","max","mean","std.dev")
#write.table(ass_summ,"results\\summary of assets data.csv",sep=",")

#generate translog version of data:
data_TL=gen_translog(assets[,-1:-5],tech_type="full",x_names=names(assets[,-1:-5]))
data_TL=data.frame(assets[,1:5],data_TL)


#test for functional form in both periods - expect Cobb-Douglas is fine under FE both ways treatment.
#GET r-squared for twoways model:
eq_CD=(LIVELIHOODS)~as.factor(HHID)+as.factor(TIME)+TIME*SML_LIVESTOCK_PC+TIME*LRG_LIVESTOCK_PC+TIME*BUILDINGS_PC+TIME*DRYLAND_PC+
                    TIME*IRRIGATED_PC+TIME*MOBILE_CAPITAL_PC+TIME*SAVINGS_PC+TIME*EDUC_PC+HHSIZE

regCD=lm(eq_CD,data=data_TL)

#calculate asset index:
IDs = levels(as.factor(data$HHID))
pred_asset = predict(regCD)

data_pred=data.frame("HHID"=data_TL$HHID,
                     "TIME"=data_TL$TIME,
                     "PRED"=pred_asset)


#############################################
#  CONVERGENCE GRAPH
#############################################

#convergence graph first and second epochs 
pred_asset_true=pred_asset
live_mat = data.frame("HHID"=assets$HHID,"TIME"=assets$TIME,"LIVE"=pred_asset_true)
t=as.numeric(levels(as.factor(TIME)))
dat_=list()
for(i in 1:length(t)){
  dat_[[i]]=density(subset((live_mat$LIVE),TIME==t[i]),na.rm=TRUE)
  }
fac <- 5.2  # A factor to make the densities overlap
plot(1, type = "n", xlim = c(-5,25), ylim = c(1, 14 + 3),
     axes=FALSE,xlab="", ylab = "")
axis(side=1)

YEARS=as.numeric(levels(as.factor(data$YEAR)))
# Add each density, shifted by i and scaled by fac
for(i in 1:14){
  lines(  dat_[[i]]$x, fac*dat_[[i]]$y + i)
  polygon(dat_[[i]]$x, fac*dat_[[i]]$y + i, col = rgb(0, 0, 0, 0.4), border = NA)
  abline(h = i, lwd = 0.5)
  text(-0.4,i+0.5,YEARS[i],cex=0.8)
}
#livelihoods densities across time indicate increasing inequality


#############################################
#  CLUSTERING TO INITIALISE LATENT GROUPING
#############################################

#load required libraries
require(mclust)

#latent class analysis
n_classes = 7 #number of groups tested
n_polynomials = 5 #max polynomial level tested

classvec=seq(1,n_classes)
polyvec=seq(1,n_polynomials)
#first get priors for class membership - need a single vector of observations for each individual
#    so that allocate to a single cohort only
#    Baseline regressions show that IRRIGATED, BUILDINGS, MOBILE_CAPITAL, EDUC, LRG AND SMALL LIVESTOCK AND HHSIZE ARE ALL IMPORTANT.
#    But need to limit the number of variables for effective clustering.
#    Use weighted average (Weights as CD coefs) across three epochs in the 2001-2014 period.


IDs=levels(as.factor(assets$HHID))
mat=matrix(NA,nrow=length(IDs),ncol=6) #ncol = number of variables for clustering 
coefs=regCD$coef[-1:-227]
coefs=coefs[c(1,2,4,5,6,7,8)] #retains the variables not interacted with time
for(i in 1:length(IDs)){
  dat1 = subset(assets,assets$HHID==IDs[i])
  ASS1 = subset(pred_asset,assets$HHID==IDs[i]) #finds average of fitted values from cobb-douglas asset index
  ass1 = mean(subset(ASS1,dat1$YEAR>2000 & dat1$YEAR<2005))
  ass2 = mean(subset(ASS1,dat1$YEAR>2004 & dat1$YEAR<2010))
  ass3 = mean(subset(ASS1,dat1$YEAR>2009))
  dataTL1 = subset(cbind(data_TL$SML_LIVESTOCK_PC,data_TL$LRG_LIVESTOCK_PC,data_TL$DRYLAND_PC,data_TL$IRRIGATED_PC,data_TL$MOBILE_CAPITAL_PC,data_TL$SAVINGS_PC,data_TL$EDUC_PC),data_TL$HHID==IDs[i])
  INP=dataTL1%*%coefs  #finds the weighted average of asset levels using coefficients as the weights
  inp1 = mean(subset(INP,dat1$YEAR>2000 & dat1$YEAR<2005))
  inp2 = mean(subset(INP,dat1$YEAR>2004 & dat1$YEAR<2010))
  inp3 = mean(subset(INP,dat1$YEAR>2009))
  mat[i,]=c(ass1,ass2,ass3,inp1,inp2,inp3)
}

classmat=matrix(NA,nrow=length(IDs),ncol=length(classvec)+1)
classmat[,1]=IDs
for(i in 1:length(classvec)){
  if(i==1){
    classmat[,(i+1)]=rep(1,nrow(classmat))
    next
    }
  classmat[,(i+1)]=Mclust(mat,G=i)$classification
}

classmat1=sapply(data.frame(classmat[,-1]), function(x) as.numeric(as.character(x)))#turns group cohorts from character strings to recognisable elements
classmat=data.frame(classmat[,1],classmat1) #returns HH ids

#############################################
#  LATENT GROUPING ESTIMATION
#############################################

#first get a full length vector of classes - duplicates household observations across all observations 
classes=matrix(NA,nrow=nrow(data),ncol=length(classvec))
for(r in 1:nrow(data)){
  classes[r,]=c(1,unlist(classmat[classmat[,1]==assets$HHID[r],-1:-2]))
  }

HHID=assets$HHID
class1=c()
classes1=classmat[,-1]
classmat1=data.frame("ID"=data$HHID,classes)
modlist=reglist=list()

#Retrieves demeaned data
y=data_TL$LIVELIHOOD
eq_cd_basic=y~SML_LIVESTOCK_PC+LRG_LIVESTOCK_PC+BUILDINGS_PC+DRYLAND_PC+IRRIGATED_PC+MOBILE_CAPITAL_PC+SAVINGS_PC+EDUC_PC+HHSIZE
xmat_cd=model.matrix(eq_cd_basic,data=data_TL)[,-1]
xmat_for_means=xmat_demeaned=cbind(y,xmat_cd)
means=matrix(NA,nrow=length(IDs),ncol=ncol(xmat_for_means))
for(rr in 1:length(IDs)){
  for(cc in 1:ncol(xmat_for_means)){
    means[rr,cc]=mean(xmat_for_means[data_TL$HHID==IDs[rr],cc])
    xmat_demeaned[data_TL$HHID==IDs[rr],cc]=xmat_for_means[data_TL$HHID==IDs[rr],cc]-means[rr,cc]
    }
  }

data_demeaned=data.frame(HHID=data$HHID,TIME=TIME,xmat_demeaned)
eqCD_demeaned=y~as.factor(TIME)+SML_LIVESTOCK_PC+LRG_LIVESTOCK_PC+BUILDINGS_PC+DRYLAND_PC+IRRIGATED_PC+MOBILE_CAPITAL_PC+SAVINGS_PC+EDUC_PC+HHSIZE
regCD_demeaned=lm(eqCD_demeaned,data=data_demeaned)


#build demeaned matrices for loop
yd=xmat_demeaned[,1] #demeaned livelihoods
xmat_demeaned=model.matrix(regCD_demeaned) #remove the single TIME series put in there by R's silly lm function:
#xmat_demeaned=xmat_demeaned[,-15]


nmin=15 #minimum in a group
for(cc in 1:length(classvec)){
  cl=classes1[,cc]
  repeat{
    loglID=matrix(NA,nrow=length(IDs),ncol=max(cl))
    logls=matrix(NA,nrow=nrow(data_TL),ncol=max(cl))
    for(gg in 1:max(cl)){
      reglist[[gg]]=out=lm(yd~0+xmat_demeaned,subset=classmat1[,(cc+1)]==gg) #already has intercept...
      sigma_reg=summary(out)$sigma
      preds=xmat_demeaned%*%out$coef  #replace the NAs from the time interactions
      errors=yd-preds
      logls[,gg]=dnorm(errors,mean=0,sd=sigma_reg,log=TRUE)
      }
      for(i in 1:length(IDs)){
        loglID[i,]=apply(subset(logls,HHID==IDs[i]),2,sum,na.rm=TRUE)
        }
      pmat=exp(loglID)/apply(exp(loglID),1,sum)

    for(rr in 1:nrow(loglID)){
      class1[rr]=which(pmat[rr,]==max(pmat[rr,]))
      }
    modlist[[cc]]=reglist

    #check numbers are at least nmin and reallocate:
    counts=c()
    for(gg in 1:cc){
      counts[gg]=sum(class1==gg)
      }
    if(sum(counts<nmin)>0){
    #there is at least one group with less than nmin households
      permlist=list()
      for(gg in 1:max(class1)){
        if(counts[gg]>=nmin){next}
        perm=order(pmat[,gg],decreasing=TRUE)
        #now ensure no swapping
        class1_perm=class1[perm]
        for(i in 1:length(IDs)){
          #repeat gathering households into this group without taking from a group with not enough
          #start with most probable households
          if(class1_perm[i]==gg){next}    #don't take someone from their own group
          if(sum(class1==class1_perm[i])<=nmin){next}   #don't take someone from another small group
          if(0.2<runif(1,0,1)){next}  #random selection component to avoid a never-ending loop
          class1[perm[i]]=gg   #else change the group
          if((sum(class1)==gg)==nmin){break}  #if numbers are now enough break from this level and check remaining groups
          }
        }
      }

    #model has converged - no difference in class allocations between iterations
    if(sum(class1==classes1[,cc])==nrow(loglID)){break}

    #model has not converged - difference in class allocations, reset to new allocation.
    classes1[,cc]=class1
    for(i in 1:nrow(data)){
      classmat1[i,(cc+1)]=class1[classmat1[i,1]==IDs]
      }
    }
  }


#####################################################
#  LATENT GROUPING POST ESTIMATION ANALYSIS (Table 2)
#####################################################

#counts the final numbers in each group 
countsmat=matrix(NA,nrow=max(classvec),ncol=max(classvec))
for(r in 1:n_classes){
  for(c in 1:n_classes){
    countsmat[r,c]=sum(classes1[,c]==r)
  }
}

#Tests the number of clusters to determine best  based on aic, bic and log likelihood
len=length(modlist)
k=ncol(xmat_demeaned)
logl=bic=aic=rep(0,len)
#do chi-squared tests for models as well (they are nested) - need to do this
for(cc in 1:len){
  logl_=0
  for(gg in 1:cc){
    logl_=logLik(modlist[[cc]][[gg]])+logl_
    }
  logl[cc]=logl_
  aic[cc]=-2*logl_+2*k*cc
  bic[cc]=-2*logl_+log(nrow(xmat_demeaned))*k*cc
  }
models_summary=rbind(k,logl,aic,bic)
rownames(models_summary)=c("K","Log-likelihood","AIC","BIC")
write.table(models_summary,"model summaries at class level.csv",sep=",")

#manually pick cluster number - should be based on minimum aic
model_index=4  #check if still the minimum from line above
model=modlist[[model_index]]
classvec1=classes1[,model_index]

#############################################
#  ESTABLISH ORDERING OF GROUPS FIRST
#############################################

#For each person get their actual predicted value with their mean value
pred_vals=rep(NA,nrow=nrow(data_TL))
for(i in 1:length(IDs)){
  meanval=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
  group=classvec1[classmat[,1]==IDs[i]]
  xmat1=xmat_demeaned[data_TL$HHID==IDs[i],]
  pred=xmat1%*%model[[group]]$coef
  pred_vals[data_TL$HHID==IDs[i]]=meanval+pred
  }

boxplot(pred_vals ~ classmat1[,(model_index+1)])
#order is 3, 1, 2, 4 in ascending

#re-label groups so they are ascending
classvec = classmat1[,model_index+1]
classvec[classvec==3]=0
classvec[classvec==1]=1
classvec[classvec==2]=2
classvec[classvec==4]=3
data$COHORT = classvec

#check:
boxplot(pred_vals ~ classvec,outline = FALSE)

##########################################################
#  TABLE 1 - SUMMARY STATISTICS OF ASSET LEVELS (BY GROUP)
##########################################################

mobile_capital = data$DURABLES+data$FARM_EQUIP+data$STOCK
datasummary_group = cbind(data, mobile_capital)

#Subset data for each group
datasummary_group1 = subset(datasummary_group, data$COHORT==0)
datasummary_group2 = subset(datasummary_group, data$COHORT==1)
datasummary_group3 = subset(datasummary_group, data$COHORT==2)
datasummary_group4 = subset(datasummary_group, data$COHORT==3)

#Summary statistics table (mean and standard error for each group)
assetsummary_groups=matrix(NA,ncol=10,nrow=ncol(datasummary_group[,-1:-5]))
for(r in 1:nrow(assetsummary_groups)){
  assetsummary_groups[r,1] = mean((datasummary_group[,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,2] = sqrt(var((datasummary_group [,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_group ))
  assetsummary_groups[r,3] = mean((datasummary_group1[,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,4] = sqrt(var((datasummary_group1[,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_group1))
  assetsummary_groups[r,5] = mean((datasummary_group2[,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,6] = sqrt(var((datasummary_group2[,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_group2))
  assetsummary_groups[r,7] = mean((datasummary_group3[,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,8] = sqrt(var((datasummary_group3[,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_group3))
  assetsummary_groups[r,9] = mean((datasummary_group4[,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,10] = sqrt(var((datasummary_group4[,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_group4))
}
rownames(assetsummary_groups) = colnames(datasummary_group[,-1:-5])
colnames(assetsummary_groups) = c("mean","std.err","mean","std.err","mean","std.err","mean","std.err","mean","std.err")
assetsummary_groups = assetsummary_groups[!is.na(assetsummary_groups[1,]),] #removes NA rows where variable is a character string and not numeric (e.g. Caste, occupation etc)
write.table(assetsummary_groups, file = "summarystats_byGroup.csv", sep=",")

#T-Tests for significance of asset levels for each group compared to the mean
assetsummary_ttest=matrix(NA,ncol=4,nrow=ncol(datasummary_group[,-1:-5]))
for(r in 1:nrow(assetsummary_ttest)){
  if (class(datasummary_group[,-1:-5][,r])=="character"){
    assetsummary_ttest[r,] = NA
  } else {
  assetsummary_ttest[r,1] = t.test(datasummary_group[,-1:-5][,r], datasummary_group1[,-1:-5][,r])$p.value
  assetsummary_ttest[r,2] = t.test(datasummary_group[,-1:-5][,r], datasummary_group2[,-1:-5][,r])$p.value
  assetsummary_ttest[r,3] = t.test(datasummary_group[,-1:-5][,r], datasummary_group3[,-1:-5][,r])$p.value
  assetsummary_ttest[r,4] = t.test(datasummary_group[,-1:-5][,r], datasummary_group4[,-1:-5][,r])$p.value
  }
}
rownames(assetsummary_ttest) = colnames(datasummary_group[,-1:-5])
colnames(assetsummary_ttest) = c("G1", "G2", "G3", "G4")
assetsummary_ttest = assetsummary_groups[!is.na(assetsummary_ttest[1,]),] #removes NA rows where variable is a character string and not numeric (e.g. Caste, occupation etc)
write.table(assetsummary_ttest, file = "summarystats_byGroupttests.csv", sep=",")


#############################################
#  ORDINAL PROBIT AND INVERSE MILLS RATIO TO RE-ESTIMATE GROUP-LEVEL TECHNOLOGIES
#############################################

#equation for ordinal logit:
require(ordinal)

eq_ord = as.factor(COHORT) ~ CASTE_CATEG + HH_HEAD_AGE + HH_HEAD_EDUCYR + DEPENDENCY_RAT + MAIN_OCCUPATION + HH_CATEGORY_ORIG + HHSIZE + I(EDUC_TOT-HH_HEAD_EDUCYR) + CHILDREN 
reg_ord = clm(eq_ord,data=data,link="probit")
anova(reg_ord)

coefs = reg_ord$coef
intercepts = matrix(c(-Inf, coefs[1],coefs[1],coefs[2],coefs[2],coefs[3],coefs[3], Inf),ncol=2,byrow=TRUE)
colnames(intercepts) = c("Lower Threshold", "Upper Threshold")
rownames(intercepts) = c("Cohort 1", "Cohort 2", "Cohort 3", "Cohort 4")

xmat_ord = model.matrix(~ CASTE_CATEG + HH_HEAD_AGE + HH_HEAD_EDUCYR + DEPENDENCY_RAT + MAIN_OCCUPATION + HH_CATEGORY_ORIG + HHSIZE + I(EDUC_TOT-HH_HEAD_EDUCYR) + CHILDREN  ,data=data)[,-1]
pred_ordinal = xmat_ord%*%coefs[-1:-3]

#Percent-correctly predicted
PCP_mat = cbind(data$COHORT, pred_ordinal, rep(NA, length(pred_ordinal)), rep(NA, length(pred_ordinal)))
for (i in 1:length(data$COHORT)){
  PCP_mat[i,3] = ifelse(PCP_mat[i,2]<=intercepts[1,2], 0, ifelse(PCP_mat[i,2]>intercepts[1,2]&&PCP_mat[i,2]<=intercepts[2,2], 1, ifelse(PCP_mat[i,2]>intercepts[3,2], 3, 2)))
  PCP_mat[i,4] = ifelse(PCP_mat[i,1]==PCP_mat[i,3], 1, 0)
}
PCP_ordinal = sum(PCP_mat[,4])/length(PCP_mat[,4]) 

#Pseudo R2 (McFadden R2 and McFadden adjusted R2)
eq_ord_int = as.factor(COHORT)~ 1
reg_ord_int = clm(eq_ord_int,data=data,link="probit")
PseudoR2 = 1-(reg_ord$logLik/reg_ord_int$logLik)
adj_PseudoR2 = 1-((reg_ord$logLik-(reg_ord$nobs-reg_ord$df.residual))/reg_ord_int$logLik)

#Likelihood ratio test for weak instruments
chi2 = (2* log(reg_ord$logLik-(reg_ord_int$logLik)) - qchisq(0.99,1))

## Function to calculate mills ratio from ordinal probit
## From the outline in Chiburis, R. and Lokshin, M. 2007. Maximum likelihood and two-step estimation of
##  an ordered-probit selection model. The Stata Journal 7(2): 167-182

millsCalcOrdinal = function (pred, int0, int1){
  numerator = (dnorm(int0 - pred) - dnorm(int1 - pred))
  denominator = (pnorm(int1 - pred) - pnorm(int0 - pred))
  mills = numerator/denominator
  return(mills)
  }

millsvec = rep(0,nrow(data))
for(i in 1:nrow(data)){
  cohort = data$COHORT[i]
  millsvec[i] = millsCalcOrdinal(pred_ordinal[i],intercepts[(data$COHORT[i]+1),1],intercepts[(data$COHORT[i]+1),2])
  }

###############################################
#  RE-ESTIMATE MODELS WITH MILLS VEC INCLUDED - without region effects. 
###############################################
model1 = lm(yd~0+xmat_demeaned + millsvec ,subset=data$COHORT == 0)
model2 = lm(yd~0+xmat_demeaned + millsvec ,subset=data$COHORT == 1)
model3 = lm(yd~0+xmat_demeaned + millsvec ,subset=data$COHORT == 2)
model4 = lm(yd~0+xmat_demeaned + millsvec ,subset=data$COHORT == 3)

modelsList = list(model1, model2, model3, model4)
model0 = modlist[[1]][[1]] #Homogeneous technology estimation

#demeaned lm models report the within R2. Below we estimate the overall R2 by adding the levels back in. The predicted values are also needed for the trajectory functions. 

#For each person get their actual predicted value with their mean value - this replaces the pred_vals calculated without the inverse mills ratio. 

xmat_demeaned_IMR = model.matrix(yd~0+xmat_demeaned + millsvec)
pred_vals=rep(NA,nrow=nrow(data_TL))
for(i in 1:length(IDs)){
  meanval=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
  group = as.numeric(levels(as.factor(data$COHORT[data_TL$HHID==IDs[i]])))
  group=group+1
  xmat1=xmat_demeaned_IMR[data_TL$HHID==IDs[i],]
  pred=xmat1%*%modelsList[[group]]$coef
  pred_vals[data_TL$HHID==IDs[i]]=meanval+pred
}

OverallR2_vec = as.data.frame(cbind(pred_vals, data_TL$LIVELIHOOD, data$COHORT))
colnames(OverallR2_vec) = c("pred_vals", "LIVELIHOOD", "COHORT")

#repeat this for the full sample
pred_vals0=rep(NA,nrow=nrow(data_TL))
for(i in 1:length(IDs)){
  meanval0=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
  xmat10=xmat_demeaned[data_TL$HHID==IDs[i],]
  pred0=xmat10%*%model0$coef
  pred_vals0[data_TL$HHID==IDs[i]]=meanval0+pred0
}

OverallR2 = Adj_OverallR2 = rep(NA, 5)

#First find Overall R2 for whole group
OverallR2[1] = cor(pred_vals0, data_TL$LIVELIHOOD)^2
N = length(model0$fitted.values)   #number of observations
K = ncol(xmat_demeaned)+1   #number of coefficients plus inverse mills ratio
Adj_OverallR2[1] = 1-( ((1-OverallR2[1]) * (N-1))/(N-K-1) )


#Now for group technologies
for (g in 1:model_index){
  R2_vec = subset(OverallR2_vec,OverallR2_vec$COHORT==g-1)
  OverallR2[g+1] = cor(R2_vec$pred_vals, R2_vec$LIVELIHOOD)^2
  N = length(modelsList[[g]]$fitted.values) #number of observations
  K = ncol(xmat_demeaned)+1  #number of coefficients plus inverse mills ratio
  Adj_OverallR2[g+1] = 1-( ((1-OverallR2[g+1]) * (N-1))/(N-K-1) )
}

###############################################
#  PAIRWISE DIFFERENCES IN TECHNOLOGIES (Z-SCORES)
###############################################

#test pairwise differences in coefficients between the groups using z-scores

modelsList2 = modelsList 
modelsList2[[5]] = model0


#Function to compare beta coefficients between groups using the z scores
z.score = function(beta1, beta2, se1, se2){
  z.score = (beta1 - beta2)/(sqrt(se1^2 +se2^2))
  return(z.score)
}


Coeftests = list() #stores pairwise comparisons in a list, each element in the list is a different coefficient
for (l in 1:length(modelsList2[[5]]$coefficients)){ #ignore inverse mills ratio as not relevant for whole group 
  coeftestmat = matrix(NA, ncol = 5, nrow = 5)
  rownames(coeftestmat) = c(colnames(xmat_demeaned)[l],1,2,3,4)
  colnames(coeftestmat) = c("WS",1,2,3,4)
  for (r in 1:5){
    beta1 = modelsList2[[r]]$coefficients[[l]] 
    se1 = coef(summary(modelsList2[[r]]))[l, "Std. Error"]
    for (c in 1:5){
      beta2 = modelsList2[[c]]$coefficients[[l]] 
      se2 = coef(summary(modelsList2[[c]]))[l, "Std. Error"]
      coeftestmat[r,c] = z.score(beta1, beta2, se1, se2)
    }
  }
  Coeftests[[l]] = coeftestmat
}

###############################################
#  RE_ESTIMATE MODELS WITH Household and village dummies instead of demeaned data. 
# This has been done post-hoc following reviewer comments. 
###############################################

# generate region dummies instead of village dummies.
Region = data$VILLAGE
Region[Region=="AUREPALLE"]= "Mahbubnagor"
Region[Region=="DOKUR"]= "Mahbubnagor"
Region[Region=="SHIRAPUR"]= "Solapur"
Region[Region=="KALMAN"]= "Solapur"
Region[Region=="KANZARA"]= "Akola"
Region[Region=="KINKHED"]= "Akola"

data_TL_Vil = cbind(data_TL, Region, data$COHORT)

eq_V_FE = LIVELIHOOD ~ as.factor(TIME) + as.factor(Region) + (as.factor(TIME)*as.factor(Region)) + SML_LIVESTOCK_PC + LRG_LIVESTOCK_PC + BUILDINGS_PC + DRYLAND_PC +
  IRRIGATED_PC + MOBILE_CAPITAL_PC + SAVINGS_PC + EDUC_PC + HHSIZE + millsvec + as.factor(HHID)

model1 = lm(eq_V_FE, data=data_TL_Vil, subset=data_TL_Vil$`data$COHORT` == 0)
model2 = lm(eq_V_FE, data=data_TL_Vil, subset=data_TL_Vil$`data$COHORT` == 1)
model3 = lm(eq_V_FE, data=data_TL_Vil, subset=data_TL_Vil$`data$COHORT` == 2)
model4 = lm(eq_V_FE, data=data_TL_Vil, subset=data_TL_Vil$`data$COHORT` == 3)

modelsList = list(model1, model2, model3, model4) #Turn off if want version without region fixed effects. 
model0 = lm(eq_V_FE, data=data_TL_Vil) #Homogeneous technology estimation

#predicted values with village fixed effects included. 
predict1 = predict(model1,  data=data_TL_Vil, subset=data_TL_Vil$`data$COHORT` == 0)
predict1 = cbind(subset(data_TL_Vil,data_TL_Vil$`data$COHORT` == 0), predict1)
colnames(predict1[73])=c("predict")
predict2 = predict(model2,  data=data_TL_Vil, subset=data_TL_Vil$`data$COHORT` == 1)
predict2 = cbind(subset(data_TL_Vil,data_TL_Vil$`data$COHORT` == 1), predict2)
colnames(predict2[73])=c("predict")
predict3 = predict(model3,  data=data_TL_Vil, subset=data_TL_Vil$`data$COHORT` == 2)
predict3 = cbind(subset(data_TL_Vil,data_TL_Vil$`data$COHORT` == 2), predict3)
colnames(predict3[73])=c("predict")
predict4 = predict(model4,  data=data_TL_Vil, subset=data_TL_Vil$`data$COHORT` == 3)
predict4 = cbind(subset(data_TL_Vil,data_TL_Vil$`data$COHORT` == 3), predict4)
colnames(predict4[73])=c("predict")

pred_list = list(predict1, predict2, predict3, predict4)

pred_vals=rep(NA,nrow=nrow(data_TL))

for(i in 1:nrow(data_TL)){
  hh = data_TL_Vil$HHID[i]
  y = data_TL_Vil$TIME[i]
  cc = data_TL_Vil$`data$COHORT`[i]
  data_predicted = pred_list[[cc+1]]
  data_predicted = subset(data_predicted, data_predicted$HHID == hh)
  data_predicted = subset(data_predicted, data_predicted$TIME == y)
  pred_vals[i] = data_predicted$predict
}
  
###############################################
#  LIVELIHOOD DYNAMICS - retrieve predicted livelihoods
###############################################

#First retrieve the predicted livelihood values and lagged livelihood values from the estimated livelihood technologies. 

pred_vals = pred_vals #already calculated for the overall R2. This reintroduces mean levels from demeaned data

#Retrieve lagged predicted livelihoods from predicted values. 
current=pred_vals
lagged=rep(NA,length(current))
IDvec=data_TL$HHID
for(i in 2:length(current)){
  if(IDvec[i]==IDvec[i-1]){
    lagged[i]=current[i-1]
    } else {
    lagged[i]=NA
  }
}

#Save predicted livelihood values in a data.frame for analysis. 
polydat=data.frame("HHID"=data_TL$HHID,"YEAR"=data$YEAR,"CLASS"=data$COHORT,"CURRENT"=(current),"LAGGED"=(lagged))

#remove NAs from lagged and current:
polydat=(polydat[!is.na(polydat$LAGGED),])

#check for outliers using:
plot(polydat$CURRENT~polydat$LAGGED,subset=polydat$CLASS==0)
abline(0,1)
title(main = "G1 predicted livelihood values")
plot(polydat$CURRENT~polydat$LAGGED,subset=polydat$CLASS==1)
abline(0,1)
title(main = "G2 predicted livelihood values")
plot(polydat$CURRENT~polydat$LAGGED,subset=polydat$CLASS==2)
abline(0,1)
title(main = "G3 predicted livelihood values")
plot(polydat$CURRENT~polydat$LAGGED,subset=polydat$CLASS==3)
abline(0,1)
title(main = "G4 predicted livelihood values")
#Some outliers are present at the maximum values of predicted livelihoods. We restrict the domain and range of the livelihood trajectories to to omit the bottom and top 5 percent of observations.  

#Repeat the same steps for the whole sample (to compare with sub-technology results)

#Predicted livelihood values. AIC and BIC statistics indicate that estimating liveihoods for the whole sample, a fixed effects TRANSLOG technology is most appropriate. (For the groups cobb-douglas is preferred to maintain power in the estimation)
regTL=lm(data_TL$LIVELIHOOD~as.factor(data_TL$HHID)+as.factor(data_TL$TIME)+as.matrix(data_TL[,-1:-5]))
pred0=predict(regTL)

#now get the lagged livelihoods
current0=pred0
lagged0=rep(NA,length(current0))
IDvec=data_TL$HHID
for(i in 2:length(current0)){
  if(IDvec[i]==IDvec[i-1]){
    lagged0[i]=current0[i-1]
  } else {
    lagged0[i]=NA
  }
}
alldat=data.frame("CURRENT"=current0,"LAGGED"=lagged0)


###############################################
#  LIVELIHOOD DYNAMICS - estimate livelihood dynamics
###############################################

reglist=list()

#Parametric estimation of livelihood trajectory for each group
for(mm in 1:(model_index)){
  reglist[[mm]]=lm(CURRENT~LAGGED+I(LAGGED^2)+I(LAGGED^3)+I(LAGGED^4)+I(LAGGED^5)+I(LAGGED^6),data=polydat,subset=polydat$CLASS==mm-1)
  }

#Parametric estimation of whole sample livelihood trajectory
reglist[[(mm+1)]]=lm(CURRENT~LAGGED+I(LAGGED^2)+I(LAGGED^3)+I(LAGGED^4)+I(LAGGED^5)+I(LAGGED^6),data=alldat)

#Print results
stargazer(reglist[[4]], out = "Trajectory_results.html")

#Durbin Watson tests to test for serial correlation
library (lmtest)    #for Durbin-Watson tests
Dwtests = list()
for(t in 1:(model_index+1)){
  Dwtests[[t]]= dwtest(reglist[[t]])
}

AICtests_traj = list()
for(t in 1:(model_index+1)){
  AICtests_traj[[t]]= AIC(reglist[[t]])
}
#AIC tests can be repeated for different polynomial lengths. Larger polynomials are shown to be preferred. 

#sets relevant domain and range to what we observe in the data
rangemat=matrix(NA,nrow=2,ncol=model_index)
for(mm in 1:model_index){
  rangemat[1,mm]=quantile(polydat$LAGGED[polydat$CLASS==mm-1],0.05)
  rangemat[2,mm]=quantile(polydat$LAGGED[polydat$CLASS==mm-1],0.95)
  }

#generates new data in range to predict values
newdata1=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,1],rangemat[2,1],length=1000),"LAGGED"=seq(rangemat[1,1],rangemat[2,1],length=1000))
newdata2=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,2],rangemat[2,2],length=1000),"LAGGED"=seq(rangemat[1,2],rangemat[2,2],length=1000))
newdata3=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,3],rangemat[2,3],length=1000),"LAGGED"=seq(rangemat[1,3],rangemat[2,3],length=1000))
newdata4=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,4],rangemat[2,4],length=1000),"LAGGED"=seq(rangemat[1,4],rangemat[2,4],length=1000))

newdata0=data.frame("CURRENT"=seq(min(rangemat),max(rangemat),length=3000),"LAGGED"=seq(min(rangemat),max(rangemat),length=3000))

#finds predicted value using new generated data, including 5% confidence intervals. 
pred1=(predict(reglist[[1]],newdata=newdata1, interval = 'confidence', level = 0.95))
pred2=(predict(reglist[[2]],newdata=newdata2, interval = 'confidence', level = 0.95))
pred3=(predict(reglist[[3]],newdata=newdata3, interval = 'confidence', level = 0.95))
pred4=(predict(reglist[[4]],newdata=newdata4, interval = 'confidence', level = 0.95))

pred0=(predict(reglist[[(mm+1)]],newdata=newdata0, interval = 'confidence', level = 0.95))

#Plot trajectories - with homogeneous estimation and without confidence bands. 
plot((pred1[,1])~(seq(0,rangemat[2,1],length=1000)),type="l",lwd=2,col=grey(0.2),lty=2,ylim=c(0,25),xlim=c(0,25), ylab="Predicted log livelihood",xlab="Lagged predicted log livelihood")

lines((pred2[,1])~(seq(0,rangemat[2,2],length=1000)),type="l",lwd=2,col=grey(0.4),lty=3)
lines((pred3[,1])~(seq(0,rangemat[2,3],length=1000)),type="l",lwd=2,col=grey(0.4),lty=4)
lines((pred4[,1])~(seq(0,rangemat[2,4],length=1000)),type="l",lwd=2,col=grey(0.8),lty=5)
lines((pred0[,1])~(seq(0,max(rangemat),length=3000)),type="l",lwd=2,lty=1)
axis(1, at = seq(0, 25, by = 1))
abline(a=0,b=1,lwd=1)
legend(x=8, y=5,c("G1","G2","G3","G4","ALL"),bty="n",lty=c(2,3,4,5,1),lwd=2,
       col=c(grey(0.2),grey(0.4),grey(0.6), grey(0.8),"black"))
dev.copy(png, "dynamics.png") #saves chart into working directory
dev.off()

#Plot trajectories - without homogeneous estimation and with confidence bands. 
plot((pred1[,1])~(seq(0,rangemat[2,1],length=1000)),type="l",lwd=2,col=grey(0.2),lty=2,ylim=c(0,25),xlim=c(0,25), ylab="Predicted log livelihood",xlab="Lagged predicted log livelihood")
polygon(c(rev(seq(0,rangemat[2,1],length=1000)),(seq(0,rangemat[2,1],length=1000))),c(rev(pred1[,2]),pred1[,3]),col = "grey75", border = F)
lines(pred1[,1]~seq(0,rangemat[2,1],length=1000),type="l",lwd=2,col=grey(0.2),lty=2)

polygon(c(rev(seq(0,rangemat[2,2],length=1000)),(seq(0,rangemat[2,2],length=1000))),c(rev(pred2[,2]),pred2[,3]),col = "grey75", border = F)
lines((pred2[,1])~(seq(0,rangemat[2,2],length=1000)),type="l",lwd=2,col=grey(0.4),lty=3)

polygon(c(rev(seq(0,rangemat[2,4],length=1000)),(seq(0,rangemat[2,4],length=1000))),c(rev(pred4[,2]),pred4[,3]),col = "grey95", border = F)
lines((pred4[,1])~(seq(0,rangemat[2,4],length=1000)),type="l",lwd=2,col=grey(0.8),lty=5)

polygon(c(rev(seq(0,rangemat[2,3],length=1000)),(seq(0,rangemat[2,3],length=1000))),c(rev(pred3[,2]),pred3[,3]),col = "grey75", border = F)
lines((pred3[,1])~(seq(0,rangemat[2,3],length=1000)),type="l",lwd=2,col=grey(0.4),lty=4)

abline(a=0,b=1,lwd=1)
legend(x=20, y=5,c("G1","G2","G3","G4"),bty="n",lty=c(2,3,4,5,1),lwd=2,
       col=c(grey(0.2),grey(0.4),grey(0.6), grey(0.8),"black"))
dev.copy(png, "dynamics_conf_int.png") #saves chart into working directory
dev.off()

#Plot trajectories - just homogeneous estimation and with confidence bands. 
plot((pred0[,1])~(seq(0,max(rangemat),length=3000)),type="l",lwd=2,col=grey(0.2),lty=2,ylim=c(0,25),xlim=c(0,25), ylab="Predicted log livelihood",xlab="Lagged predicted log livelihood")
polygon(c(rev(seq(0,max(rangemat),length=3000)),(seq(0,max(rangemat),length=3000))),c(rev(pred0[,2]),pred0[,3]),col = "grey75", border = F)
lines(pred0[,1]~seq(0,max(rangemat),length=3000),type="l",lwd=2,col=grey(0.2),lty=2)
abline(a=0,b=1,lwd=1)
dev.copy(png, "dynamics_conf_int_ALL.png") #saves chart into working directory
dev.off()


###############################################
#  Robustness check = heterogeneity in livelihood technology but non-parametric linear regression. 
###############################################

#This tests whether a non-parametric regression can pick up  multiple equilibria when we estimate livelihoods through the latent technology. 
library(np)

#set range for non-parametric estimation. 
range_nonpara = rep(NA,4)
range_nonpara[1]= quantile(polydat$CURRENT,0.95)
range_nonpara[2]= quantile(polydat$LAGGED,0.95)
range_nonpara[3] = quantile(polydat$LAGGED,0.05)
range_nonpara[4] = quantile(polydat$CURRENT,0.05)

polydat2 = subset(polydat,polydat$LAGGED<=range_nonpara[2]) 
polydat2 = subset(polydat,polydat$CURRENT<=range_nonpara[1])
polydat2 = subset(polydat,polydat$LAGGED>=range_nonpara[4]) 
polydat2 = subset(polydat,polydat$CURRENT>=range_nonpara[3])

reg_nonpara = npreg(CURRENT~LAGGED, data = polydat2)
plot(reg_nonpara, plot.errors.method="bootstrap", xlim = c(0,16), ylim = c(min(polydat2$LAGGED), max(polydat2$LAGGED)))
axis(2, at = seq(0, 16, by = 2))
abline(a=0,b=1,lwd=1)
#Indicates a single equilibrium point 

###############################################
#  End of main script
###############################################


###############################################
#  APPENDIX A - Model selection - Cobb-Douglas or Translog function for livelihood technologies
###############################################

#Test cobb-douglas function form and print results. 

Region = data$VILLAGE
Region[Region=="AUREPALLE"]= "Mahbubnagor"
Region[Region=="DOKUR"]= "Mahbubnagor"
Region[Region=="SHIRAPUR"]= "Solapur"
Region[Region=="KALMAN"]= "Solapur"
Region[Region=="KANZARA"]= "Akola"
Region[Region=="KINKHED"]= "Akola"


#First retrieve appropriate data from script above
data_TL = cbind(data_TL, Region) #uses Data_TL and the inverse hyperbolic income. This is used as an input for data_demeaned so should be most appropriate data
data_demeaned = data_demeaned #Demeaned data for household fixed effects

#Set up equations
eqCD_OLS = y ~ SML_LIVESTOCK_PC+LRG_LIVESTOCK_PC+BUILDINGS_PC+DRYLAND_PC+IRRIGATED_PC+MOBILE_CAPITAL_PC+SAVINGS_PC+EDUC_PC+HHSIZE #Use standard data for OLS
eqCD_fixedTIME = y ~ SML_LIVESTOCK_PC+LRG_LIVESTOCK_PC+BUILDINGS_PC+DRYLAND_PC+IRRIGATED_PC+MOBILE_CAPITAL_PC+SAVINGS_PC+EDUC_PC+HHSIZE+as.factor(TIME)#Use standard data for time 
eqCD_fixedHH = y ~ SML_LIVESTOCK_PC+LRG_LIVESTOCK_PC+BUILDINGS_PC+DRYLAND_PC+IRRIGATED_PC+MOBILE_CAPITAL_PC+SAVINGS_PC+EDUC_PC+HHSIZE #Use demeaned data for twoway fixed effects
eqCD_fixedTWOWAY = y ~ SML_LIVESTOCK_PC+LRG_LIVESTOCK_PC+BUILDINGS_PC+DRYLAND_PC+IRRIGATED_PC+MOBILE_CAPITAL_PC+SAVINGS_PC+EDUC_PC+HHSIZE+as.factor(TIME) #Use demeaned data for twoway fixed effects
eqCD_fixedTWOWAY_TIME = y ~ SML_LIVESTOCK_PC+LRG_LIVESTOCK_PC+BUILDINGS_PC+DRYLAND_PC+IRRIGATED_PC+MOBILE_CAPITAL_PC+SAVINGS_PC+EDUC_PC+ HHSIZE+as.factor(TIME) + as.factor(HHID) + (Region) + as.factor(Region)*as.factor(TIME)

Modelselectionlist = list()

#Pooled OLS model
Modelselectionlist[[1]] = lm(eqCD_OLS, data=data_TL)               #Pooled OLS

#Fixed effects models
Modelselectionlist[[2]]= lm(eqCD_fixedHH,data=data_demeaned)       #household fixed effects, demeaned data
Modelselectionlist[[3]] = lm(eqCD_fixedTIME, data = data_TL)       #time fixed effects (dummy method)
Modelselectionlist[[4]] = lm(eqCD_fixedTWOWAY, data=data_demeaned) #twoway fixed effects, demeaned data with time dummies
Modelselectionlist[[5]] = lm(eqCD_fixedTWOWAY_TIME, data=data_TL) #twoway fixed effects, demeaned data with time dummies

#Random effects models
Modelselectionlist[[6]] = plm(eqCD_OLS,data=data_TL,  model = "random", effect = "individual") #random effects - individual
Modelselectionlist[[7]] = plm(eqCD_OLS,data=data_TL,  model = "random", effect = "time")       #random effects - time
Modelselectionlist[[8]] = plm(eqCD_OLS,data=data_TL,  model = "random", effect = "twoway")     #random effects - twoway

#plm regressions for plm tests
plmregCD_hh = plm(eqCD_fixedHH,data=data_TL,  model = "within")
plmregCD_time = plm(eqCD_OLS, data = data_TL, model = "within", effect = "time")
plmregCD_twoway = plm(eqCD_OLS, data=data_TL, model = "within", effect = "twoway")

#R-squared will be different for each one as some report overall R2 (OLS, time) and the demeaned data just reports within R squared. We retrieve the overall R2 from the LSDV method. 
#The coefficients between the two are the same (as expected)

#print results 
# stargazer(Modelselectionlist[[1]], Modelselectionlist[[2]], Modelselectionlist[[3]], Modelselectionlist[[4]], Modelselectionlist[[5]], Modelselectionlist[[6]], Modelselectionlist[[7]], column.labels = c("Pooled OLS", "Fixed-HH", "Fixed-Time", "Fixed-Twoway", "Random-HH", "Random-Time", "Random-Twoway"), out = "CD model selection table.html")

#Now for tests
plmtests=phtests=Ftests=AIClist=list()

#Test random vs pooled OLS
Ftests[[1]] = pFtest(plmregCD_hh, Modelselectionlist[[1]]) 
Ftests[[2]] = pFtest(plmregCD_time, Modelselectionlist[[1]])
Ftests[[3]] = pFtest(plmregCD_twoway, Modelselectionlist[[1]]) 

#Breusch Pagan (LM) tests for fixed vs random effects
plmtests[[1]]=plmtest(plmregCD_hh, type="bp", effect = "individual")
plmtests[[2]]=plmtest(plmregCD_time, type = "bp", effect = "time")
plmtests[[3]]=plmtest(plmregCD_twoway, type="bp", effect = "twoways")

#Hausmann tests for fixed vs random effects
phtests[[1]]=phtest(plmregCD_hh, Modelselectionlist[[5]])     #hausmann test for which model to use, small p value is for a fixed effects model
phtests[[2]]=phtest(plmregCD_time, Modelselectionlist[[6]])   #hausmann test for which model to use, small p value is for a fixed effects model
phtests[[3]]=phtest(plmregCD_twoway, Modelselectionlist[[7]]) #hausmann test for which model to use, small p value is for a fixed effects model

#AIC and BIC for each model.
AICList = list() 
for (m in 1:7){
  AICvec = rep(NA,2)
  n.obs = nrow(Modelselectionlist[[m]]$model)
  resid = residuals(Modelselectionlist[[m]])
  s.sq = log((sum(resid^2)/(n.obs)))
  k = length(coef(Modelselectionlist[[m]])) + 1
  AICvec[1] = 2*k  +  n.obs * (  log(2*pi) + s.sq  + 1 ) #AIC statistic
  AICvec[2] = log(n.obs)*k +  n.obs * (  log(2*pi) + s.sq  + 1 )  #BIC statistic
  AICList[[m]] = AICvec
}

#Fixed effects models provide within R2 - so we must retrieve the overall R2. 

pred_val_mat = matrix(NA,nrow=nrow(data_TL), ncol = 2)

#For HH fixed effects model
 for(i in 1:length(IDs)){
   meanval=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
   xmat1=xmat_demeaned[data_TL$HHID==IDs[i],]
   xmat1 = xmat1[,c(-2:-14)] #remove time dummies
   coefs = Modelselectionlist[[2]]$coef
   pred=xmat1%*%coefs
   pred_vals[data_TL$HHID==IDs[i]]=meanval+pred
 }
pred_val_mat[,1] = pred_vals

#Two-way model
for(i in 1:length(IDs)){
    meanval=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
    xmat1=xmat_demeaned[data_TL$HHID==IDs[i],]
    time_dta = xmat1[,c(2:14)]
    xmat1 = xmat1[,c(-2:-14)]
    xmat1 = cbind(xmat1, time_dta) #swap time dummies to end to align with coefficients in model
    coefs = Modelselectionlist[[4]]$coef
    pred=xmat1%*%coefs
    pred_vals[data_TL$HHID==IDs[i]]=meanval+pred
}
pred_val_mat[,2] = pred_vals
  

OverallR2_vec_CD = as.data.frame(cbind(pred_val_mat[,1],pred_val_mat[,2], data_TL$LIVELIHOOD))
colnames(OverallR2_vec_CD) = c("pred_vals HH", "pred_vals Twoway", "LIVELIHOOD")

OverallR2_CD = Adj_OverallR2_CD = rep(NA, 2)

for (mm in 1:length(OverallR2_CD)){
  R2_vec = OverallR2_vec_CD[,mm]
  OverallR2_CD[mm] = cor(R2_vec, OverallR2_vec_CD$LIVELIHOOD)^2
  obs = Modelselectionlist[[2]]$fitted.values
  N = length(obs) #number of observations
  coef = Modelselectionlist[[mm*2]]$coefficients
  K = length(coef)-1  #number of coefficients for the model minus the intercept. *2 so that it selects the right model from the list 
  Adj_OverallR2_CD[mm] = 1-( ((1-OverallR2_CD[mm]) * (N-1))/(N-K-1) )
}

#Test Translog function form

#Retrieves demeaned data
dataTL_eq = data.frame(data_TL[,-1:-4])
eq_TL = as.formula(paste("y",paste(colnames(dataTL_eq), collapse = " + "), sep = " ~ "))

xmat_TL=model.matrix(eq_TL,data=data_TL)[,-1]
xmat_for_means_TL=xmat_demeaned_TL=cbind(y,xmat_TL)
means_TL=matrix(NA,nrow=length(IDs),ncol=ncol(xmat_for_means_TL))
for(rr in 1:length(IDs)){
  for(cc in 1:ncol(xmat_for_means_TL)){
    means_TL[rr,cc]=mean(xmat_for_means_TL[data_TL$HHID==IDs[rr],cc])
    xmat_demeaned_TL[data_TL$HHID==IDs[rr],cc]=xmat_for_means_TL[data_TL$HHID==IDs[rr],cc]-means_TL[rr,cc]
  }
}

data_demeaned_TL=data.frame(HHID=HHID,TIME=TIME,xmat_demeaned_TL)
data_TL = data_TL #uses Data_TL and the inverse hyperbolic income. This is used as an input for data_demeaned so should be most appropriate data
eq_TL =eq_TL

#Set up equations

eqTL_OLS = y  ~   SML_LIVESTOCK_PC + LRG_LIVESTOCK_PC + 
                  BUILDINGS_PC + DRYLAND_PC + IRRIGATED_PC + MOBILE_CAPITAL_PC + 
                  SAVINGS_PC + EDUC_PC + HHSIZE + SML_LIVESTOCK_PC2 + 
                  LRG_LIVESTOCK_PC2 + BUILDINGS_PC2 + DRYLAND_PC2 + IRRIGATED_PC2 + 
                  MOBILE_CAPITAL_PC2 + SAVINGS_PC2 + EDUC_PC2 + HHSIZE2 + 
                  SML_LIVESTOCK_PCxLRG_LIVESTOCK_PC + SML_LIVESTOCK_PCxBUILDINGS_PC + 
                  SML_LIVESTOCK_PCxDRYLAND_PC + SML_LIVESTOCK_PCxIRRIGATED_PC + 
                  SML_LIVESTOCK_PCxMOBILE_CAPITAL_PC + SML_LIVESTOCK_PCxSAVINGS_PC + 
                  SML_LIVESTOCK_PCxEDUC_PC + SML_LIVESTOCK_PCxHHSIZE + LRG_LIVESTOCK_PCxBUILDINGS_PC + 
                  LRG_LIVESTOCK_PCxDRYLAND_PC + LRG_LIVESTOCK_PCxIRRIGATED_PC + 
                  LRG_LIVESTOCK_PCxMOBILE_CAPITAL_PC + LRG_LIVESTOCK_PCxSAVINGS_PC + 
                  LRG_LIVESTOCK_PCxEDUC_PC + LRG_LIVESTOCK_PCxHHSIZE + BUILDINGS_PCxDRYLAND_PC + 
                  BUILDINGS_PCxIRRIGATED_PC + BUILDINGS_PCxMOBILE_CAPITAL_PC + 
                  BUILDINGS_PCxSAVINGS_PC + BUILDINGS_PCxEDUC_PC + BUILDINGS_PCxHHSIZE + 
                  DRYLAND_PCxIRRIGATED_PC + DRYLAND_PCxMOBILE_CAPITAL_PC + 
                  DRYLAND_PCxSAVINGS_PC + DRYLAND_PCxEDUC_PC + DRYLAND_PCxHHSIZE + 
                  IRRIGATED_PCxMOBILE_CAPITAL_PC + IRRIGATED_PCxSAVINGS_PC + 
                  IRRIGATED_PCxEDUC_PC + IRRIGATED_PCxHHSIZE + MOBILE_CAPITAL_PCxSAVINGS_PC + 
                  MOBILE_CAPITAL_PCxEDUC_PC + MOBILE_CAPITAL_PCxHHSIZE + SAVINGS_PCxEDUC_PC + 
                  SAVINGS_PCxHHSIZE + EDUC_PCxHHSIZE  #Use this for HH fixed effects but with demeaned data, and use data_TL for pooled OLS

eqTL_TIME = y ~   SML_LIVESTOCK_PC + LRG_LIVESTOCK_PC + 
                  BUILDINGS_PC + DRYLAND_PC + IRRIGATED_PC + MOBILE_CAPITAL_PC + 
                  SAVINGS_PC + EDUC_PC + HHSIZE  + SML_LIVESTOCK_PC2 + 
                  LRG_LIVESTOCK_PC2 + BUILDINGS_PC2 + DRYLAND_PC2 + IRRIGATED_PC2 + 
                  MOBILE_CAPITAL_PC2 + SAVINGS_PC2 + EDUC_PC2 + HHSIZE2 + 
                  SML_LIVESTOCK_PCxLRG_LIVESTOCK_PC + SML_LIVESTOCK_PCxBUILDINGS_PC + 
                  SML_LIVESTOCK_PCxDRYLAND_PC + SML_LIVESTOCK_PCxIRRIGATED_PC + 
                  SML_LIVESTOCK_PCxMOBILE_CAPITAL_PC + SML_LIVESTOCK_PCxSAVINGS_PC + 
                  SML_LIVESTOCK_PCxEDUC_PC + SML_LIVESTOCK_PCxHHSIZE + LRG_LIVESTOCK_PCxBUILDINGS_PC + 
                  LRG_LIVESTOCK_PCxDRYLAND_PC + LRG_LIVESTOCK_PCxIRRIGATED_PC + 
                  LRG_LIVESTOCK_PCxMOBILE_CAPITAL_PC + LRG_LIVESTOCK_PCxSAVINGS_PC + 
                  LRG_LIVESTOCK_PCxEDUC_PC + LRG_LIVESTOCK_PCxHHSIZE + BUILDINGS_PCxDRYLAND_PC + 
                  BUILDINGS_PCxIRRIGATED_PC + BUILDINGS_PCxMOBILE_CAPITAL_PC + 
                  BUILDINGS_PCxSAVINGS_PC + BUILDINGS_PCxEDUC_PC + BUILDINGS_PCxHHSIZE + 
                  DRYLAND_PCxIRRIGATED_PC + DRYLAND_PCxMOBILE_CAPITAL_PC + 
                  DRYLAND_PCxSAVINGS_PC + DRYLAND_PCxEDUC_PC + DRYLAND_PCxHHSIZE + 
                  IRRIGATED_PCxMOBILE_CAPITAL_PC + IRRIGATED_PCxSAVINGS_PC + 
                  IRRIGATED_PCxEDUC_PC + IRRIGATED_PCxHHSIZE + MOBILE_CAPITAL_PCxSAVINGS_PC + 
                  MOBILE_CAPITAL_PCxEDUC_PC + MOBILE_CAPITAL_PCxHHSIZE + SAVINGS_PCxEDUC_PC + 
                  SAVINGS_PCxHHSIZE + EDUC_PCxHHSIZE +
                  as.factor(TIME)


Modelselectionlist_TL = list()

#Pooled OLS model
Modelselectionlist_TL[[1]] = lm(eqTL_OLS, data=data_TL) #Pooled OLS

#Fixed effects models
Modelselectionlist_TL[[2]] = lm(eqTL_OLS,data=data_demeaned_TL)   #household fixed effects, demeaned data
Modelselectionlist_TL[[3]]= lm(eqTL_TIME, data = data_TL)         #time fixed effects (dummy method)
Modelselectionlist_TL[[4]] = lm(eqTL_TIME, data=data_demeaned_TL) #twoway fixed effects, demeaned data with time dummies

#Random effects models
Modelselectionlist_TL[[5]] = plm(eqTL_OLS,data=data_TL,  model = "random", effect = "time", random.method = "walhus")    #random effects - individual 
Modelselectionlist_TL[[6]] = plm(eqTL_OLS,data=data_TL,  model = "random", effect = "individual")                        #random effects - time
Modelselectionlist_TL[[7]] = plm(eqTL_OLS,data=data_TL,  model = "random", effect = "twoway", random.method = "walhus")  #random effects - twoway

#plm regressions for plm tests
plmregTL_time = plm(eqTL_OLS, data = data_TL, model = "within", effect = "time")
plmregTL_hh = plm(eqTL_OLS,data=data_TL,  model = "within", effect = "individual")
plmregTL_twoway = plm(eqTL_OLS, data=data_TL, model = "within", effect = "twoway")

#print results in CSV file
# stargazer(Modelselectionlist_TL[[1]], Modelselectionlist_TL[[2]], Modelselectionlist_TL[[3]], Modelselectionlist_TL[[4]], Modelselectionlist_TL[[5]], Modelselectionlist_TL[[6]], Modelselectionlist_TL[[7]],column.labels = c("Pooled OLS", "Fixed-HH", "Fixed-Time", "Fixed-Twoway", "Random-HH", "Random-Time", "Random-Twoway"), out = "TL model selection table.html")

#Now for tests
plmtests_TL=phtests_TL=Ftests_TL=list()

#Test random vs pooled OLS
Ftests_TL[[1]] = pFtest(plmregTL_hh, Modelselectionlist_TL[[1]]) 
Ftests_TL[[2]] = pFtest(plmregTL_time, Modelselectionlist_TL[[1]])
Ftests_TL[[3]] = pFtest(plmregTL_twoway, Modelselectionlist_TL[[1]]) 

#Breusch Pagan (LM) tests for fixed vs random effects
plmtests_TL[[1]]=plmtest(plmregTL_hh, type="bp", effect = "individual")
plmtests_TL[[2]]=plmtest(plmregTL_time, type = "bp", effect = "time")
plmtests_TL[[3]]=plmtest(plmregTL_twoway, type="bp", effect = "twoways")

#Hausmann tests for fixed vs random effects
phtests_TL[[1]]=phtest(plmregTL_hh, Modelselectionlist_TL[[5]])     #hausmann test for which model to use, small p value is for a fixed effects model
phtests_TL[[2]]=phtest(plmregTL_time, Modelselectionlist_TL[[6]])   #hausmann test for which model to use, small p value is for a fixed effects model
phtests_TL[[3]]=phtest(plmregTL_twoway, Modelselectionlist_TL[[7]]) #hausmann test for which model to use, small p value is for a fixed effects model

#AIC and BIC for each model.
AICList_TL = list() 
for (m in 1:7){
  AICvec = rep(NA,2)
  n.obs = nrow(Modelselectionlist_TL[[m]]$model)
  resid = residuals(Modelselectionlist_TL[[m]])
  s.sq = log((sum(resid^2)/(n.obs)))
  k = length(coef(Modelselectionlist_TL[[m]])) + 1
  AICvec[1] = 2*k  +  n.obs * (  log(2*pi) + s.sq  + 1 ) #AIC statistic
  AICvec[2] = log(n.obs)*k +  n.obs * (  log(2*pi) + s.sq  + 1 )  #BIC statistic
  AICList_TL[[m]] = AICvec
}


#Fixed effects models provide within R2 - so we must retrieve the overall R2. 

pred_val_mat_TL = matrix(NA,nrow=nrow(data_TL), ncol = 2)

#For HH fixed effects model
TL_model_HH = model.matrix(eqTL_OLS, data = data_TL)
for(i in 1:length(IDs)){
  meanval=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
  xmat1= TL_model_HH[data_TL$HHID==IDs[i],]
  coefs = Modelselectionlist_TL[[2]]$coefficients
  pred=xmat1%*%coefs
  pred_vals[data_TL$HHID==IDs[i]]=meanval+pred
}
pred_val_mat_TL[,1] = pred_vals

#Two-way model
TL_model_TW = model.matrix(eqTL_TIME, data = data_TL)
for(i in 1:length(IDs)){
  meanval=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
  xmat1=TL_model_TW[data_TL$HHID==IDs[i],]
  coefs = Modelselectionlist_TL[[4]]$coefficients
  pred=xmat1%*%coefs
  pred_vals[data_TL$HHID==IDs[i]]=meanval+pred
}
pred_val_mat_TL[,2] = pred_vals


OverallR2_vec_TL = as.data.frame(cbind(pred_val_mat_TL[,1],pred_val_mat_TL[,2], data_TL$LIVELIHOOD))
colnames(OverallR2_vec_TL) = c("pred_vals HH", "pred_vals Twoway", "LIVELIHOOD")

OverallR2_TL = Adj_OverallR2_TL = rep(NA, 2)

for (mm in 1:length(OverallR2_TL)){
  R2_vec = OverallR2_vec_TL[,mm]
  OverallR2_TL[mm] = cor(R2_vec, OverallR2_vec_TL$LIVELIHOOD)^2
  obs = Modelselectionlist[[2]]$fitted.values
  N = length(obs) #number of observations
  coef = Modelselectionlist[[mm*2]]$coefficients
  K = length(coef)-1  #number of coefficients for the model minus the intercept. *2 so that it selects the right model from the list 
  Adj_OverallR2_TL[mm] = 1-( ((1-OverallR2_TL[mm]) * (N-1))/(N-K-1) )
}



###############################################
#  APPENDIX B - Ex-ante grouping by village
###############################################
#This section follows reviewer comments to seperate the sample based on the village cohorts. 
#We run the same script as for the latent groupings to get both technology estimates and trajectories. 

n_classes = 6 #number of villages

classvec=seq(1,n_classes)
polyvec=seq(1,n_polynomials)

classmat=matrix(NA,nrow=nrow(data),ncol=3)
classmat[,1]=data$HHID
classmat[,2] = data$VILLAGE

#Allocate group number for each village
for(i in 1:nrow(classmat)){
    if(classmat[i,2]=="AUREPALLE"){
      classmat[i,3]=1
    } else if (classmat[i,2]=="DOKUR") {
      classmat[i,3]=2
    } else if (classmat[i,2]=="KALMAN") {
      classmat[i,3]=3
    } else if (classmat[i,2]=="KANZARA") {
      classmat[i,3]=4
    } else if (classmat[i,2]=="KINKHED") {
      classmat[i,3]=5
    } else if (classmat[i,2]=="SHIRAPUR"){
      classmat[i,3]=6 
    }
}

classmat1=sapply(data.frame(classmat[,3]), function(x) as.numeric(as.character(x)))#turns group cohorts from character strings to recognisable elements
classmat=data.frame(classmat[,1],classmat1) #returns HH ids
data$COHORT = classes = classmat[,2]

#SUMMARY STATISTICS OF ASSET LEVELS (by village)

#Subset data for each group
datasummary_aurepelle = subset(datasummary_group, data$COHORT==1)
datasummary_dokur = subset(datasummary_group, data$COHORT==2)
datasummary_kalman = subset(datasummary_group, data$COHORT==3)
datasummary_kanzara = subset(datasummary_group, data$COHORT==4)
datasummary_kinkhed = subset(datasummary_group, data$COHORT==5)
datasummary_shirapur = subset(datasummary_group, data$COHORT==6)

#Summary statistics table (mean and standard error for each group)
assetsummary_groups=matrix(NA,ncol=12,nrow=ncol(datasummary_group[,-1:-5]))
for(r in 1:nrow(assetsummary_groups)){
  assetsummary_groups[r,1] = mean((datasummary_aurepelle[,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,2] = sqrt(var((datasummary_aurepelle [,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_aurepelle))
  assetsummary_groups[r,3] = mean((datasummary_dokur[,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,4] = sqrt(var((datasummary_dokur[,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_dokur))
  assetsummary_groups[r,5] = mean((datasummary_kalman [,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,6] = sqrt(var((datasummary_kalman [,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_kalman))
  assetsummary_groups[r,7] = mean((datasummary_kanzara[,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,8] = sqrt(var((datasummary_kanzara[,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_kanzara))
  assetsummary_groups[r,9] = mean((datasummary_kinkhed [,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,10] = sqrt(var((datasummary_kinkhed [,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_kinkhed))
  assetsummary_groups[r,11] = mean((datasummary_shirapur[,-1:-5][,r]),na.rm=TRUE)
  assetsummary_groups[r,12] = sqrt(var((datasummary_shirapur[,-1:-5][,r]),na.rm=TRUE))/sqrt(nrow(datasummary_shirapur))
}

rownames(assetsummary_groups) = colnames(datasummary_group[,-1:-5])
colnames(assetsummary_groups) = c("mean","std.err","mean","std.err","mean","std.err","mean","std.err","mean","std.err", "mean","std.err")
assetsummary_groups = assetsummary_groups[!is.na(assetsummary_groups[1,]),] #removes NA rows where variable is a character string and not numeric (e.g. Caste, occupation etc)
write.table(assetsummary_groups, file = "summarystats_byGroup.csv", sep=",")


#  ESTIMATE MODELS for villages
model1 = lm(yd~0+xmat_demeaned  ,subset=data$COHORT == 1)
model2= lm(yd~0+xmat_demeaned   ,subset=data$COHORT == 2)
model3 = lm(yd~0+xmat_demeaned  ,subset=data$COHORT == 3)
model4 = lm(yd~0+xmat_demeaned  ,subset=data$COHORT == 4)
model5 = lm(yd~0+xmat_demeaned  ,subset=data$COHORT == 5)
model6 = lm(yd~0+xmat_demeaned  ,subset=data$COHORT == 6)

modelsList = list(model1, model2, model3, model4, model5, model6)

#print results in excel table, including relevant statistics
len = length(modelsList[[1]]$coefficients)
namesX = names(modelsList[[1]]$coefficients)
namesX = substr(namesX,14,100)
modelcoefficients = matrix(NA, nrow = len+3, ncol=5)
colnames(modelcoefficients) = c("", "Aurepalle", "Dokur", "Kalman", "Kanzara", "Kinkhed", "Shirapur")
for (ccc in 1:model_index){
  coefs = summary(modelsList[[ccc]])$coef
  for(rrr in 1:nrow(coefs)){
    est = round(coefs[rrr,1],4)
    p_stars = ifelse(coefs[rrr,4]>0.1,"",ifelse(coefs[rrr,4]>0.05,"*",ifelse(coefs[rrr,4]>0.01,"**","***")))
    modelcoefficients[rrr,ccc+1] = paste(est,p_stars,sep="")
  }
  modelcoefficients[c(len+1),ccc+1] = sum((classvec+1)==ccc)
  modelcoefficients[c(len+2),ccc+1] = summary(modelsList[[ccc]])$r.squared
  modelcoefficients[c(len+3),ccc+1] = summary(modelsList[[ccc]])$adj.r.squared
}
modelcoefficients[,1] = c(namesX, "HH in group", "R2", "adjusted R2" )
write.table(modelcoefficients, "model results table.csv", sep=",")

#demeaned lm models report the within R2. Below we estimate the overall R2 by adding the levels back in. The predicted values are also needed for the trajecotry functions. 

#For each person get their actual predicted value with their mean value - this replaces the pred_vals calculated without the inverse mills ratio. 

xmat_demeaned = model.matrix(yd~0+xmat_demeaned) 
pred_vals=rep(NA,nrow=nrow(data_TL))
for(i in 1:length(IDs)){
  meanval=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
  group = as.numeric(levels(as.factor(data$COHORT[data_TL$HHID==IDs[i]])))
  group=group
  xmat1=xmat_demeaned[data_TL$HHID==IDs[i],]
  pred=xmat1%*%modelsList[[group]]$coef
  pred_vals[data_TL$HHID==IDs[i]]=meanval+pred
}

OverallR2_vec = as.data.frame(cbind(pred_vals, data_TL$LIVELIHOOD, data$COHORT))
colnames(OverallR2_vec) = c("pred_vals", "LIVELIHOOD", "village")

#repeat this for the full sample
pred_vals0=rep(NA,nrow=nrow(data_TL))
for(i in 1:length(IDs)){
  meanval0=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
  xmat10=xmat_demeaned[data_TL$HHID==IDs[i],]
  pred0=xmat10%*%model0$coef
  pred_vals0[data_TL$HHID==IDs[i]]=meanval0+pred0
}

OverallR2 = Adj_OverallR2 = rep(NA, 5)


for (g in 1:model_index){
  R2_vec = subset(OverallR2_vec,OverallR2_vec$COHORT==g-1)
  OverallR2[g+1] = cor(R2_vec$pred_vals, R2_vec$LIVELIHOOD)^2
  N = length(modelsList[[g]]$fitted.values) #number of observations
  K = ncol(xmat_demeaned)+1  #number of coefficients plus inverse mills ratio
  Adj_OverallR2[g+1] = 1-( ((1-OverallR2[g+1]) * (N-1))/(N-K-1) )
}

#test pairwise differences in coefficients between the groups using z-scores

modelsList2 = modelsList 
modelsList2[[7]] = model0

Coeftests = list() #stores pairwise comparisons in a list, each element in the list is a different coefficient
for (l in 1:length(modelsList2[[7]]$coefficients)){ #ignore inverse mills ratio as not relevant for village analysis
  coeftestmat = matrix(NA, ncol = 7, nrow = 7)
  rownames(coeftestmat) = c(colnames(xmat_demeaned)[l],1,2,3,4,5,6)
  colnames(coeftestmat) = c("WS",1,2,3,4,5,6)
  for (r in 1:7){
    beta1 = modelsList2[[r]]$coefficients[[l]] 
    se1 = coef(summary(modelsList2[[r]]))[l, "Std. Error"]
    for (c in 1:7){
      beta2 = modelsList2[[c]]$coefficients[[l]] 
      se2 = coef(summary(modelsList2[[c]]))[l, "Std. Error"]
      coeftestmat[r,c] = z.score(beta1, beta2, se1, se2)
    }
  }
  Coeftests[[l]] = coeftestmat
}


#  LIVELIHOOD DYNAMICS - retrieve predicted livelihoods
#First retrieve the predicted livelihood values and lagged livelihood values from the estimated livelihood technologies. 

pred_vals = pred_vals #already calculated for the overall R2. This reintroduces mean levels from demeaned data

#Retrieve lagged predicted livelihoods from predicted values. 
current=pred_vals
lagged=rep(NA,length(current))
IDvec=data_TL$HHID
for(i in 2:length(current)){
  if(IDvec[i]==IDvec[i-1]){
    lagged[i]=current[i-1]
  } else {
    lagged[i]=NA
  }
}

#Save predicted livelihood values in a data.frame for analysis. 
polydat=data.frame("HHID"=data_TL$HHID,"YEAR"=data$YEAR,"CLASS"=data$COHORT,"CURRENT"=(current),"LAGGED"=(lagged))

#remove NAs from lagged and current:
polydat=(polydat[!is.na(polydat$LAGGED),])

#  LIVELIHOOD DYNAMICS - estimate livelihood dynamics

reglist=list()

#Parametric estimation of livelihood trajectory for each group
for(mm in 1:6){
  print(mm)
  reglist[[mm]]=lm(CURRENT~LAGGED+I(LAGGED^2)+I(LAGGED^3)+I(LAGGED^4)+I(LAGGED^5)+I(LAGGED^6),data=polydat,subset=polydat$CLASS==mm)
}

#Parametric estimation of whole sample livelihood trajectory
reglist[[(mm+1)]]=lm(CURRENT~LAGGED+I(LAGGED^2)+I(LAGGED^3)+I(LAGGED^4)+I(LAGGED^5)+I(LAGGED^6),data=alldat)

#Print results
# stargazer(reglist[[5]], reglist[[1]], reglist[[2]], reglist[[3]], reglist[[4]], reglist[[5]], reglist[[6]], out = "Trajectory results.html", column.labels = c("Whole sample", "G1", "G2", "G3", "G4". "G5", "G6"))

#sets relevant domain and range to what we observe in the data
rangemat=matrix(NA,nrow=2,ncol=6)
for(mm in 1:6){
  rangemat[1,mm]=quantile(polydat$LAGGED[polydat$CLASS==mm],0.05)
  rangemat[2,mm]=quantile(polydat$LAGGED[polydat$CLASS==mm],0.95)
}

#generates new data in range to predict values
newdata1=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,1],rangemat[2,1],length=1000),"LAGGED"=seq(rangemat[1,1],rangemat[2,1],length=1000))
newdata2=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,2],rangemat[2,2],length=1000),"LAGGED"=seq(rangemat[1,2],rangemat[2,2],length=1000))
newdata3=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,3],rangemat[2,3],length=1000),"LAGGED"=seq(rangemat[1,3],rangemat[2,3],length=1000))
newdata4=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,4],rangemat[2,4],length=1000),"LAGGED"=seq(rangemat[1,4],rangemat[2,4],length=1000))
newdata5=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,4],rangemat[2,5],length=1000),"LAGGED"=seq(rangemat[1,4],rangemat[2,5],length=1000))
newdata6=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,4],rangemat[2,6],length=1000),"LAGGED"=seq(rangemat[1,4],rangemat[2,6],length=1000))

#finds predicted value using new generated data, including 5% confidence intervals. 
pred1=(predict(reglist[[1]],newdata=newdata1, interval = 'confidence', level = 0.95))
pred2=(predict(reglist[[2]],newdata=newdata2, interval = 'confidence', level = 0.95))
pred3=(predict(reglist[[3]],newdata=newdata3, interval = 'confidence', level = 0.95))
pred4=(predict(reglist[[4]],newdata=newdata4, interval = 'confidence', level = 0.95))
pred5=(predict(reglist[[5]],newdata=newdata4, interval = 'confidence', level = 0.95))
pred6=(predict(reglist[[6]],newdata=newdata4, interval = 'confidence', level = 0.95))


#Plot trajectories - with homogeneous estimation and without confidence bands. 
plot((pred1[,1])~(seq(0,rangemat[2,1],length=1000)),type="l",lwd=2,col=grey(0.2),lty=2,ylim=c(0,25),xlim=c(0,25), ylab="Predicted log livelihood",xlab="Lagged predicted log livelihood")
lines((pred2[,1])~(seq(0,rangemat[2,2],length=1000)),type="l",lwd=2,col=grey(0.4),lty=3)
lines((pred3[,1])~(seq(0,rangemat[2,3],length=1000)),type="l",lwd=2,col=grey(0.4),lty=4)
lines((pred4[,1])~(seq(0,rangemat[2,4],length=1000)),type="l",lwd=2,col=grey(0.8),lty=5)
lines((pred5[,1])~(seq(0,rangemat[2,5],length=1000)),type="l",lwd=2,col=grey(0.8),lty=6)
lines((pred6[,1])~(seq(0,rangemat[2,6],length=1000)),type="l",lwd=2,col=grey(0.8),lty=7)

abline(a=0,b=1,lwd=1)
legend(x=8, y=5,c("AU","DOK", "KAL","KAN","KIN","SHI"),bty="n",lty=c(2,3,4,5,6,7),lwd=2,
       col=c(grey(0.2),grey(0.4),grey(0.6), grey(0.8), grey(0.8), grey(0.8)))
dev.copy(png, "dynamicsVIL.png") #saves chart into working directory
dev.off()

#Plot trajectories - without homogeneous estimation and with confidence bands. 
plot((pred1[,1])~(seq(0,rangemat[2,1],length=1000)),type="l",lwd=2,col=grey(0.2),lty=2,ylim=c(0,25),xlim=c(0,25), ylab="Predicted log livelihood",xlab="Lagged predicted log livelihood")
polygon(c(rev(seq(0,rangemat[2,1],length=1000)),(seq(0,rangemat[2,1],length=1000))),c(rev(pred1[,2]),pred1[,3]),col = "grey75", border = F)
lines(pred1[,1]~seq(0,rangemat[2,1],length=1000),type="l",lwd=2,col=grey(0.2),lty=2)

polygon(c(rev(seq(0,rangemat[2,2],length=1000)),(seq(0,rangemat[2,2],length=1000))),c(rev(pred2[,2]),pred2[,3]),col = "grey75", border = F)
lines((pred2[,1])~(seq(0,rangemat[2,2],length=1000)),type="l",lwd=2,col=grey(0.4),lty=3)

polygon(c(rev(seq(0,rangemat[2,4],length=1000)),(seq(0,rangemat[2,4],length=1000))),c(rev(pred4[,2]),pred4[,3]),col = "grey95", border = F)
lines((pred4[,1])~(seq(0,rangemat[2,4],length=1000)),type="l",lwd=2,col=grey(0.8),lty=5)

polygon(c(rev(seq(0,rangemat[2,3],length=1000)),(seq(0,rangemat[2,3],length=1000))),c(rev(pred3[,2]),pred3[,3]),col = "grey75", border = F)
lines((pred3[,1])~(seq(0,rangemat[2,3],length=1000)),type="l",lwd=2,col=grey(0.4),lty=4)

polygon(c(rev(seq(0,rangemat[2,5],length=1000)),(seq(0,rangemat[2,5],length=1000))),c(rev(pred5[,2]),pred5[,3]),col = "grey85", border = F)
lines((pred5[,1])~(seq(0,rangemat[2,5],length=1000)),type="l",lwd=2,col=grey(0.8),lty=6)

polygon(c(rev(seq(0,rangemat[2,6],length=1000)),(seq(0,rangemat[2,6],length=1000))),c(rev(pred6[,2]),pred6[,3]),col = "grey75", border = F)
lines((pred6[,1])~(seq(0,rangemat[2,6],length=1000)),type="l",lwd=2,col=grey(0.8),lty=7)

abline(a=0,b=1,lwd=1)
legend(x=8, y=5,c("AU","DOK", "KAL","KAN","KIN","SHI"),bty="n",lty=c(2,3,4,5,6,7),lwd=2,
       col=c(grey(0.2),grey(0.4),grey(0.6), grey(0.8), grey(0.8), grey(0.8)))
dev.copy(png, "dynamicsVIL_conf_int.png") #saves chart into working directory
dev.off()

#At a village level there are not many differences between the equilibria. This supports the use of latent estimations. 


###############################################
#  APPENDIX C - Ex-ante grouping by landholdings
###############################################
#This section follows reviewer comments to seperate the sample based on the village cohorts. 
#We run the same script as for the latent groupings to get both technology estimates and trajectories. 

n_classes = 4 

classvec=seq(1,n_classes)
polyvec=seq(1,n_polynomials)

#Landholding quantiles


#Allocate group number for each land holding in 2001

data2001 = subset(data, data$YEAR==2001)
land = data2001$IRRIGATED + data2001$DRYLAND

classmat=matrix(NA,nrow=nrow(data2001),ncol=3)
classmat[,1]= data2001$HHID
classmat[,2] = land

for(i in 1:length(land)){
  if(classmat[i,2]<=quantile(land,0.25)){
    classmat[i,3]=1
  } else if (classmat[i,2]<=quantile(land,0.5) & classmat[i,2] > quantile(land,0.25)) {
    classmat[i,3]=2
  } else if (classmat[i,2]<=quantile(land,0.75) & classmat[i,2] > quantile(land,0.5)) {
    classmat[i,3]=3
  } else {
    classmat[i,3]=4
  }
}

classmat1=sapply(data.frame(classmat[,3]), function(x) as.numeric(as.character(x)))#turns group cohorts from character strings to recognisable elements
classmat=data.frame(classmat[,1],classmat1) #returns HH ids

#allocated classes to rest of data

for (i in 1:nrow(classmat)){
  data[data$HHID==classmat[i,1],]$COHORT = classmat[i,2]
}
classes = data$COHORT

#  ESTIMATE MODELS for villages
model1 = lm(yd~0+xmat_demeaned  ,subset=data$COHORT == 1)
model2= lm(yd~0+xmat_demeaned   ,subset=data$COHORT == 2)
model3 = lm(yd~0+xmat_demeaned  ,subset=data$COHORT == 3)
model4 = lm(yd~0+xmat_demeaned  ,subset=data$COHORT == 4)

modelsList = list(model1, model2, model3, model4)

#print results in excel table, including relevant statistics
len = length(modelsList[[1]]$coefficients)
namesX = names(modelsList[[1]]$coefficients)
namesX = substr(namesX,14,100)
modelcoefficients = matrix(NA, nrow = len+3, ncol=5)
colnames(modelcoefficients) = c("", "25%", "50%", "75%", "100%")
for (ccc in 1:model_index){
  coefs = summary(modelsList[[ccc]])$coef
  for(rrr in 1:nrow(coefs)){
    est = round(coefs[rrr,1],4)
    p_stars = ifelse(coefs[rrr,4]>0.1,"",ifelse(coefs[rrr,4]>0.05,"*",ifelse(coefs[rrr,4]>0.01,"**","***")))
    modelcoefficients[rrr,ccc+1] = paste(est,p_stars,sep="")
  }
  modelcoefficients[c(len+1),ccc+1] = sum((classvec+1)==ccc)
  modelcoefficients[c(len+2),ccc+1] = summary(modelsList[[ccc]])$r.squared
  modelcoefficients[c(len+3),ccc+1] = summary(modelsList[[ccc]])$adj.r.squared
}
modelcoefficients[,1] = c(namesX, "HH in group", "R2", "adjusted R2" )
write.table(modelcoefficients, "model results table.csv", sep=",")

#demeaned lm models report the within R2. Below we estimate the overall R2 by adding the levels back in. The predicted values are also needed for the trajecotry functions. 

#For each person get their actual predicted value with their mean value - this replaces the pred_vals calculated without the inverse mills ratio. 

xmat_demeaned = model.matrix(yd~0+xmat_demeaned) 
pred_vals=rep(NA,nrow=nrow(data_TL))
for(i in 1:length(IDs)){
  meanval=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
  group = as.numeric(levels(as.factor(data$COHORT[data_TL$HHID==IDs[i]])))
  xmat1 =xmat_demeaned[data_TL$HHID==IDs[i],]
  pred = xmat1%*%modelsList[[group]]$coef
  pred_vals[data_TL$HHID == IDs[i]]=meanval+pred
}

OverallR2_vec = as.data.frame(cbind(pred_vals, data_TL$LIVELIHOOD, data$COHORT))
colnames(OverallR2_vec) = c("pred_vals", "LIVELIHOOD", "land cohort")

#repeat this for the full sample
pred_vals0=rep(NA,nrow=nrow(data_TL))
for(i in 1:length(IDs)){
  meanval0=mean(data_TL$LIVELIHOOD[data_TL$HHID==IDs[i]])
  xmat10=xmat_demeaned[data_TL$HHID==IDs[i],]
  pred0=xmat10%*%model0$coef
  pred_vals0[data_TL$HHID==IDs[i]]=meanval0+pred0
}

OverallR2 = Adj_OverallR2 = rep(NA, 5)


for (g in 1:model_index){
  R2_vec = subset(OverallR2_vec,OverallR2_vec$COHORT==g-1)
  OverallR2[g+1] = cor(R2_vec$pred_vals, R2_vec$LIVELIHOOD)^2
  N = length(modelsList[[g]]$fitted.values) #number of observations
  K = ncol(xmat_demeaned)+1  #number of coefficients plus inverse mills ratio
  Adj_OverallR2[g+1] = 1-( ((1-OverallR2[g+1]) * (N-1))/(N-K-1) )
}

#test pairwise differences in coefficients between the groups using z-scores

modelsList2 = modelsList 
modelsList2[[7]] = model0

Coeftests = list() #stores pairwise comparisons in a list, each element in the list is a different coefficient
for (l in 1:length(modelsList2[[7]]$coefficients)){ #ignore inverse mills ratio as not relevant for village analysis
  coeftestmat = matrix(NA, ncol = 7, nrow = 7)
  rownames(coeftestmat) = c(colnames(xmat_demeaned)[l],1,2,3,4,5,6)
  colnames(coeftestmat) = c("WS",1,2,3,4,5,6)
  for (r in 1:7){
    beta1 = modelsList2[[r]]$coefficients[[l]] 
    se1 = coef(summary(modelsList2[[r]]))[l, "Std. Error"]
    for (c in 1:7){
      beta2 = modelsList2[[c]]$coefficients[[l]] 
      se2 = coef(summary(modelsList2[[c]]))[l, "Std. Error"]
      coeftestmat[r,c] = z.score(beta1, beta2, se1, se2)
    }
  }
  Coeftests[[l]] = coeftestmat
}


#  LIVELIHOOD DYNAMICS - retrieve predicted livelihoods
#First retrieve the predicted livelihood values and lagged livelihood values from the estimated livelihood technologies. 

pred_vals = pred_vals #already calculated for the overall R2. This reintroduces mean levels from demeaned data

#Retrieve lagged predicted livelihoods from predicted values. 
current=pred_vals
lagged=rep(NA,length(current))
IDvec=data_TL$HHID
for(i in 2:length(current)){
  if(IDvec[i]==IDvec[i-1]){
    lagged[i]=current[i-1]
  } else {
    lagged[i]=NA
  }
}

#Save predicted livelihood values in a data.frame for analysis. 
polydat=data.frame("HHID"=data_TL$HHID,"YEAR"=data$YEAR,"CLASS"=data$COHORT,"CURRENT"=(current),"LAGGED"=(lagged))

#remove NAs from lagged and current:
polydat=(polydat[!is.na(polydat$LAGGED),])

#  LIVELIHOOD DYNAMICS - estimate livelihood dynamics

reglist=list()

#Parametric estimation of livelihood trajectory for each group
for(mm in 1:4){
  reglist[[mm]]=lm(CURRENT~LAGGED+I(LAGGED^2)+I(LAGGED^3)+I(LAGGED^4)+I(LAGGED^5)+I(LAGGED^6),data=polydat,subset=polydat$CLASS==mm)
}

#Parametric estimation of whole sample livelihood trajectory
reglist[[(mm+1)]]=lm(CURRENT~LAGGED+I(LAGGED^2)+I(LAGGED^3)+I(LAGGED^4)+I(LAGGED^5)+I(LAGGED^6),data=alldat)

#sets relevant domain and range to what we observe in the data
rangemat=matrix(NA,nrow=2,ncol=4)
for(mm in 1:4){
  rangemat[1,mm]=quantile(polydat$LAGGED[polydat$CLASS==mm],0.05)
  rangemat[2,mm]=quantile(polydat$LAGGED[polydat$CLASS==mm],0.95)
}

#generates new data in range to predict values
newdata1=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,1],rangemat[2,1],length=1000),"LAGGED"=seq(rangemat[1,1],rangemat[2,1],length=1000))
newdata2=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,2],rangemat[2,2],length=1000),"LAGGED"=seq(rangemat[1,2],rangemat[2,2],length=1000))
newdata3=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,3],rangemat[2,3],length=1000),"LAGGED"=seq(rangemat[1,3],rangemat[2,3],length=1000))
newdata4=data.frame("HHID"=NA,"CLASS"=rep(1,1000),"CURRENT"=seq(rangemat[1,4],rangemat[2,4],length=1000),"LAGGED"=seq(rangemat[1,4],rangemat[2,4],length=1000))

#finds predicted value using new generated data, including 5% confidence intervals. 
pred1=(predict(reglist[[1]],newdata=newdata1, interval = 'confidence', level = 0.95))
pred2=(predict(reglist[[2]],newdata=newdata2, interval = 'confidence', level = 0.95))
pred3=(predict(reglist[[3]],newdata=newdata3, interval = 'confidence', level = 0.95))
pred4=(predict(reglist[[4]],newdata=newdata4, interval = 'confidence', level = 0.95))


#Plot trajectories - with homogeneous estimation and without confidence bands. 
plot((pred1[,1])~(seq(0,rangemat[2,1],length=1000)),type="l",lwd=2,col=grey(0.2),lty=2,ylim=c(0,25),xlim=c(0,25), ylab="Predicted log livelihood",xlab="Lagged predicted log livelihood")
lines((pred2[,1])~(seq(0,rangemat[2,2],length=1000)),type="l",lwd=2,col=grey(0.4),lty=3)
lines((pred3[,1])~(seq(0,rangemat[2,3],length=1000)),type="l",lwd=2,col=grey(0.4),lty=4)
lines((pred4[,1])~(seq(0,rangemat[2,4],length=1000)),type="l",lwd=2,col=grey(0.8),lty=5)


abline(a=0,b=1,lwd=1)
legend(x=8, y=5,c("25%","50%","75%","100%"),bty="n",lty=c(2,3,4,5),lwd=2,
       col=c(grey(0.2),grey(0.4),grey(0.6), grey(0.8)))
dev.copy(png, "dynamicsLand.png") #saves chart into working directory
dev.off()

#Plot trajectories - without homogeneous estimation and with confidence bands. 
plot((pred1[,1])~(seq(0,rangemat[2,1],length=1000)),type="l",lwd=2,col=grey(0.2),lty=2,ylim=c(0,25),xlim=c(0,25), ylab="Predicted log livelihood",xlab="Lagged predicted log livelihood")
polygon(c(rev(seq(0,rangemat[2,1],length=1000)),(seq(0,rangemat[2,1],length=1000))),c(rev(pred1[,2]),pred1[,3]),col = "grey75", border = F)
lines(pred1[,1]~seq(0,rangemat[2,1],length=1000),type="l",lwd=2,col=grey(0.2),lty=2)

polygon(c(rev(seq(0,rangemat[2,2],length=1000)),(seq(0,rangemat[2,2],length=1000))),c(rev(pred2[,2]),pred2[,3]),col = "grey75", border = F)
lines((pred2[,1])~(seq(0,rangemat[2,2],length=1000)),type="l",lwd=2,col=grey(0.4),lty=3)

polygon(c(rev(seq(0,rangemat[2,4],length=1000)),(seq(0,rangemat[2,4],length=1000))),c(rev(pred4[,2]),pred4[,3]),col = "grey95", border = F)
lines((pred4[,1])~(seq(0,rangemat[2,4],length=1000)),type="l",lwd=2,col=grey(0.8),lty=5)

polygon(c(rev(seq(0,rangemat[2,3],length=1000)),(seq(0,rangemat[2,3],length=1000))),c(rev(pred3[,2]),pred3[,3]),col = "grey75", border = F)
lines((pred3[,1])~(seq(0,rangemat[2,3],length=1000)),type="l",lwd=2,col=grey(0.4),lty=4)


abline(a=0,b=1,lwd=1)
legend(x=8, y=5,c("25%","50%","75%","100%"),bty="n",lty=c(2,3,4,5),lwd=2,
       col=c(grey(0.2),grey(0.4),grey(0.6), grey(0.8)))
dev.copy(png, "dynamicsLAND_conf_int.png") #saves chart into working directory
dev.off()

#At a landholdings there are not many differences between the equilibria. This supports the use of latent estimations. 

###############################################
#  END OF SCRIPT
###############################################
