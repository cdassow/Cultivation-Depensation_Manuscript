# CJD 2.28.22
# REMEMBER TO COMMENT OUT OR DELETE FILE PATHS BEFORE SUBMISSION

##### CHUNK OF DATA PREP BEFORE COPYING VIGNETTE #####
setwd("~/OAS_git_repos/Cultivation-Depensation_Manuscript")
rm(list=ls())
library(dplyr)
library(vegan)
library(ggplot2)
library(glmulti) # It is strongly recommended that the user familiarizes themselves with the associated paper and user manual for this package in order to understand how it is used here. Calcagno, Vincent and Claire de Mazancourt. 2010. "glmulti: An R package for easy automated model selection with (generalized) linear models." Journal of Statistical Software May 2010, Volume 34, Issue 12. This is also cited in the main manuscript text.
library(ModelMetrics)

# see Open Research statement for details on accessing the data used below
dat=read.csv("qbAndCovariates.csv",stringsAsFactors = F) # q values from Sass et al. (2021)
fishComm=read.csv("Hansen_WI_fishCommunity.csv",stringsAsFactors = F) # fish community data used to conduct PCoA
riparian=read.csv("Hansen_WI_riparian.csv",stringsAsFactors = F) # riparian land cover data used to conduct PCA
watershed=read.csv("Hansen_WI_watershed.csv",stringsAsFactors = F) # watershed land cover data used to conduct PCA
ralphGDD=read.csv("NLDAS_mean_temperatures_WBIC.csv",stringsAsFactors = F) # growing degree days information from Winslow et al. (2017) <- DOUBLE CHECK THIS IS CORRECT
d1=read.csv('BASS_SEII_CPE.csv',stringsAsFactors = F) # largemouth bass relative abundance data from WDNR Fisheries Management Database

##### PCA/PCoA SECTION #####
# calculating number of species present for each lake and doing pca analysis on all lakes so I can use the same results for the prediction lakes
#creating one df with all pca columns
colnames(riparian)[2:9]=paste(rep("rip",length(colnames(riparian))-1),colnames(riparian)[2:9], sep="")
colnames(watershed)[3:10]=paste(rep("wshed",length(colnames(watershed))-2),colnames(watershed)[3:10], sep="")

bassCPEs=d1%>%
  group_by(WBIC)%>%
  summarise(bassCPE=mean(CPEmile,na.rm=T))

multivars=watershed%>%
  full_join(riparian, by="wbic")%>%
  full_join(fishComm, by="wbic")%>%
  full_join(bassCPEs, by=c("wbic"="WBIC"))%>%
  full_join(ralphGDD, by=c("wbic"="WBIC"))%>%
  select(wbic,starts_with("wshed"),starts_with("rip"),Panfish,Catfish,Largemouth.Bass,Musky,Northern.Pike,Rainbow.Smelt,Smallmouth.Bass,Sturgeon,Trout,Cisco,bassCPE,gdd5)
multivars=multivars[complete.cases(multivars),]
wshed=multivars[,c(1,2:9)]
comm=multivars[,c(1,18:27)]
rip=multivars[,c(1,10:17)]
wshed.dist=dist(wshed[-1])
rip.dist=dist(rip[-1])
wshed.pca=cmdscale(wshed.dist)
wshed.ef=envfit(wshed.pca,wshed)
rip.pca=cmdscale(rip.dist)
rip.ef=envfit(rip.pca,rip)

# plot(wshed.pca[,1],wshed.pca[,2]) # not necessary for analysis but can be used to visualize PCA results
# plot(wshed.ef,add=TRUE) # not necessary for analysis but can be used to visualize PCA results
# cor(cbind(wshed.pca,wshed))[1:2,] # not necessary for analysis but can be used to visualize PCA results
# notes on loadings for each axis below
# wshed axis 1 = +pasture +cultivated 
# wshed axis 2 = + Forest -Wetlands

# plot(rip.pca[,1],rip.pca[,2]) # not necessary for analysis but can be used to visualize PCA results
# plot(rip.ef,add=TRUE) # not necessary for analysis but can be used to visualize PCA results
# cor(cbind(rip.pca,rip))[1:2,] # not necessary for analysis but can be used to visualize PCA results
# notes on loadings for each axis below
# rip axis 1 = -Forest + Wetlands
# rip axis 2 = -Developed

comm.dist=vegdist(comm,method="bray")
comm.dist[is.na(comm.dist)]=0
comm.pcoa=cmdscale(comm.dist)
comm.ef=envfit(comm.pcoa,comm)
# plot(comm.pcoa[,1],comm.pcoa[,2]) # not necessary for analysis but can be used to visualize PCoA results
# text(comm.pcoa[,1],comm.pcoa[,2],1:nrow(comm)) # not necessary for analysis but can be used to visualize PCoA results
# plot(comm.ef,add=TRUE) # not necessary for analysis but can be used to visualize PCoA results
# new plot to zoom in on non-outlier lakes (lakes with no fish or only  cisco) # not necessary for analysis but can be used to visualize PCoA results
# plot(comm.pcoa[,1],comm.pcoa[,2],xlim=c(-0.4,0.4)) # not necessary for analysis but can be used to visualize PCoA results
# plot(comm.ef,add=TRUE) # not necessary for analysis but can be used to visualize PCoA results
# cor(cbind(comm.pcoa,comm))[2:3,]# not necessary for analysis but can be used to visualize PCoA results
# notes on loadings for each axis below
# fcomm axis 1 = + panfish + lmb
# fcomm axis 2 = + musky + smb

# linkng pca/pcoa results back to lake ID numbers for later use
comm.pcoa=as.data.frame(cbind(comm$wbic,comm.pcoa));colnames(comm.pcoa)=c("wbic","fcomm1","fcomm2")
rip.pca=as.data.frame(cbind(rip$wbic,rip.pca));colnames(rip.pca)=c("wbic","rip1","rip2")
wshed.pca=as.data.frame(cbind(wshed$wbic,wshed.pca));colnames(wshed.pca)=c("wbic","wshed1","wshed2")

#adding the pca/pcoa axes to the df
dat=dat%>%
  left_join(comm.pcoa, by=c("WBIC"="wbic"))%>%
  left_join(rip.pca, by=c("WBIC"="wbic"))%>%
  left_join(wshed.pca, by=c("WBIC"="wbic"))

#### SETTING UP DATA FRAME FOR INFERENCE SET OF LAKES ####
dat$overallMeanWLYDens=rowMeans(dat[,c(22,24)],na.rm = T)
dat$overallSDWLYDens=rowMeans(dat[,c(23,25)],na.rm = T)

# PICKING VARIABLES
#subsetting to our inference lakes
q=dat[dat$parm=="q",]
q=q[!is.na(q$bassMCPE),]
q$bassMCPE[q$bassMCPE==0]=1e-5 # one bass CPE is 0 and can't be logged, I made it 0.00001 so I could add that lake to the dataset, now 28 lakes from Sass et al. (2021) can be used
repl=11.157 # mean bass CPE sd from the inference lakes
replWLY=1.430889 # mean WLYsd from the inference lakes
# replacing NA bassCPE SDs with mean SD so that I can take random draws later on for bootstrapping
for(i in 1:nrow(q)) {
  q$bassSDCPE[i]=ifelse(is.na(q$bassSDCPE[i]),repl,q$bassSDCPE[i])
  q$overallSDWLYDens[i]=ifelse(is.na(q$overallSDWLYDens[i]),replWLY,q$overallSDWLYDens[i])
}
q$logBassCPE=log(q$bassMCPE)
q$logBassSDCPE=log(q$bassSDCPE)
q$logN=log(q$N)
vars=c("WBIC","Mean","fcomm1","wshed1","wshed2","rip1","rip2","gdd5","overallMeanWLYDens","logBassCPE","logN") # taking out fcomm2 since it doesn't contribute to explaining variation in the data
bootStrapvars=c("WBIC","Mean","fcomm1","wshed1","wshed2","rip1","rip2","gdd5","overallMeanWLYDens", "overallSDWLYDens","logBassCPE", "logBassSDCPE","logN")# taking out fcomm2 since it doesn't contribute to explaining variation in the data

# creating data frame with just inference lakes and just the 9 candidate predictor variables
q2=q[q$bassMCPE>0,colnames(q)%in%vars]
# one last filter to make sure to only use lakes that have observations for all 9 candidate predictor variables.
q3=q2[complete.cases(q2),-2]

#### GLMULTI FIT ####

# # first figuing out how many models would have to be fit to do an exhaustive search
# 
# 2^(ncol(q3)-1) # single effects only
# 
# (2^(ncol(q3)-1)^2) # interaction effects with everything
# 
# 2^((ncol(q3)-1)*2) # interaction effects with logBassCPE only
# 
# genetic algorithm approach, first working with single effects to get the mechanics down
# working with full model now

#first to get an estimate of how long 1 run might take
Sys.time()
gaF1=glmulti(y=Mean~ .,
             data = q3, level = 2,
             method = "g", marginality = F,
             crit = aicc, confsetsize = 100,
             popsize = 200, mutrate = 10e-3,
             sexrate = 0.15, imm = 0.2,
             deltaM = .01, deltaB = 0, conseq = 5,
             plotty = F,report = F,
             maxsize = 28)
Sys.time()
# takes ~ 1min per run
#quick look at the output
plot(gaF1,type = "s")
coef(gaF1)
# setting up loop to run 20 fits, should take about 20 mins
modNames=paste(rep("ga",20),1:20,sep = "_")
gsOut=list()
Sys.time()
for(i in 1:20){
  gsOut[[i]]=glmulti(y=Mean~., data = q3, level = 2,
                     method = "g", marginality = F,
                     crit = aicc, confsetsize = 100,
                     popsize = 200, mutrate = 10e-3,
                     sexrate = 0.15, imm = 0.2,
                     deltaM = .01, deltaB = 0, conseq = 5,
                     plotty = F,report = F,
                     maxsize = 28,
                     name=modNames[i])
}
Sys.time()

#no log N - This is a check to make sure the number of stock-recruitment observations doesn't explain variation in q. If the results of this genetic algorithm differed from 'gsOut' above then we would be worried about sampling bias effect on our results. Turns out this is not the case and whether or not logN is included doesn't impact the results.
#make sure df has only vars you want first
q5=q3[,-10]
modNames=paste(rep("ga",20),1:20,sep = "_")
gsOutN=list()
Sys.time()
for(i in 1:20){
  gsOutN[[i]]=glmulti(y=Mean~., data = q5, level = 2,
                     method = "g", marginality = F,
                     crit = aicc, confsetsize = 100,
                     popsize = 200, mutrate = 10e-3,
                     sexrate = 0.15, imm = 0.2,
                     deltaM = .01, deltaB = 0, conseq = 5,
                     plotty = F,report = F,
                     maxsize = 28,
                     name=modNames[i])
}
Sys.time()

#saving the output so that I'm always working with the same set of results. Because of the stochasticity in the genetic algorithm each one is slightly different (hence the 20 replicate runs of the algorithm).
saveRDS(gsOut,file = "gsOut_3.29.22.RData")
saveRDS(gsOutN,file = "gsOut_noLogN_3.29.22.RData")

#### WORKING WITH MODEL OUTPUT ####
gsOut=readRDS("gsOut_3.29.22.RData")
gsOutN=readRDS("gsOut_noLogN_3.29.22.RData")
# spot checking some model output to see if anything has gone wrong - all looks good. If something has gone wrong, like if the model is having a hard time identifying important predictor variables then large numbers (likely more than a dozen or so) of points on the following plots will fall below the red line. This signals that a lot of models are performing about the same instead of a clear subset separating themselves from the pack.
plot(gsOut[[1]],type = "p")
plot(gsOut[[5]],type = "p")
plot(gsOut[[10]],type = "p")
plot(gsOut[[15]],type = "p")
plot(gsOut[[20]],type = "p")

plot(gsOut[[1]],type = "s")
plot(gsOut[[5]],type = "s")
plot(gsOut[[10]],type = "s")
plot(gsOut[[15]],type = "s")
plot(gsOut[[20]],type = "s")

#creating a consensus model of the top 100 models from across the 20 replicate runs of the genetic algorithm
con=consensus(gsOut,confsetsize = 100)
conN=consensus(gsOutN, confsetsize = 100)

plot(con, type = "p")
plot(1:100,con@crits,pch=16,ylab="AICc",xlab = "Models")
abline(h=(min(con@crits)+2),col="red")
plot(con, type = "s")
plot(conN, type = "p")
plot(conN, type = "s")
plot(1:100,conN@crits,pch=16,ylab="AICc",xlab = "Models")
abline(h=(min(conN@crits)+2),col="red")

# looking at the important variables from each consensus set of models
parms=as.data.frame(coef(con))
row.names(parms[parms$Importance>0.80,]) # those with > 80% importance, not used in manuscript
wts=weightable(con)
wts$model[wts$aicc==min(wts$aicc)] # pulling out the model with the minimum AICc value

parmsN=as.data.frame(coef(conN))
row.names(parmsN[parmsN$Importance>0.80,]) # those with > 80% importance, not used in manuscript
wtsN=weightable(conN)
wtsN$model[wtsN$aicc==min(wtsN$aicc)]# pulling out the model with the minimum AICc value

# the fact that the models identified on lines 217 and 222 are the same further demonstrates the lack of influence that number of stock-recruitment observations (logN) has on the outcome which is good reassurance that any sampling bias among the 28 inference lakes that might exists is not effecting the results.

#### BUILDING A FEW MODELS TO CROSS VALIDATE FROM GENETIC ALGORITHM  ####

# this is simple leave-one-out cross-validation (LOOCV)
# models within the 2 aicc cutoff of the best model suggested by Calcagno & de Mazancourt (2010) glmuti paper
bestAIC=min(wts$aicc)
mForms=wts$model[(wts$aicc-bestAIC)<2] # model forms to cross validate
topMods=list()
# fitting those top models to the inference set of lakes 
for(i in 1:length(mForms)){
  topMods[[i]]=glm(mForms[i],data = q3)
}
names(topMods)=paste(rep("m",length(mForms)),1:length(mForms),sep = "")

plot(topMods$m1$y,topMods$m1$fitted.values,pch=16,ylim = c(0,1.5),xlim=c(0,1.5), xlab="Observed q values", ylab = "Model predicted q value")
abline(b=1,a=0,lty=2) # 1:1 line to see how well the example model (m1=first model form to cross validate) is predicts the data it's fit to.
for(i in 2:length(mForms)){ # adding points from the other model form fits to the same plot, this is just a vizualizatio exercise to make sure all the code is working propoerly
  points(topMods[[i]]$y,topMods[[i]]$fitted.values,pch=i)
}
rmsesFullFit=data.frame(mod=names(topMods),rmse=numeric(length = length(mForms)))
for(i in 1:length(mForms)){
  rmsesFullFit$rmse[i]=rmse(topMods[[i]]$y,topMods[[i]]$fitted.values)
}

# LOOCV
loocvTMods=list()
for(j in 1:length(mForms)){
  loocvTMods[[j]]=data.frame(heldOut=numeric(nrow(q3)),pred=numeric(nrow(q3)))
  for(i in 1:nrow(q3)){
    tempQ=q3[-i,]
    tempFit=glm(mForms[j],data=tempQ)
    loocvTMods[[j]]$pred[i]=predict(tempFit,q3[i,])
    loocvTMods[[j]]$heldOut[i]=q3$Mean[i]
  }
}
#calculate root-mean-squared-error for the LOOCV points
rmseLOOCV=data.frame(mod=names(topMods),rmse=numeric(length = length(mForms)))
for(i in 1:length(mForms)){
  rmseLOOCV$rmse[i]=rmse(loocvTMods[[i]]$heldOut,loocvTMods[[i]]$pred)
}

# k-fold cross-validation using 5 groups of data to iteratively hold out.
# split data into 5 groups - these are used again below for the 'no logN' cross-validation to again check that any sampling bias that could exist in our inference set of lakes has no effect on our results.
gNums=c(rep(1,5),rep(2,5),rep(3,6),rep(4,6),rep(5,6)) #group labels to randomly draw from and assign to rows of q3
overallKFCV=list(q3,q3,q3,q3,q3,q3,q3,q3,q3,q3)
for(i in 1:10){
  overallKFCV[[i]]$kfoldGroup=sample(gNums,size = nrow(q3),replace = F)
}
#looping to do k fold cross validation 10 times to make sure the random group assignment isn't having a large effect on rmse
overallRMSE=data.frame(mod=character(),rmse=numeric())
for(k in 1:length(overallKFCV)){
  kFCV=list()
  for(j in 1:length(mForms)){
    tempDF=data.frame(groupLabel=overallKFCV[[k]][["kfoldGroup"]],heldOut=q3$Mean,pred=numeric(nrow(q3)))
    for(i in 1:5){
      tempQ=overallKFCV[[k]][overallKFCV[[k]]$kfoldGroup!=i,] 
      tempFit=glm(mForms[j],data = tempQ)
      tempDF$pred[tempDF$groupLabel==i]=predict(tempFit,overallKFCV[[k]][overallKFCV[[k]]$kfoldGroup==i,])
    }
    kFCV[[j]]=list(mForms[j],tempDF)
  }
  rmseKFCV=data.frame(mod=names(topMods),rmse=numeric(length = length(mForms)))
  for(i in 1:length(mForms)){
    rmseKFCV$rmse[i]=rmse(kFCV[[i]][[2]][["heldOut"]],kFCV[[i]][[2]][["pred"]])
  }
  overallRMSE=rbind(overallRMSE,rmseKFCV)
}

avgRMSE=overallRMSE%>%
  group_by(mod)%>%
  summarise(meanRMSE=mean(rmse))

#### BUILDING A FEW MODELS TO CROSS VALIDATE FROM FIT W/OUT LOG_N ####
# this is simple leave-one-out cross-validation
# models within the 2 aicc cutoff of the best model suggested by glmulti paper
bestAICN=min(wtsN$aicc)
mForms=wtsN$model[(wtsN$aicc-bestAICN)<2]
topN=list()
# fits
for(i in 1:length(mForms)){
  topN[[i]]=glm(mForms[i],data = q3)
}
names(topN)=paste(rep("m",length(mForms)),1:length(mForms),sep = "")

# calculate root-mean-squared-error
rmsesFullFitN=data.frame(mod=names(topN),rmse=numeric(length = length(mForms)))
for(i in 1:length(mForms)){
  rmsesFullFitN$rmse[i]=rmse(topN[[i]]$y,topN[[i]]$fitted.values)
}

# LOOCV
loocvTN=list()
for(j in 1:length(mForms)){
  loocvTN[[j]]=data.frame(heldOut=numeric(nrow(q3)),pred=numeric(nrow(q3)))
  for(i in 1:nrow(q3)){
    tempQ=q3[-i,]
    tempFit=glm(mForms[j],data=tempQ)
    loocvTN[[j]]$pred[i]=predict(tempFit,q3[i,])
    loocvTN[[j]]$heldOut[i]=q3$Mean[i]
  }
}
rmseLOOCVN=data.frame(mod=names(topN),rmse=numeric(length = length(mForms)))
for(i in 1:length(mForms)){
  rmseLOOCVN$rmse[i]=rmse(loocvTN[[i]]$heldOut,loocvTN[[i]]$pred)
}

# 5-fold k-fold cross-validation
overallKFCV.N=list(q3,q3,q3,q3,q3,q3,q3,q3,q3,q3)
for(i in 1:10){
  overallKFCV.N[[i]]$kfoldGroup=sample(gNums,size = nrow(q3),replace = F)
}
#looping to do k fold cross validation 10 times to make sure the random group assignment isn't having a large effect on rmse
overallRMSE.N=data.frame(mod=character(),rmse=numeric())
for(k in 1:length(overallKFCV.N)){
  kFCV.N=list()
  for(j in 1:length(mForms)){
    tempDF=data.frame(groupLabel=overallKFCV.N[[k]][["kfoldGroup"]],heldOut=q3$Mean,pred=numeric(nrow(q3)))
    for(i in 1:5){
      tempQ=overallKFCV.N[[k]][overallKFCV.N[[k]]$kfoldGroup!=i,] 
      tempFit=glm(mForms[j],data = tempQ)
      tempDF$pred[tempDF$groupLabel==i]=predict(tempFit,overallKFCV.N[[k]][overallKFCV.N[[k]]$kfoldGroup==i,])
    }
    kFCV.N[[j]]=list(mForms[j],tempDF)
  }
  rmseKFCV.N=data.frame(mod=names(topN),rmse=numeric(length = length(mForms)))
  for(i in 1:length(mForms)){
    rmseKFCV.N$rmse[i]=rmse(kFCV.N[[i]][[2]][["heldOut"]],kFCV.N[[i]][[2]][["pred"]])
  }
  overallRMSE.N=rbind(overallRMSE.N,rmseKFCV.N)
}

avgRMSE.N=overallRMSE.N%>%
  group_by(mod)%>%
  summarise(meanRMSE=mean(rmse))


##### CONCLUDING NOTE #####

# The end result of this script is that we've identified the model structure that is best supported via model fit and k-fold cross validation. This model structure is presented in table 1 and the results of this cross validation are also presented in the Appendix S1. From this script the chosen model structure is used in the next script (Dassow.et.al._tablesFigures.R) to create tables, figures, and run other analyses presented in the manuscript. The decision to separate the model fitting and selection process here from the generation of tables, figures, and other analyses is so that this analysis and model selection code doesn't have to be re-run each time someone wants to look at the analyses, figures, and tables presented in the next script. This should make it easier for someone to decide which part of the whole analysis (model selection & validation vs. statistical analysis and management application) they want to explore.