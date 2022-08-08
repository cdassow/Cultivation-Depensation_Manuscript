# CJD 2.28.22
# same script as 'glmultiModelSelection.R' script but using pcoa's that don't include walleye in the fish community data

##### CHUNK OF DATA PREP BEFORE COPYING VIGNETTE #####
setwd("~/NotreDame/NDstuff/Dissertation/5")
rm(list=ls())
library(dplyr)
library(vegan)
library(ggplot2)
library(glmulti)
library(ModelMetrics)

# see Hansen_WI_Co.variates.xlsx for metadata on covariates
# see Sassetal._Fisheries_SupplementalTables_Final.xlsx for meta data on depensation parms
dat=read.csv("qbAndCovariates.csv",stringsAsFactors = F)
fishComm=read.csv("Hansen_WI_fishCommunity.csv",stringsAsFactors = F)
riparian=read.csv("Hansen_WI_riparian.csv",stringsAsFactors = F)
watershed=read.csv("Hansen_WI_watershed.csv",stringsAsFactors = F)
ralphGDD=read.csv("~/NotreDame/NDstuff/Dissertation/4/Conversion_file_MedResTemp_to_WBIC/Conversion_file_MedResTemp_to_WBIC/output/NLDAS_mean_temperatures_WBIC.csv",stringsAsFactors = F)
d1=read.csv('~/NotreDame/NDstuff/CNH/BASS_SEII_CPE.csv',stringsAsFactors = F)

# calculating number of species present for each lake and doing pca analysis on all lakes so I can use the same resutls for the prediction lakes
#dat$sppRichness=rowSums(dat[,c(40:50)])
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
# plot(wshed.pca[,1],wshed.pca[,2])
# plot(wshed.ef,add=TRUE)
# cor(cbind(wshed.pca,wshed))[1:2,]
# wshed axis 1 = +pasture +cultivated 
# wshed axis 2 = + Forest -Wetlands

# plot(rip.pca[,1],rip.pca[,2])
# plot(rip.ef,add=TRUE)
# cor(cbind(rip.pca,rip))[1:2,]
# rip axis 1 = -Forest + Wetlands
# rip axis 2 = -Developed

comm.dist=vegdist(comm,method="bray")
comm.dist[is.na(comm.dist)]=0
comm.pcoa=cmdscale(comm.dist)
comm.ef=envfit(comm.pcoa,comm)
# plot(comm.pcoa[,1],comm.pcoa[,2])
# text(comm.pcoa[,1],comm.pcoa[,2],1:nrow(comm))
# plot(comm.ef,add=TRUE)
# new plot to zoom in on non-outlier lakes (lakes with no fish or only  cisco)
# plot(comm.pcoa[,1],comm.pcoa[,2],xlim=c(-0.4,0.4))
# plot(comm.ef,add=TRUE)
# cor(cbind(comm.pcoa,comm))[2:3,]
# fcomm axis 1 = + panfish + lmb
# fcomm axis 2 = + musky + smb


comm.pcoa=as.data.frame(cbind(comm$wbic,comm.pcoa));colnames(comm.pcoa)=c("wbic","fcomm1","fcomm2")
rip.pca=as.data.frame(cbind(rip$wbic,rip.pca));colnames(rip.pca)=c("wbic","rip1","rip2")
wshed.pca=as.data.frame(cbind(wshed$wbic,wshed.pca));colnames(wshed.pca)=c("wbic","wshed1","wshed2")

#adding the useful axes to the df
dat=dat%>%
  left_join(comm.pcoa, by=c("WBIC"="wbic"))%>%
  left_join(rip.pca, by=c("WBIC"="wbic"))%>%
  left_join(wshed.pca, by=c("WBIC"="wbic"))

#### FORWARD MODEL SELECTION ####
dat$overallMeanWLYDens=rowMeans(dat[,c(22,24)],na.rm = T)
dat$overallSDWLYDens=rowMeans(dat[,c(23,25)],na.rm = T)
# replacing NA bassCPE with mean SD so that I can take random draws later on for bootstrapping

# PICKING VARIABLES
#subsetting to our inference lakes
q=dat[dat$parm=="q",]
q=q[!is.na(q$bassMCPE),]
q$bassMCPE[q$bassMCPE==0]=1e-5 # one bass CPE is 0 and can't be logged, I made it 0.00001 so I could add an extra point, now 28 can be used
repl=11.157 # mean bass CPE sd from the inference lakes
replWLY=1.430889 # mean WLYsd from the inference lakes
for(i in 1:nrow(q)) {
  q$bassSDCPE[i]=ifelse(is.na(q$bassSDCPE[i]),repl,q$bassSDCPE[i])
  q$overallSDWLYDens[i]=ifelse(is.na(q$overallSDWLYDens[i]),replWLY,q$overallSDWLYDens[i])
}
q$logBassCPE=log(q$bassMCPE)
q$logBassSDCPE=log(q$bassSDCPE)
q$logN=log(q$N)
vars=c("WBIC","Mean","fcomm1","wshed1","wshed2","rip1","rip2","gdd5","overallMeanWLYDens","logBassCPE","logN") # taking out fcomm2 since it doesn't really contribute to explaining variation in the data
bootStrapvars=c("WBIC","Mean","fcomm1","wshed1","wshed2","rip1","rip2","gdd5","overallMeanWLYDens", "overallSDWLYDens","logBassCPE", "logBassSDCPE","logN")# taking out fcomm2 since it doesn't really contribute to explaining variation in the data


q2=q[q$bassMCPE>0,colnames(q)%in%vars]
q3=q2[complete.cases(q2),-2]

#### GLMULTI FIT ####

# MARGINAL = T to make sure that if an interaction effect is included then both main effects are included. Doing this for now but there may be good reasons not to do this. It does help reduce the number of models that need to be searched

# METHOD ="g" so improve computation time because of the large number of models I have to search through Genetic algorithm approach will help do this more efficiently
  # popsize - the number of candidate models to evaluate at each generation
  # mutrate - rate at which model parms are turned on or off (locus) randomly between generations
  # sexrate - proportion of new models in the next generation that are produced by combining two randomly chose models from the current generation (model probability of being chosen is a function of their ability to fit the data)
  # imm - new randomly generated model with each locus turned on or off. This prevents the algorithm from getting stuck in local minima allowing for better global convergence.
  # Stopping rules for the algorithm
    # deltaB - minimum improvement in the IC for the best model
    # deltaM - minimum improvement in the average IC for the population of models
    # conseq - number of consecutive failed improvements that has to happen before stopping

# CRIT = "bic" bayesian information criterion. The IC method for finding best model doesn't allow for cross validation RMSE as the scoring criteria for finding best model so I'll have to do that manually with the model average model that I end up with. 



# # first figuing out how many models would have to be fit to do an exhaustive search
# 
# 2^(ncol(q3)-1) # single effects only
# 
# (2^(ncol(q3)-1)^2) # interaction effects with everything
# 
# 2^((ncol(q3)-1)*2) # interaction effects with logBassCPE only
# 
# # genetic algorithm approach, first working with single effects to get the mechanics down
# 
# #number of models to search
# nMod=glmulti(Mean~logBassCPE+gdd5+fcomm1+fcomm2,data=q3,
#             method = "d", marginality = T)
# 
# #genetic algorithm
# Sys.time()
# gaS1=glmulti(Mean~logBassCPE+gdd5+fcomm1+fcomm2,data=q3,
#            method = "g", marginality = F,
#            crit = aicc, confsetsize = 500, 
#            popsize = 500, mutrate = 10^-3,
#            sexrate = 0.1, imm = 0.6,
#            deltaM = 1, deltaB = 2, conseq = 10,
#            plotty = F,report = F)
# Sys.time()
# gaS2=glmulti(Mean~logBassCPE+gdd5+fcomm1+fcomm2,data=q3,
#              method = "g", marginality = F,
#              crit = aicc, confsetsize = 100, 
#              popsize = 100, mutrate = 10^-3,
#              sexrate = 0.1, imm = 0.6,
#              deltaM = 1, deltaB = 2, conseq = 10)
# gaS3=glmulti(Mean~logBassCPE+gdd5+fcomm1+fcomm2,data=q3,
#              method = "g", marginality = F,
#              crit = aicc, confsetsize = 100, 
#              popsize = 100, mutrate = 10^-3,
#              sexrate = 0.1, imm = 0.6,
#              deltaM = 1, deltaB = 2, conseq = 10)
# gaCon=consensus(list(gaS1,gaS2,gaS3),confetsize=10)
# plot(gaCon,type = "s")
# #exhaustive search
# exS=glmulti(Mean~logBassCPE+gdd5+fcomm1+fcomm2,data=q3,
#             method = "h", marginality = F,
#             crit = bic, confsetsize = 100)
# plot(exS,type = "s")

# still working with small model to get an idea of how differnt parm values effect convergence
# pops=seq(10,200,length.out=20)
# muts=seq(10e-4,10e-2,length.out=20)
# sexR=seq(.01, .5, length.out=20)
# immR=seq(0.1,.5,length.out=20)
# db=seq(0.01, 5,length.out=20)
# dm=seq(0.01,5,length.out=20)
# gsEx=list()
# Sys.time()
# for(i in 1:20){
#   gsEx[[i]]=glmulti(y=Mean~logBassCPE+gdd5+logN+rip1, data = q3, level = 2,
#                      method = "g", marginality = F,
#                      crit = aicc, confsetsize = 500, 
#                      popsize = 100, mutrate = 10e-3,
#                      sexrate = 0.3, imm = 10e-3,
#                      deltaM = dm[i], deltaB = .5, conseq = 5,
#                      plotty = F,report = F,
#                      maxsize = 28)
# }
# Sys.time()
# 

# tried using the wrapper function he provides and it doesn't seem to work so I'm moving on for now
# myglm=function(y, data) {
#   # we get the terms for the formula
#   if (missing(data)) data<-environment(y)
#   termz = terms(y,data=data)
#   # we get the order of all terms
#   orderz = attr(termz,"order")
#   # we get all pairwise interactions, if any. Otherwise the formula is okay
#   intz = which(orderz==2)
#   # we locate the row corresponding to the desired variable (as it is not constant)
#   # HERE I WANT VARIABLE M (change accordingly)
#   index=which(dimnames(attr(termz,"factors"))[[1]] == "logBassCPE")
#   if (length(index)>0) { # the desired effect must be present
#     if(length(intz)>0) {
#       # we simply test that all interactions include the desired effect
#       # otherwise we return the crappy null model
#       if (min(attr(termz,"factors")[index,intz])==0) return(nullos)
#     }
#   } else return(nullos);
#   # if all is good we just call glm as usual
#   return(glm(formula=y, data=data)) 
# }
# nullos=glm(Mean~1,data=q3)
# Sys.time()
# gaF1=glmulti(y=Mean~ .,
#              data = q3, level = 2,
#              method = "g", marginality = F,
#              crit = aicc, confsetsize = 100, 
#              popsize = 200, mutrate = 10e-3,
#              sexrate = 0.15, imm = 0.2,
#              deltaM = .01, deltaB = 0, conseq = 5,
#              plotty = F,report = F,
#              maxsize = 28,fitfunction = myglm)
# Sys.time()
# working with full model now

#first to get an estimate of how long 1 run might take
qmodel=glm(Mean~.+.:logBassCPE,data = q3)
qmodel2=glm(Mean~.*., data = q3)
qmodel3=glm(Mean~ gdd5 + fcomm1 + fcomm2 + wshed1 + wshed2 + rip1 + rip2 + overallMeanWLYDens + logBassCPE + gdd5:logBassCPE + fcomm1:logBassCPE + fcomm2:logBassCPE + wshed1:logBassCPE + wshed2:logBassCPE + rip1:logBassCPE + rip2:logBassCPE + overallMeanWLYDens:logBassCPE - gdd5:fcomm1 - gdd5:fcomm2 - gdd5:wshed1 - gdd5:wshed2 - gdd5:rip1 - gdd5:rip2 - gdd5:overallMeanWLYDens - fcomm1:fcomm2 - fcomm1:wshed1 - fcomm1:wshed2 - fcomm1:rip1 - fcomm1:rip2 - fcomm1:overallMeanWLYDens - fcomm2:wshed1 - fcomm2:wshed2 - fcomm2:rip1 - fcomm2:rip2 - fcomm2:overallMeanWLYDens - wshed1:wshed2 - wshed1:rip1 - wshed1:rip2 - wshed1:overallMeanWLYDens - wshed2:rip1 - wshed2:rip2 - wshed2:overallMeanWLYDens - rip1:rip2 - rip1:overallMeanWLYDens - rip2:overallMeanWLYDens,
            data = q3)
pred1=names(qmodel$coefficients)[-1]
pred2=names(coef(qmodel2))[-1]
Pred=colnames(q3)[-1]
excl=pred2[!pred2%in%pred1]
preds=paste(pred1,collapse = " + ")
form=formula(paste("Mean",preds, sep = "~"))

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

#no log N
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

#no log N and now WLY
#make sure df has only vars you want first
q4=q3[,-c(8,10)]
modNames=paste(rep("ga",20),1:20,sep = "_")
gsOutNW=list()
Sys.time()
for(i in 1:20){
  gsOutNW[[i]]=glmulti(y=Mean~., data = q4, level = 2,
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
#saving the output
saveRDS(gsOut,file = "gsOut_3.29.22.RData")
saveRDS(gsOutN,file = "gsOut_noLogN_3.29.22.RData")
saveRDS(gsOutNW,file = "gsOut_noLogNWLY_3.29.22.RData")

#### WORKING WITH MODEL OUTPUT ####
gsOut=readRDS("gsOut_5.20.21.RData")
gsOutN=readRDS("gsOut_noLogN_3.29.22.RData")
gsOutNW=readRDS("gsOut_noLogNWLY_5.20.21.RData")
# spot checking some model output to see if anything has gone wrong
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
#creating a consensus model

con=consensus(gsOut,confsetsize = 100)
conN=consensus(gsOutN, confsetsize = 100)
conNW=consensus(gsOutNW, confsetsize = 100)

plot(con, type = "p")
plot(1:100,con@crits,pch=16,ylab="AICc",xlab = "Models")
abline(h=(min(con@crits)+2),col="red")
plot(con, type = "s")
plot(conN, type = "p")
plot(conN, type = "s")
plot(1:100,conN@crits,pch=16,ylab="AICc",xlab = "Models")
abline(h=(min(conN@crits)+2),col="red")
plot(conNW, type = "p")
plot(conNW, type = "s")

parms=as.data.frame(coef(con))
row.names(parms[parms$Importance>0.80,])
wts=weightable(con)
wts$model[wts$aicc==min(wts$aicc)]

parmsN=as.data.frame(coef(conN))
row.names(parmsN[parmsN$Importance>0.80,])
wtsN=weightable(conN)
wtsN$model[wtsN$aicc==min(wtsN$aicc)]

parmsNW=as.data.frame(coef(conNW))
row.names(parmsNW[parmsNW$Importance>0.80,])
wtsNW=weightable(conNW)
wtsNW$model[wtsNW$aicc==min(wtsNW$aicc)]

#### BUILDING A FEW MODELS TO CROSS VALIDATE FROM FIT W/ ALL 15 COVARIATES ####

# best model based on aicc and model weight
best=glm(wts$model[wts$aicc==min(wts$aicc)],data = q3)
plot(best$fitted.values,best$y,pch=16,ylim=c(0,1.5),xlim = c(0,1.5))
abline(b=1,a=0,lty=2)
rmse(best$y,best$fitted.values) #rmse for model fit to whole data set

#   LOOCV 
outbest=data.frame(heldOut=numeric(nrow(q3)),pred=numeric(nrow(q3)))
for(i in 1:nrow(q3)){
  tempQ=q3[-i,]
  tempFit=glm(best$formula,data=tempQ)
  outbest$pred[i]=predict(tempFit,q3[i,])
  outbest$heldOut[i]=q3$Mean[i]
}
rmse.best=rmse(outbest$heldOut,outbest$pred) #rmse for loocv
plot(outbest$heldOut,outbest$pred,pch=16,ylim = c(0,1.5),xlim=c(0,1.5))
abline(b=1,a=0,lty=2)

# model using only parms with >80% importance
vars=paste(row.names(parms[parms$Importance>0.80,])[-2],collapse = "+")
m80Form=as.formula(paste("Mean",vars,sep = "~"))
m80=glm(m80Form,data = q3)
plot(m80$fitted.values,m80$y,pch=16,ylim = c(0,1.5),xlim = c(0,1.5))
abline(b=1,a=0,lty=2)
rmse(m80$y,m80$fitted.values)


#   LOOCV 
out80=data.frame(heldOut=numeric(nrow(q3)),pred=numeric(nrow(q3)))
for(i in 1:nrow(q3)){
  tempQ=q3[-i,]
  tempFit=glm(m80$formula,data=tempQ)
  out80$pred[i]=predict(tempFit,q3[i,])
  out80$heldOut[i]=q3$Mean[i]
}
rmse.80=rmse(out80$heldOut,out80$pred) #rmse for loocv
plot(out80$heldOut,out80$pred,pch=16,ylim = c(0,1.5),xlim=c(0,1.5))
abline(b=1,a=0,lty=2)


# models within the 2 aicc cutoff of the best model suggested by glmulti paper
bestAIC=min(wts$aicc)
mForms=wts$model[(wts$aicc-bestAIC)<2]
topMods=list()
# fits
for(i in 1:length(mForms)){
  topMods[[i]]=glm(mForms[i],data = q3)
}
names(topMods)=paste(rep("m",length(mForms)),1:length(mForms),sep = "")

plot(topMods$m1$y,topMods$m1$fitted.values,pch=16,ylim = c(0,1.5),xlim=c(0,1.5))
abline(b=1,a=0,lty=2)
for(i in 2:length(mForms)){
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
rmseLOOCV=data.frame(mod=names(topMods),rmse=numeric(length = length(mForms)))
for(i in 1:length(mForms)){
  rmseLOOCV$rmse[i]=rmse(loocvTMods[[i]]$heldOut,loocvTMods[[i]]$pred)
}

# 5-fold k-fold cross-validation
# split data into 5 groups - these are used again below for the N & NW kfold cross-validation
gNums=c(rep(1,5),rep(2,5),rep(3,6),rep(4,6),rep(5,6)) #group labels to randomly draw and assign to rows of q3
ss=sample(gNums,size = nrow(q3),replace = F) #practice sample, no longer used
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

# best model based on aicc and model weight
bestN=glm(wtsN$model[wtsN$aicc==min(wtsN$aicc)],data = q3)
# plot(bestN$fitted.values,bestN$y,pch=16,ylim=c(0,1.5),xlim = c(0,1.5))
# abline(b=1,a=0,lty=2)
rmse(bestN$y,bestN$fitted.values) #rmse for model fit to whole data set

#   LOOCV 
outbestN=data.frame(heldOut=numeric(nrow(q3)),pred=numeric(nrow(q3)))
for(i in 1:nrow(q3)){
  tempQ=q3[-i,]
  tempFit=glm(bestN$formula,data=tempQ)
  outbestN$pred[i]=predict(tempFit,q3[i,])
  outbestN$heldOut[i]=q3$Mean[i]
}
rmse.bestN=rmse(outbestN$heldOut,outbestN$pred) #rmse for loocv
# plot(outbestN$heldOut,outbestN$pred,pch=16,ylim = c(0,1.5),xlim=c(0,1.5))
# abline(b=1,a=0,lty=2)


# model using only parms with >80% importance
vars=paste(row.names(parmsN[parmsN$Importance>0.80,])[-4],collapse = "+")
m80Form=as.formula(paste("Mean",vars,sep = "~"))
m80N=glm(m80Form,data = q3)
# plot(m80N$fitted.values,m80N$y,pch=16,ylim = c(0,1.5),xlim = c(0,1.5))
# abline(b=1,a=0,lty=2)
rmse(m80N$y,m80N$fitted.values)


#   LOOCV 
out80N=data.frame(heldOut=numeric(nrow(q3)),pred=numeric(nrow(q3)))
for(i in 1:nrow(q3)){
  tempQ=q3[-i,]
  tempFit=glm(m80N$formula,data=tempQ)
  out80N$pred[i]=predict(tempFit,q3[i,])
  out80N$heldOut[i]=q3$Mean[i]
}
rmse.80N=rmse(out80N$heldOut,out80N$pred) #rmse for loocv
# plot(out80N$heldOut,out80N$pred,pch=16,ylim = c(0,1.5),xlim=c(0,1.5))
# abline(b=1,a=0,lty=2)


# models within the 2 aicc cutoff of the best model suggested by glmulti paper
bestAICN=min(wtsN$aicc)
mForms=wtsN$model[(wtsN$aicc-bestAICN)<2]
topN=list()
# fits
for(i in 1:length(mForms)){
  topN[[i]]=glm(mForms[i],data = q3)
}
names(topN)=paste(rep("m",length(mForms)),1:length(mForms),sep = "")

# plot(topN$m1$y,topN$m1$fitted.values,pch=16,ylim = c(0,1.5),xlim=c(0,1.5))
# abline(b=1,a=0,lty=2)
# for(i in 2:length(mForms)){
#   points(top6[[i]]$y,top6[[i]]$fitted.values,pch=i+1)
# }
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
# split data into 5 groups - these are used again below for the NW kfold cross-validation
# gNums=c(rep(1,5),rep(2,5),rep(3,6),rep(4,6),rep(5,6)) #group labels to randomly draw and assign to rows of q3
# ss=sample(gNums,size = nrow(q3),replace = F) #practice sample, no longer used
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

#### BUILDING A FEW MODELS TO CROSS VALIDATE FROM FIT W/OUT LOG_N AND WLY_DENS ####

# best model based on aicc and model weight
bestNW=glm(wtsNW$model[wtsNW$aicc==min(wtsNW$aicc)],data = q4)
# plot(bestN$fitted.values,bestN$y,pch=16,ylim=c(0,1.5),xlim = c(0,1.5))
# abline(b=1,a=0,lty=2)
rmse(bestN$y,bestN$fitted.values) #rmse for model fit to whole data set

#   LOOCV 
outbestNW=data.frame(heldOut=numeric(nrow(q3)),pred=numeric(nrow(q4)))
for(i in 1:nrow(q4)){
  tempQ=q4[-i,]
  tempFit=glm(bestNW$formula,data=tempQ)
  outbestNW$pred[i]=predict(tempFit,q4[i,])
  outbestNW$heldOut[i]=q4$Mean[i]
}
rmse.bestNW=rmse(outbestNW$heldOut,outbestNW$pred) #rmse for loocv
# plot(outbestN$heldOut,outbestN$pred,pch=16,ylim = c(0,1.5),xlim=c(0,1.5))
# abline(b=1,a=0,lty=2)


# model using only parms with >80% importance
vars=paste(row.names(parmsNW[parmsNW$Importance>0.80,])[-4],collapse = "+")
m80Form=as.formula(paste("Mean",vars,sep = "~"))
m80NW=glm(m80Form,data = q4)
# plot(m80N$fitted.values,m80N$y,pch=16,ylim = c(0,1.5),xlim = c(0,1.5))
# abline(b=1,a=0,lty=2)
rmse(m80NW$y,m80NW$fitted.values)


#   LOOCV 
out80NW=data.frame(heldOut=numeric(nrow(q4)),pred=numeric(nrow(q4)))
for(i in 1:nrow(q4)){
  tempQ=q4[-i,]
  tempFit=glm(m80NW$formula,data=tempQ)
  out80NW$pred[i]=predict(tempFit,q4[i,])
  out80NW$heldOut[i]=q4$Mean[i]
}
rmse.80NW=rmse(out80NW$heldOut,out80NW$pred) #rmse for loocv
# plot(out80N$heldOut,out80N$pred,pch=16,ylim = c(0,1.5),xlim=c(0,1.5))
# abline(b=1,a=0,lty=2)


# models within the 2 aicc cutoff of the best model suggested by glmulti paper
bestAICNW=min(wtsNW$aicc)
mForms=wtsNW$model[(wtsNW$aicc-bestAICNW)<2]
topNw=list()
# fits
for(i in 1:length(mForms)){
  topNw[[i]]=glm(mForms[i],data = q3)
}
names(topNw)=paste(rep("m",length(mForms)), 1:length(mForms), sep="")

plot(topNw$m1$y,topNw$m1$fitted.values,pch=16,ylim = c(0,1.5),xlim=c(0,1.5))
abline(b=1,a=0,lty=2)
for(i in 2:length(mForms)){
  points(topNw[[i]]$y,topNw[[i]]$fitted.values,pch=i+1)
}
points(topNw[[7]]$y,topNw[[7]]$fitted.values,pch=16,col="red")
rmsesFullFitNw=data.frame(mod=names(topNw),rmse=numeric(length = length(mForms)))
for(i in 1:length(mForms)){
  rmsesFullFitNw$rmse[i]=rmse(topNw[[i]]$y,topNw[[i]]$fitted.values)
}

# LOOCV
loocvTNW=list()
for(j in 1:length(mForms)){
  loocvTNW[[j]]=data.frame(heldOut=numeric(nrow(q4)),pred=numeric(nrow(q4)))
  for(i in 1:nrow(q4)){
    tempQ=q4[-i,]
    tempFit=glm(mForms[j],data=tempQ)
    loocvTNW[[j]]$pred[i]=predict(tempFit,q4[i,])
    loocvTNW[[j]]$heldOut[i]=q4$Mean[i]
  }
}
rmseLOOCVNW=data.frame(mod=names(topNw),rmse=numeric(length = length(mForms)))
for(i in 1:length(mForms)){
  rmseLOOCVNW$rmse[i]=rmse(loocvTNW[[i]]$heldOut,loocvTNW[[i]]$pred)
}

# 5-fold k-fold cross-validation
# split data into 5 groups - I'm using the same splits as I used above.

overallKFCV.NW=overallKFCV.N

#looping to do k fold cross validation 10 times to make sure the random group assignment isn't having a large effect on rmse
overallRMSE.NW=data.frame(mod=character(),rmse=numeric())
for(k in 1:length(overallKFCV.NW)){
  kFCV.NW=list()
  for(j in 1:length(mForms)){
    tempDF=data.frame(groupLabel=overallKFCV.NW[[k]][["kfoldGroup"]],heldOut=q3$Mean,pred=numeric(nrow(q3)))
    for(i in 1:5){
      tempQ=overallKFCV.NW[[k]][overallKFCV.NW[[k]]$kfoldGroup!=i,]
      tempFit=glm(mForms[j],data = tempQ)
      tempDF$pred[tempDF$groupLabel==i]=predict(tempFit,overallKFCV.NW[[k]][overallKFCV.NW[[k]]$kfoldGroup==i,])
    }
    kFCV.NW[[j]]=list(mForms[j],tempDF)
  }
  rmseKFCV.NW=data.frame(mod=names(topNw),rmse=numeric(length = length(mForms)))
  for(i in 1:length(mForms)){
    rmseKFCV.NW$rmse[i]=rmse(kFCV.NW[[i]][[2]][["heldOut"]],kFCV.NW[[i]][[2]][["pred"]])
  }
  overallRMSE.NW=rbind(overallRMSE.NW,rmseKFCV.NW)
}

avgRMSE.NW=overallRMSE.NW%>%
  group_by(mod)%>%
  summarise(meanRMSE=mean(rmse))


#### set of models that make biological sense to compare to GA models ####
m1=(Mean~logBassCPE + gdd5 + logBassCPE:gdd5) # essentially hansen et al. formula
m2=(Mean~logBassCPE + logBassCPE:gdd5) #variation of hansen et al. 
m3=(Mean~logBassCPE + rip2 + logBassCPE:rip2) # bass + development
m4=(Mean~rip1*rip2*wshed1*wshed2) # land use effects only
m5=(Mean~fcomm1*fcomm2) # fish community effect only, presence/absence
mForms=list(m1,m2,m3,m4,m5)
specMods=data.frame(model=c("m1","m2","m3","m4","m5"), rmse=numeric(length = 5),loocvRMSE=numeric(length = 5))

for(i in 1:length(mForms)){
  tempGLM=glm(mForms[[i]],data = q4)
  specMods$rmse[i]=rmse(tempGLM$y,tempGLM$fitted.values)
  
  cvTemp=data.frame(heldOut=numeric(nrow(q4)),pred=numeric(nrow(q4)))
  for(j in 1:nrow(q4)){
    tempQ=q4[-j,]
    tempFit=glm(mForms[[i]],data=tempQ)
    cvTemp$pred[j]=predict(tempFit,q4[j,])
    cvTemp$heldOut[j]=q4$Mean[j]
  }
  specMods$loocvRMSE[i]=rmse(cvTemp$heldOut,cvTemp$pred)
  
}



#### PREDICTING OUT OF SET ####
### MOVE THIS TO NEW CLEANED UP SCRIPT ONCE FINAL MODEL HAS BEEN CHOSEN ###

#reading in lake characteristic data to make predictions
fishComm=read.csv("Hansen_WI_fishCommunity.csv",stringsAsFactors = F)
riparian=read.csv("Hansen_WI_riparian.csv",stringsAsFactors = F)
watershed=read.csv("Hansen_WI_watershed.csv",stringsAsFactors = F)
ralphGDD=read.csv("C:/Users/jones/BoxSync/NDstuff/Dissertation/4/Conversion_file_MedResTemp_to_WBIC/Conversion_file_MedResTemp_to_WBIC/output/NLDAS_mean_temperatures_WBIC.csv",stringsAsFactors = F)
d1=read.csv('C:/Users/jones/BoxSync/NDstuff/CNH/BASS_SEII_CPE.csv',stringsAsFactors = F)
wlyDens=read.csv("Hansen_WI_adultWLYdens.csv",stringsAsFactors = F)
stock=read.csv("stocking_receipts.csv",stringsAsFactors = F)
d4=read.csv('C:/Users/jones/BoxSync/NDstuff/CNH/TREATY_WALLEYE_PE.csv',stringsAsFactors = F)

#combine all the data and subset to the columns in our best fitting model
#I need lakes with walleye that aren't included in the original 82
keepAIC=c("wbic", "wshed2","fcomm2","rip1","overallMeanWLYDens","logBassCPE","stoDens","meanNSto","Area","meanLBsto","sumNSto","sumLBSto","stoLBDens")

#getting overallMeanWLYDens per acre for the prediction lakes

allWLY=d4%>%
  full_join(wlyDens, by=c("WBIC"="wbic"))%>%
  group_by(WBIC)%>%
  summarize(treatyWLY=mean(numberPerAcre,na.rm=T),
            hansenWLY=mean(adults.acre,na.rm=T))
allWLY$overallMeanWLYDens=rowMeans(allWLY[,2:3],na.rm=T)

sto=stock%>%
  filter(speciesCode=="X22")%>%
  group_by(WBIC)%>%
  summarise(meanNSto=mean(fishStockedNumber,na.rm=T),
            stoLen=mean(averageLength,na.rm=T),
            meanLBsto=mean(fishStockedLbs,na.rm=T),
            sumNSto=sum(fishStockedNumber,na.rm=T),
            sumLBSto=sum(fishStockedLbs,na.rm=T))
pred1=fishComm%>%
  full_join(riparian,by="wbic")%>%
  full_join(watershed,by="wbic")%>%
  full_join(ralphGDD,by=c("wbic"="WBIC"))%>%
  rename(ripDeveloped=Developed.x,ripBarren=Barren.x,ripForest=Forest.x,ripShrubScrub=Shrub.Scrub.x,
         ripGrassland=Grassland.x,ripPasture=Pasture.hay.x,ripCultivated=CultivatedCrops,
         ripWetlands=Wetlands.x,
         wshedDeveloped=Developed.y,wshedBarren=Barren.y,wshedForest=Forest.y,
         wshedShrubScrub=Shrub.Scrub.y,wshedGrassland=Grassland.y,wshedPasture=Pasture.hay.y,
         wshedCultivated=Cultivated,wshedWetlands=Wetlands.y)%>%
  filter(Walleye==1,!(wbic%in%dat$WBIC))
pred2=d1%>%
  group_by(WBIC)%>%
  summarise(bassCPE=mean(CPEmile,na.rm=T))%>%
  full_join(allWLY,by=c("WBIC"))%>%
  rename(wbic=WBIC)%>%
  filter(!(wbic%in%dat$WBIC))
pred3=pred2%>%
  full_join(pred1)%>%
  full_join(sto, by=c("wbic"="WBIC"))%>%
  mutate(stoDens=meanNSto/(Area*4047),
         stoLBDens=(meanLBsto*454)/(Area*4047)) #converting to g/m^2 and N/m^2

## RIPARIAN PCA
ripP=pred3[,c(1,17:24)]
ripP=ripP[!is.na(ripP$ripDeveloped),]
ripP.dist=dist(ripP[,-1])
ripP.pca=cmdscale(ripP.dist)
ripP.pca=as.data.frame(cbind(ripP$wbic,ripP.pca))
colnames(ripP.pca)=c("wbic","rip1","rip2")
#ripP.ef=envfit(ripP.pca,ripP)

pred4=pred3%>%
  inner_join(ripP.pca,by="wbic")%>%
  filter(!is.na(bassCPE))%>%
  mutate(logBassCPE=log(bassCPE+0.01))

predAIC=pred4[,colnames(pred4)%in%keepAIC]
