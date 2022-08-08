## CJD 5.17.21
## taking the resulting model from the selection process in 'glmultiModelSelection.R' script and making new out of set predictions
## also creating new figures for the paper

setwd("C:/Users/dassocju/OneDrive - State of Wisconsin/NotreDame/NDstuff/Dissertation/5")
rm(list=ls())
set.seed(1)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)

#### DATA PRE-PROCESSING ####

# see Hansen_WI_Co.variates.xlsx for metadata on covariates
# see Sassetal._Fisheries_SupplementalTables_Final.xlsx for meta data on depensation parms
dat=read.csv("qbAndCovariates.csv",stringsAsFactors = F)
fishComm=read.csv("Hansen_WI_fishCommunity.csv",stringsAsFactors = F)
riparian=read.csv("Hansen_WI_riparian.csv",stringsAsFactors = F)
watershed=read.csv("Hansen_WI_watershed.csv",stringsAsFactors = F)
ralphGDD=read.csv("C:/Users/dassocju/OneDrive - State of Wisconsin/NotreDame/NDstuff/Dissertation/4/Conversion_file_MedResTemp_to_WBIC/Conversion_file_MedResTemp_to_WBIC/output/NLDAS_mean_temperatures_WBIC.csv",stringsAsFactors = F)
d1=read.csv('C:/Users/dassocju/OneDrive - State of Wisconsin/NotreDame/NDstuff/CNH/BASS_SEII_CPE.csv',stringsAsFactors = F)

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
#% variance explained =
wshed.pca$eig[1]/sum(wshed.pca$eig) + wshed.pca$eig[2]/sum(wshed.pca$eig)

# plot(rip.pca[,1],rip.pca[,2])
# plot(rip.ef,add=TRUE)
# cor(cbind(rip.pca,rip))[1:2,]
# rip axis 1 = -Forest + Wetlands
# rip axis 2 = -Developed
#% variance explained =
rip.pca$eig[1]/sum(rip.pca$eig) + rip.pca$eig[2]/sum(rip.pca$eig)


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
# cor(cbind(comm.pcoa,comm))[1:2,]
# fcomm axis 1 = + panfish + lmb
# fcomm axis 2 = + musky + smb
#% variance explained =
comm.pcoa$eig[1]/sum(comm.pcoa$eig) + comm.pcoa$eig[2]/sum(comm.pcoa$eig)


comm.pcoa=as.data.frame(cbind(comm$wbic,comm.pcoa));colnames(comm.pcoa)=c("wbic","fcomm1","fcomm2")
rip.pca=as.data.frame(cbind(rip$wbic,rip.pca));colnames(rip.pca)=c("wbic","rip1","rip2")
wshed.pca=as.data.frame(cbind(wshed$wbic,wshed.pca));colnames(wshed.pca)=c("wbic","wshed1","wshed2")

# plots of PCAs for supplement
ggplot(data = wshed.pca,aes(x=wshed1,y=wshed2))+theme_classic()+
  geom_point()+
  labs(x="Watershed Land Use PC1",y="Watershed Land Use PC2")

ggplot(data=rip.pca,aes(x=rip1, y=rip2))+theme_classic()+
  geom_point()+
  labs(x="Riparian Land Use PC1", y="Riparian Land Use PC2")

ggplot(data=comm.pcoa,aes(x=fcomm1, y=fcomm2))+theme_classic()+
  geom_point()+
  labs(x="Fish Community PCo1",y="Fish Community PCo2")

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
vars=c("WBIC","Mean","fcomm1","wshed1","wshed2","rip1","rip2","gdd5","overallMeanWLYDens","logBassCPE","logN") # removing fcomm2 because it doesn't contribute much to explaining variation in the fish community data
bootStrapvars=c("WBIC","Mean","fcomm1","wshed1","wshed2","rip1","rip2","gdd5","overallMeanWLYDens", "overallSDWLYDens","logBassCPE", "logBassSDCPE","logN")# removing fcomm2 because it doesn't contribute much to explaining variation in the fish community data


q2=q[q$bassMCPE>0,colnames(q)%in%vars]
q3=q2[complete.cases(q2),-2]


#### BEST MODEL ####

# pulling best model structure from the averaged RMSE for 5 fold cross-validation of the -logN & -meanWLYDens set of predictors. This model is also the same one that you get from the -LogN only GA 5 fold cross validation. Same RMSE too

bestMod=glm(Mean ~ 1 + rip1:fcomm1 + logBassCPE:gdd5 + fcomm1:logBassCPE + rip1:logBassCPE + logBassCPE:wshed2, data=q3)

# model fit plot
ggplot()+theme_classic()+
  geom_point(aes(x=bestMod$y, y=bestMod$fitted.values))+
  labs(x="Sass et al. (2021)",y="Model Predicted")+
  geom_abline(slope = 1, intercept = 0)+
  coord_cartesian(ylim=c(0,1.5),xlim=c(0,1.5))+
  geom_vline(xintercept = 1,linetype=2)+
  geom_hline(yintercept = 1,linetype=2)

#### OUT OF SAMPLE PREDICTIONS ####
#reading in lake characteristic data to make predictions
wlyDens=read.csv("Hansen_WI_adultWLYdens.csv",stringsAsFactors = F)
stock=read.csv("stocking_receipts.csv",stringsAsFactors = F)
d4=read.csv('C:/Users/dassocju/OneDrive - State of Wisconsin/NotreDame/NDstuff/CNH/TREATY_WALLEYE_PE.csv',stringsAsFactors = F)

#combine all the data and subset to the columns in our best fitting model
#I need lakes with walleye that aren't included in the original 82
keepAIC=c("wbic", "wshed2","fcomm1","rip1","overallMeanWLYDens","logBassCPE","stoDens","meanNSto","Area","meanLBsto","sumNSto","sumLBSto","stoLBDens","gdd5")

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

pred4=pred3%>%
  inner_join(rip.pca,by="wbic")%>%
  inner_join(wshed.pca,by="wbic")%>%
  inner_join(comm.pcoa,by="wbic")%>%
  filter(!is.na(bassCPE))%>%
  mutate(logBassCPE=log(bassCPE+0.01))

predAIC=pred4[,colnames(pred4)%in%keepAIC]
predAIC=predAIC[complete.cases(predAIC),]
## PREDICTING Q
qAIC=predict.glm(bestMod,newdata = predAIC,se.fit = T,type = "response")
# comparison of predicted qs to the measured qs
mq=data.frame(obs=rep("Inference",nrow(q)),q=q$Mean,wbic=q$WBIC,se=q$SE,LCL=q$X2.50.,UCL=q$X97.50.)
pq1=data.frame(obs=rep("Predicted",length(qAIC$fit)),q=qAIC$fit,wbic=predAIC$wbic,se=qAIC$se.fit,LCL=numeric(nrow(predAIC)),UCL=numeric(nrow(predAIC)))

#### BOOTSTRAPPING CONFIDENCE INTERVALS ####
q2=q2[complete.cases(q2),] #getting just WBICS for lakes in q3 (these are the 28 inference lakes)
qCIdat=matrix(ncol=nrow(predAIC),nrow=1000)
for( i in 1:1000){
  #generate new observations of q
  q_new=rnorm(nrow(q2),mean = q$Mean[q$WBIC%in%q2$WBIC] ,sd = q$SD[q$WBIC%in%q2$WBIC])
  qFit=as.data.frame(cbind(q_new,q2))
  tempFit=glm(q_new ~1 + rip1:fcomm1 + logBassCPE:gdd5 + logBassCPE:fcomm1 + logBassCPE:rip1 + logBassCPE:wshed2,data=qFit)
  qPred=predict.glm(tempFit,newdata = predAIC,type = "response")
  qCIdat[i,]=qPred
}


for(i in 1:nrow(pq1)){
  pq1$LCL[i]=sort(qCIdat[,i])[25]
  pq1$UCL[i]=sort(qCIdat[,i])[975]
}

pQ=rbind(mq,pq1)
pQ$wbic=factor(pQ$wbic, levels = unique(pQ$wbic[order(pQ$q)]))
pQ=pQ[!is.na(pQ$q) & pQ$q>=0,]

p=ggplot(pQ,aes(x=q,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  geom_vline(xintercept = 1)+
  annotate(geom="text",x=0.35,y=2,label="Depensation")+
  annotate(geom="text",x=1.5,y=2,label="Compensation")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x=expression(italic("q")), y="density")
p
# there are 2 lakes with negative predictions, just throw them out?
   # I think so, they have extreme covariate values, doesn't suggest the model is bad, just that we have a few odd ducks and it doesn't change the takeaways from the results either.

# are the 2 distributions significantly different? 
ks.test(pQ$q[pQ$obs=="Inference"],pQ$q[pQ$obs=="Predicted"]) # not significantly different
pBass=data.frame(wbic=c(q2$WBIC,predAIC$wbic),
                 logBassCPE=c(q2$logBassCPE,predAIC$logBassCPE),
                 gdd5=c(q2$gdd5,predAIC$gdd5),
                 rip1=c(q2$rip1,predAIC$rip1),
                 wshed2=c(q2$wshed2,predAIC$wshed2),
                 fcomm1=c(q2$fcomm1,predAIC$fcomm1),
                 obs=c(rep("Inference",nrow(q2)),rep("Predicted",nrow(predAIC))))
pb=ggplot(pBass,aes(x=logBassCPE,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x="Log of Largemouth Bass Catch-per-Km")
pb
pgdd=ggplot(pBass,aes(x=gdd5,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x=expression("Growing Degree Days Above 5" *~degree*C))

prip=ggplot(pBass,aes(x=rip1,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x="Riparian Land Use PC1")

pwshed=ggplot(pBass,aes(x=wshed2,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x="Watershed Land Use PC2")

pfcomm=ggplot(pBass,aes(x=fcomm1,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x="Fish Community PCO1")

pb;t.test(pBass$logBassCPE[pBass$obs=="Inference"],pBass$logBassCPE[pBass$obs=="Predicted"]) # sig diff, data not normally distributed, do no use, use mann-whitney below
pgdd;t.test(pBass$gdd5[pBass$obs=="Inference"],pBass$gdd5[pBass$obs=="Predicted"]) # predicted observations not normal, do not use, use mann-whitney below
prip;t.test(pBass$rip1[pBass$obs=="Inference"],pBass$rip1[pBass$obs=="Predicted"]) # sig diff, predicted observations not normal, do not use, use mann-whitney below
pwshed;t.test(pBass$wshed2[pBass$obs=="Inference"],pBass$wshed2[pBass$obs=="Predicted"]) # sig diff, predicted observations not normal, do not use, use mann-whitney below
pfcomm;t.test(pBass$fcomm1[pBass$obs=="Inference"],pBass$fcomm1[pBass$obs=="Predicted"]) #  data not normally distributed, do not use, use mann-whitney below


# doing the same with mann-whitney U test to see if distributions are different since they're not normally distributed (according to Shapiro-wilk test), did same with Kolmogorov-Smirnov test since it better tests differences in distributions as opposed to differences in rank like Mann-Whitney. However, Mann-Whitney can better deal with ties in the data for non-continuous data which our pca data sort of is since it comes from non-continuous measures.

# bottom line is they generally agree except for the comparison between fcomm1 in inference and predicted sets. 

# use Mann-Whitney for rip1, wshed2, and fcomm1 comparisons and K-S for logBassCPE and gdd5 comparisons. 

pb;wilcox.test(pBass$logBassCPE[pBass$obs=="Inference"],pBass$logBassCPE[pBass$obs=="Predicted"], alternative = "two.sided") # sig diff
pgdd;wilcox.test(pBass$gdd5[pBass$obs=="Inference"],pBass$gdd5[pBass$obs=="Predicted"], alternative = "two.sided")
prip;wilcox.test(pBass$rip1[pBass$obs=="Inference"],pBass$rip1[pBass$obs=="Predicted"], alternative = "two.sided") # sig diff
pwshed;wilcox.test(pBass$wshed2[pBass$obs=="Inference"],pBass$wshed2[pBass$obs=="Predicted"], alternative = "two.sided")
pfcomm;wilcox.test(pBass$fcomm1[pBass$obs=="Inference"],pBass$fcomm1[pBass$obs=="Predicted"], alternative = "two.sided") # sig diff

pb;ks.test(pBass$logBassCPE[pBass$obs=="Inference"],pBass$logBassCPE[pBass$obs=="Predicted"], alternative = "two.sided") # sig diff
pgdd;ks.test(pBass$gdd5[pBass$obs=="Inference"],pBass$gdd5[pBass$obs=="Predicted"], alternative = "two.sided")
prip;ks.test(pBass$rip1[pBass$obs=="Inference"],pBass$rip1[pBass$obs=="Predicted"], alternative = "two.sided") # sig diff
pwshed;ks.test(pBass$wshed2[pBass$obs=="Inference"],pBass$wshed2[pBass$obs=="Predicted"], alternative = "two.sided")
pfcomm;ks.test(pBass$fcomm1[pBass$obs=="Inference"],pBass$fcomm1[pBass$obs=="Predicted"], alternative = "two.sided") 

ggarrange(p,pb,prip,pfcomm,pgdd,pwshed, ncol=2, nrow=3, legend = "top", common.legend = T, labels = "auto")

paic=ggplot(pQ,aes(x=q,y=wbic))+theme_classic()+
  geom_pointrange(aes(xmin=LCL,xmax=UCL))+
  facet_wrap(~obs,scales="free")+
  geom_vline(xintercept = 1)+
  theme(axis.text.y = element_blank())+
  xlab(expression(italic("q")))+
  ylab("")
paic
nrow(pQ[pQ$obs=="Predicted" & pQ$LCL<1 & pQ$UCL>1,]) # number of lakes who's CIs overlap 1
nrow(pQ[pQ$obs=="Predicted" & pQ$LCL>1 & pQ$UCL>1,]) # number of lakes fully above 1
nrow(pQ[pQ$obs=="Predicted" & pQ$LCL<1 & pQ$UCL<1,]) # number of lakes fully below 1

# same plot as above but horizonally oriented for ppt presentations
paicH=ggplot(pQ[pQ$obs=="Predicted",],aes(x=q,y=wbic))+theme_classic()+
  geom_pointrange(aes(xmin=LCL,xmax=UCL))+
  geom_vline(xintercept = 1)+
  coord_flip()+theme(axis.text.x = element_blank())+ylab("")
paicH
# quadrat plots of q v density and q v stocking
pWLY=data.frame(wbic=c(q2$WBIC,predAIC$wbic),
                overallMeanWLYDens=c(q2$overallMeanWLYDens,predAIC$overallMeanWLYDens),
                obs=c(rep("Inference",nrow(q2)),rep("Predicted",nrow(predAIC))))
pWLY$wbic=as.factor(pWLY$wbic)
pD=pWLY%>%
  full_join(pQ, by=c("wbic","obs"))
qualityDep=sum(pD$overallMeanWLYDens>=3&pD$q<1,na.rm=T)
qualityCom=sum(pD$overallMeanWLYDens>=3&pD$q>=1,na.rm=T)
nonqualDep=sum(pD$overallMeanWLYDens<3&pD$q<1,na.rm=T)
nonqualCom=sum(pD$overallMeanWLYDens<3&pD$q>=1,na.rm=T)

# Kolmogorov-Smirnov tests on walleye density data to see if there are any differences, not really sure what the use of this would be, I'm not putting it in the paper for now
ks.test(pD$overallMeanWLYDens[pD$q<1],pD$overallMeanWLYDens[pD$q>=1])
ks.test(pD$q[pD$overallMeanWLYDens>=3],pD$q[pD$overallMeanWLYDens<3])

pd=ggplot()+theme_classic()+ # converting density in acres to ha
  geom_rect(aes(ymax=3/.405,ymin=-1.5,xmin=-1.5,xmax=1),alpha=0.1,fill="red")+
  geom_rect(aes(ymax=27,ymin=3/.405,xmin=-1.5,xmax=1),alpha=0.1,fill="yellow")+
  geom_rect(aes(ymax=3/.405,ymin=-1.5,xmin=1,xmax=2),alpha=0.1,fill="yellow")+
  geom_rect(aes(ymax=27,ymin=3/.405,xmin=1,xmax=2),alpha=0.1,fill="green")+
  geom_point(data=pD,aes(x=q,y=(overallMeanWLYDens/.405),color=obs))+
  theme(legend.position = "bottom",legend.title = element_blank())+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 7.4)+
  labs(x=expression(italic("q")),y="Adult walleye density \n no./ha")+
  geom_text(aes(x=0.35,y=20,label=qualityDep),color="black")+
  geom_text(aes(x=0.35,y=0,label=nonqualDep),color="black")+
  geom_text(aes(x=1.25,y=20,label=qualityCom),color="black")+
  geom_text(aes(x=1.25,y=0,label=nonqualCom),color="black")+
  scale_color_manual(values = c("black","grey50"))+
  coord_cartesian(xlim=c(0,2),ylim = c(0,25))
pd

q$stoDens=q$meanNSto/(q$Area*4047)# converting to N/m^2
q$stoLBDens=(q$meanLBsto*454)/(q$Area*4047) # converting to g/m^2
pSto=data.frame(wbic=c(q2$WBIC,predAIC$wbic),
                numStocked=c(q$meanNSto[q$WBIC%in%q2$WBIC],predAIC$meanNSto),
                denStocked=c(q$stoDens[q$WBIC%in%q2$WBIC],predAIC$stoDens),
                lbsStocked=c(q$meanLBsto[q$WBIC%in%q2$WBIC],predAIC$meanLBsto),
                denLBStocked=c(q$stoLBDens[q$WBIC%in%q2$WBIC],predAIC$stoLBDens),
                sumNStocked=c(q$sumNSto[q$WBIC%in%q2$WBIC],predAIC$sumNSto),
                sumLBStocked=c(q$sumLBSto[q$WBIC%in%q2$WBIC],predAIC$sumLBSto),
                obs=c(rep("Inference",nrow(q2)),rep("Predicted",nrow(predAIC))))
pSto$wbic=as.factor(pSto$wbic)
pSto=pSto%>%
  full_join(pQ,by=c("wbic","obs"))%>%
  filter(!is.na(q))
plb=ggplot(data=pSto,aes(x=q,y=denLBStocked,color=obs))+theme_classic()+
  geom_point()+
  theme(legend.position = "bottom",legend.title = element_blank())+
  geom_vline(xintercept = 1)
plb

nrow(pSto[pSto$q<1 & !is.na(pSto$numStocked),]) # number of depensation lakes with stocking records
nrow(pSto[pSto$q<1 & is.na(pSto$numStocked),])
nrow(pSto[pSto$q>=1 & !is.na(pSto$numStocked),]) # number of compensatory lakes with stocking records
nrow(pSto[pSto$q>=1 & is.na(pSto$numStocked),])

plb2=ggplot(data=pSto,aes(x=q,y=log(denLBStocked),color=obs))+theme_classic()+
  geom_point()+
  theme(legend.position = "bottom",legend.title = element_blank())+
  geom_vline(xintercept = 1)+
  labs(x=expression(italic("q")),y = expression(paste("Log mean g/",m^{2},' Stocked')))+
  geom_text(x=0.4,y=-8, label=nrow(pSto[pSto$q<1 & !is.na(pSto$denLBStocked),]),color="black")+
  geom_text(x=1.25,y=-8, label=nrow(pSto[pSto$q>=1 & !is.na(pSto$denLBStocked),]),color="black")+
  scale_color_manual(values = c("black","grey"))
plb2

ps=ggplot(pSto,aes(x=q,color=!is.na(pSto$numStocked)))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  scale_color_manual(name=" ",labels=c("Not stocked","Stocked"),values = c("black","grey"))+
  labs(x=expression(italic("q")))
ps
pr=ggplot(pSto,aes(x=log(denLBStocked),color=q>=1))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  labs(x=expression(paste("Log mean g/",m^{2},' Stocked')))+
  scale_color_manual(name=" ",labels=c("Depensation","Compensation"),values = c("black","grey"))
pr
t.test(pSto$q[!is.na(pSto$numStocked)],pSto$q[is.na(pSto$numStocked)]) # no sig diff between q in stocked and unstocked lakes. Not normally distributed, use mann-whitney
t.test(log(pSto$denLBStocked[pSto$q>=1]),log(pSto$denLBStocked[pSto$q<1])) # no sig diff stocking density in compensatory and depensatory lakes. Not normally distributed, use mann-whitney

wilcox.test(pSto$q[!is.na(pSto$numStocked)],pSto$q[is.na(pSto$numStocked)], alternative = "two.sided") # no significant difference in q for stocked and not stocked lakes
ks.test(log(pSto$denLBStocked[pSto$q>=1]),log(pSto$denLBStocked[pSto$q<1]))# no significant difference in stocking density for compensatory and depensatory lakes

#top=ggpubr::ggarrange(plb2,labels = "auto", common.legend = T)
bot=ggpubr::ggarrange(ps,pr,labels = "auto",legend = "bottom",ncol = 1) # just put this panel in the paper, the other plot that I had grouped with this is redundant
#ggpubr::ggarrange(top,bot,ncol=1)


# how many stocking events have fish >200 mm TL = 465
sum(stock$averageLength[stock$speciesCode=="X22"]*25.4>200,na.rm = T)
# how many <200mm TL = 5120
sum(stock$averageLength[stock$speciesCode=="X22"]*25.4<200,na.rm = T)
# pWLY=data.frame(wbic=c(qstep$WBIC,predAIC$wbic),
#                  avgWLY=c(qstep$avgWLY,predAIC$avgWLY),
#                  obs=c(rep("measured",nrow(qstep)),rep("predicted",nrow(predAIC))))
# pw=ggplot(pWLY,aes(x=avgWLY,color=obs))+theme_classic()+
#   geom_density()+
#   theme(legend.position = "bottom",legend.title = element_blank())
# pw
# pRip2=data.frame(wbic=c(qstep$WBIC,predAIC$wbic),
#                  rip2=c(qstep$rip2,predAIC$rip2),
#                  obs=c(rep("measured",nrow(qstep)),rep("predicted",nrow(predAIC))))
# pr=ggplot(pRip2,aes(x=rip2,color=obs))+theme_classic()+
#   geom_density()+
#   theme(legend.position = "bottom",legend.title = element_blank())
# pr

#### plots to understand bestMod coefs ####

summary(bestMod)

# rip1:fcomm1
ggplot(data = q3, aes(y=Mean,x=rip1,color=fcomm1))+theme_classic()+
  geom_point(size=2)+scale_color_binned(breaks=seq(-0.6,0.3, by=0.1))

#logBassCPE:gdd5
ggplot(data = q3, aes(y=Mean,x=logBassCPE,size=gdd5))+theme_classic()+
  geom_point()+coord_cartesian(xlim=c(-2,6))

#fcomm1:logBassCPE
ggplot(data = q3, aes(y=Mean,x=logBassCPE,color=fcomm1))+theme_classic()+
  geom_point()+coord_cartesian(xlim=c(-2,6))+
  scale_color_binned(breaks = seq(-0.6, 0.3, by=0.1))

ggplot(data = q3, aes(y=fcomm1,x=logBassCPE,size=Mean))+theme_classic()+
  geom_point()+coord_cartesian(xlim = c(-2,6))

#rip1:logBassCPE
ggplot(data = q3, aes(y=rip1,x=logBassCPE,size=Mean))+theme_classic()+
  geom_point()+coord_cartesian(xlim=c(-2,6))
ggplot(data=q3, aes(y=rip1, x=logBassCPE))+theme_classic()+
  geom_point()+coord_cartesian(xlim=c(-2,6))

#logBassCPE:wshed2
ggplot(data = q3, aes(y=Mean,x=logBassCPE,size=wshed2))+theme_classic()+
  geom_point()+coord_cartesian(xlim=c(-2,6))
ggplot(data=q3, aes(y=wshed2, x=logBassCPE))+theme_classic()+
  geom_point()+coord_cartesian(xlim=c(-2,6))

#### EG Fingerling Stocking Stats in Discussion ####

#pull in stocking receipts
stoDat=read.csv("C:/Users/dassocju/Documents/replacementCosts/Walleye_stocking_receipts_2002_to_present.csv",stringsAsFactors=F)
stoDat$Average.Length=as.numeric(stoDat$Average.Length)
sum(stoDat$Average.Length>7.8,na.rm=T)

# plot of the distribution of q values from Sass et al. to be used in my WIAFS 2022 talk

ggplot(dat[dat$parm=="q",],aes(x=Mean))+theme_classic()+
  geom_density(size=1)+labs(x=expression(italic("q")),y=element_blank())+
  geom_vline(xintercept = 1)+theme(text = element_text(size=25))

# % of original qs overlapping or above/below 1

allq=dat[dat$parm=="q",]
nrow(allq[allq$X2.50.<1 & allq$X97.50.<1,]) #fully below 1 = 18
nrow(allq[allq$X2.50.<1 & allq$X97.50.>1,]) #overlapping 1 = 63
nrow(allq[allq$X2.50.>1 & allq$X97.50.>1,]) #fully above 1 = 1
