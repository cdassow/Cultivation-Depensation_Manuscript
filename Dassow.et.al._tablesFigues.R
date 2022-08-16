## CJD  Originate 5.17.21
## taking the resulting model from the selection process in 'Dassow.et.al._modelSelection.R' script and making new out of set predictions for management application
## also creating figures for the paper

setwd("~/OAS_git_repos/Cultivation-Depensation_Manuscript")
rm(list=ls())
set.seed(1)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)

#### DATA PRE-PROCESSING ####

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

# plots of PCAs for supplement
ggplot(data = wshed.pca,aes(x=wshed1,y=wshed2))+theme_classic()+
  geom_point()+
  labs(x="Watershed Land Use PC1",y="Watershed Land Use PC2")

ggplot(data=rip.pca,aes(x=rip1, y=rip2))+theme_classic()+
  geom_point()+
  labs(x="Riparian Land Use PC1", y="Riparian Land Use PC2")

#adding the useful axes to the df
dat=dat%>%
  left_join(comm.pcoa, by=c("WBIC"="wbic"))%>%
  left_join(rip.pca, by=c("WBIC"="wbic"))%>%
  left_join(wshed.pca, by=c("WBIC"="wbic"))

#### SETTING UP DATA FRAME FOR INFERENCE LAKE SET ####
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



#### BEST MODEL ####

# pulling best model structure from the averaged RMSE for 5 fold cross-validation.

bestMod=glm(Mean ~ 1 + rip1:fcomm1 + logBassCPE:gdd5 + fcomm1:logBassCPE + rip1:logBassCPE + logBassCPE:wshed2, data=q3)

# model fit plot - FIGURE 2
ggplot()+theme_classic()+
  geom_point(aes(x=bestMod$y, y=bestMod$fitted.values),size=3)+
  labs(x="Sass et al. (2021)",y="Model Predicted")+
  geom_abline(slope = 1, intercept = 0)+
  coord_cartesian(ylim=c(0,1.5),xlim=c(0,1.5))+
  geom_vline(xintercept = 1,linetype=2)+
  geom_hline(yintercept = 1,linetype=2)+
  theme(text = element_text(size=20,family = "serif"))

#### OUT OF SAMPLE PREDICTIONS ####
#reading in lake characteristic data to make predictions and apply results to management
wlyDens=read.csv("Hansen_WI_adultWLYdens.csv",stringsAsFactors = F)
stock=read.csv("stocking_receipts.csv",stringsAsFactors = F)
d4=read.csv('TREATY_WALLEYE_PE.csv',stringsAsFactors = F)

#combine all the data and subset to the columns in our best fitting model
#I need lakes with walleye that aren't included in the original 82 inference set of lakes from Sass et al. (2021), remember that we only used 28/82 lake given our data requirements
keepAIC=c("wbic", "wshed2","fcomm1","rip1","overallMeanWLYDens","logBassCPE","stoDens","meanNSto","Area","meanLBsto","sumNSto","sumLBSto","stoLBDens","gdd5")

#getting overallMeanWLYDens per acre for the prediction lakes from a couple different data sources that differ slightly. This allows for the maximum number of lakes to be included for prediction. In cases were walleye densities were availble for a given lake from both data sets then the average of the two is taken.

allWLY=d4%>% # adult walleye densities summarized by lake
  full_join(wlyDens, by=c("WBIC"="wbic"))%>%
  group_by(WBIC)%>%
  summarize(treatyWLY=mean(numberPerAcre,na.rm=T),
            hansenWLY=mean(adults.acre,na.rm=T))
allWLY$overallMeanWLYDens=rowMeans(allWLY[,2:3],na.rm=T)

sto=stock%>% # walleye stocking densities and lengths summarized by lake
  filter(speciesCode=="X22")%>%
  group_by(WBIC)%>%
  summarise(meanNSto=mean(fishStockedNumber,na.rm=T),
            stoLen=mean(averageLength,na.rm=T),
            meanLBsto=mean(fishStockedLbs,na.rm=T),
            sumNSto=sum(fishStockedNumber,na.rm=T),
            sumLBSto=sum(fishStockedLbs,na.rm=T))
pred1=fishComm%>% # model predictor values for the prediction lakes
  full_join(riparian,by="wbic")%>%
  full_join(watershed,by="wbic")%>%
  full_join(ralphGDD,by=c("wbic"="WBIC"))%>%
  filter(Walleye==1,!(wbic%in%dat$WBIC))
pred2=d1%>% # another set of model predictor values for the prediction lakes, because of the data structure these are handled separately and combined later
  group_by(WBIC)%>%
  summarise(bassCPE=mean(CPEmile,na.rm=T))%>%
  full_join(allWLY,by=c("WBIC"))%>%
  rename(wbic=WBIC)%>%
  filter(!(wbic%in%dat$WBIC))
pred3=pred2%>% # combining all the necessary predictor values in order to predict q values for the set of prediction lakes identified
  full_join(pred1)%>%
  full_join(sto, by=c("wbic"="WBIC"))%>%
  mutate(stoDens=meanNSto/(Area*4047),
         stoLBDens=(meanLBsto*454)/(Area*4047)) #converting to g/m^2 and N/m^2

pred4=pred3%>% # final bit of code to combine all the necessary predictor values in order to predict q values for the set of prediction lakes identified
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

#### MONTE CARLO CONFIDENCE INTERVALS ####
q2=q2[complete.cases(q2),] #getting just WBICS for lakes in q3 (these are the 28 inference lakes)
qCIdat=matrix(ncol=nrow(predAIC),nrow=1000)
for( i in 1:1000){
  #generate new observations of q
  q_new=rnorm(nrow(q2),mean = q$Mean[q$WBIC%in%q2$WBIC] ,sd = q$SD[q$WBIC%in%q2$WBIC])
  qFit=as.data.frame(cbind(q_new,q2))
  #fit model to 'new' 28 inference lake dataset
  tempFit=glm(q_new ~1 + rip1:fcomm1 + logBassCPE:gdd5 + logBassCPE:fcomm1 + logBassCPE:rip1 + logBassCPE:wshed2,data=qFit)
  #predict out-of-set lakes (these are the 115 predition lakes)
  qPred=predict.glm(tempFit,newdata = predAIC,type = "response")
  qCIdat[i,]=qPred
}

# getting confidence intervals
for(i in 1:nrow(pq1)){
  pq1$LCL[i]=sort(qCIdat[,i])[25]
  pq1$UCL[i]=sort(qCIdat[,i])[975]
}

pQ=rbind(mq,pq1)
pQ$wbic=factor(pQ$wbic, levels = unique(pQ$wbic[order(pQ$q)]))
pQ=pQ[!is.na(pQ$q) & pQ$q>=0,]

#comparing distributions of q vales for inference and predicted lake sets - APPENDIX S1 FIG. S4a
p=ggplot(pQ,aes(x=q,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  geom_vline(xintercept = 1)+
  annotate(geom="text",x=0.35,y=2,label="Depensation")+
  annotate(geom="text",x=1.5,y=2,label="Compensation")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x=expression(italic("q")), y="density")
p

# are the 2 distributions significantly different? 
ks.test(pQ$q[pQ$obs=="Inference"],pQ$q[pQ$obs=="Predicted"]) # not significantly different at an alpha of 0.05
pBass=data.frame(wbic=c(q2$WBIC,predAIC$wbic),
                 logBassCPE=c(q2$logBassCPE,predAIC$logBassCPE),
                 gdd5=c(q2$gdd5,predAIC$gdd5),
                 rip1=c(q2$rip1,predAIC$rip1),
                 wshed2=c(q2$wshed2,predAIC$wshed2),
                 fcomm1=c(q2$fcomm1,predAIC$fcomm1),
                 obs=c(rep("Inference",nrow(q2)),rep("Predicted",nrow(predAIC))))
# APPENDIX S1 FIG. S4b
pb=ggplot(pBass,aes(x=logBassCPE,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x="Log of Largemouth Bass Catch-per-Km")
pb

# APPENDIX S1 FIG. S4e
pgdd=ggplot(pBass,aes(x=gdd5,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x=expression("Growing Degree Days Above 5" *~degree*C))

# APPENDIX S1 FIG. S4c
prip=ggplot(pBass,aes(x=rip1,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x="Riparian Land Use PC1")

# APPENDIX S1 FIG. S4f
pwshed=ggplot(pBass,aes(x=wshed2,color=obs))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  scale_color_manual(values = c("black","grey"))+
  labs(x="Watershed Land Use PC2")

# APPENDIX S1 FIG. S4d
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


# doing the same with mann-whitney U test to see if distributions are different since they're not normally distributed (according to Shapiro-wilk test), did same with Kolmogorov-Smirnov test since it better tests differences in distributions as opposed to differences in rank like Mann-Whitney. However, Mann-Whitney can better deal with ties in the data for non-continuous data which our pca data sort is since it comes from non-continuous measures.

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

# APPENDIX S1 FIG. S4 - FULL PANEL FIGURE
ggarrange(p,pb,prip,pfcomm,pgdd,pwshed, ncol=2, nrow=3, legend = "top", common.legend = T, labels = "auto")

# FIGURE 3 - uncertainty in q estimates
paic=ggplot(pQ,aes(x=q,y=wbic))+theme_classic()+
  geom_pointrange(aes(xmin=LCL,xmax=UCL))+
  facet_wrap(~obs,scales="free")+
  geom_vline(xintercept = 1)+
  theme(axis.text.y = element_blank(), text = element_text(size=20, family = "serif"))+
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
# quadrat plots of q vs. density and q vs. stocking
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

# Kolmogorov-Smirnov tests on walleye density data  and whether or not a lake is stocked to see if there are any differences
ks.test(pD$overallMeanWLYDens[pD$q<1],pD$overallMeanWLYDens[pD$q>=1])
ks.test(pD$q[pD$overallMeanWLYDens>=3],pD$q[pD$overallMeanWLYDens<3])

# FIGURE 4 - quadrat plot of q value and average adult walleye density to prioritize lakes for different actions
pd=ggplot()+theme_classic()+ # converting density in acres to ha
  geom_rect(aes(ymax=3/.405,ymin=-1.5,xmin=-1.5,xmax=1),alpha=0.2,fill="#fc8d62")+
  geom_rect(aes(ymax=27,ymin=3/.405,xmin=-1.5,xmax=1),alpha=0.2,fill="#8da0cb")+
  geom_rect(aes(ymax=3/.405,ymin=-1.5,xmin=1,xmax=2),alpha=0.2,fill="#8da0cb")+
  geom_rect(aes(ymax=27,ymin=3/.405,xmin=1,xmax=2),alpha=0.2,fill="#66c2a5")+
  geom_point(data=pD,aes(x=q,y=(overallMeanWLYDens/.405),color=obs),size=2)+
  theme(legend.position = "bottom",legend.title = element_blank(), text = element_text(size=20, family = "serif"))+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 7.4)+
  labs(x=expression(italic("q")),y="Adult walleye density \n no./ha")+
  geom_text(aes(x=0.35,y=20,label=qualityDep),color="black", family="serif", size=8)+
  geom_text(aes(x=0.35,y=0,label=nonqualDep),color="black", family="serif", size=8)+
  geom_text(aes(x=1.25,y=20,label=qualityCom),color="black", family="serif", size=8)+
  geom_text(aes(x=1.25,y=0,label=nonqualCom),color="black", family="serif", size=8)+
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

nrow(pSto[pSto$q<1 & !is.na(pSto$numStocked),]) # number of depensation lakes with stocking records
nrow(pSto[pSto$q<1 & is.na(pSto$numStocked),])
nrow(pSto[pSto$q>=1 & !is.na(pSto$numStocked),]) # number of compensatory lakes with stocking records
nrow(pSto[pSto$q>=1 & is.na(pSto$numStocked),])

# a look at stocking densities compared to q values, not included in manuscript
plb2=ggplot(data=pSto,aes(x=q,y=log(denLBStocked),color=obs))+theme_classic()+
  geom_point()+
  theme(legend.position = "bottom",legend.title = element_blank())+
  geom_vline(xintercept = 1)+
  labs(x=expression(italic("q")),y = expression(paste("Log mean g/",m^{2},' Stocked')))+
  geom_text(x=0.4,y=-8, label=nrow(pSto[pSto$q<1 & !is.na(pSto$denLBStocked),]),color="black")+
  geom_text(x=1.25,y=-8, label=nrow(pSto[pSto$q>=1 & !is.na(pSto$denLBStocked),]),color="black")+
  scale_color_manual(values = c("black","grey"))
plb2

# FIGURE 5a
ps=ggplot(pSto,aes(x=q,color=!is.na(pSto$numStocked)))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  scale_color_manual(name=" ",labels=c("Not stocked","Stocked"),values = c("black","grey"))+
  labs(x=expression(italic("q")),y="Density")+
  theme(text = element_text(family = "serif", size = 20))
ps

#GIFURE 5b
pr=ggplot(pSto,aes(x=log(denLBStocked),color=q>=1))+theme_classic()+
  geom_density(size=1,key_glyph = "path")+
  labs(x=expression(paste("Log mean g/",m^{2},' Stocked')),y="Density")+
  scale_color_manual(name=" ",labels=c("Depensation","Compensation"),values = c("black","grey"))+
  theme(text = element_text(family = "serif", size = 20))
pr
t.test(pSto$q[!is.na(pSto$numStocked)],pSto$q[is.na(pSto$numStocked)]) # no sig diff between q in stocked and unstocked lakes. Not normally distributed, use mann-whitney
t.test(log(pSto$denLBStocked[pSto$q>=1]),log(pSto$denLBStocked[pSto$q<1])) # no sig diff stocking density in compensatory and depensatory lakes. Not normally distributed, use mann-whitney

wilcox.test(pSto$q[!is.na(pSto$numStocked)],pSto$q[is.na(pSto$numStocked)], alternative = "two.sided") # no significant difference in q for stocked and not stocked lakes
ks.test(log(pSto$denLBStocked[pSto$q>=1]),log(pSto$denLBStocked[pSto$q<1]))# no significant difference in stocking density for compensatory and depensatory lakes

# FIGURE 5 - full panel plot
ggpubr::ggarrange(ps,pr,labels = c("(a)","(b)"),legend = "bottom",ncol = 1) 


# how many stocking events have fish >200 mm TL = 465
sum(stock$averageLength[stock$speciesCode=="X22"]*25.4>200,na.rm = T)
# how many <200mm TL = 5120
sum(stock$averageLength[stock$speciesCode=="X22"]*25.4<200,na.rm = T)


#### plots to understand bestMod coefs ####
# this was used to help craft discussion about the effect of increasing bassCPE for example on q values when interacting with another predictor. Plots not included in manuscript
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
