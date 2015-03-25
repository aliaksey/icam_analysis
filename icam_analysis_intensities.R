rm(list=ls())
load("Icam_cells_after_correction.RDATA")
##collapse all data to image repeat
perobindall.cort<-perobindall.cor[perobindall.cor$Image_Metadata_array<9,]
library(plyr) 
icamintim<-ddply(perobindall.cort,"ImageNumber", summarise,  
                 FeatureIdx=unique(FeatureIdx),
                 Array=unique(Image_Metadata_array),
                 #Icam
                 ## noncorrected intensity 
                 #median
                 IcamIntensityMeanmed=mean(Cell_Intensity_MedianIntensity_Icam1amask,trim=0.2),
                 IcamIntensityMedmed=median(Cell_Intensity_MedianIntensity_Icam1amask),
                 ##mean
                 IcamIntensityMeanmean=mean(Cell_Intensity_MeanIntensity_Icam1amask,trim=0.2),
                 IcamIntensityMedmean=median(Cell_Intensity_MeanIntensity_Icam1amask),
                 ##integrated
                 IcamIntensityMeanint=mean(Cell_Intensity_IntegratedIntensity_Icam1amask,trim=0.2),
                 IcamIntensityMedint=median(Cell_Intensity_IntegratedIntensity_Icam1amask),
                 
                 ## corrected intensity min intensity subtracted
                 #median
                 IcamIntensityMeanmed_Nor=mean(Cell_Intensity_MedianIntensity_Icam1amask_Nor,trim=0.2),
                 IcamIntensityMedmed_Nor=median(Cell_Intensity_MedianIntensity_Icam1amask_Nor),
                 ##mean
                 IcamIntensityMeanmean_Nor=mean(Cell_Intensity_MeanIntensity_Icam1amask_Nor,trim=0.2),
                 IcamIntensityMedmean_Nor=median(Cell_Intensity_MeanIntensity_Icam1amask_Nor),
                 ##integrated
                 IcamIntensityMeanint_Nor=mean(Cell_Intensity_IntegratedIntensity_Icam1amask_Nor,trim=0.2),
                 IcamIntensityMedint_Nor=median(Cell_Intensity_IntegratedIntensity_Icam1amask_Nor),
                 
                 ## corrected intensity per batch correction
                 #median
                 IcamIntensityMeanmed_Nor_corr=mean(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr,trim=0.2),
                 IcamIntensityMedmed_Nor_corr=median(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr),
                 ##mean
                 IcamIntensityMeanmean_Nor_corr=mean(Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr,trim=0.2),
                 IcamIntensityMedmean_Nor_corr=median(Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr),
                 ##integrated
                 IcamIntensityMeanint_Nor_corr=mean(Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr,trim=0.2),
                 IcamIntensityMedint_Nor_corr=median(Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr),
                 
                 #Actin
                 ## noncorrected intensity 
                 #median
                 ActinIntensityMeanmed=mean(Cell_Intensity_MedianIntensity_Actin1amask,trim=0.2),
                 ActinIntensityMedmed=median(Cell_Intensity_MedianIntensity_Actin1amask),
                 ##mean
                 ActinIntensityMeanmean=mean(Cell_Intensity_MeanIntensity_Actin1amask,trim=0.2),
                 ActinIntensityMedmean=median(Cell_Intensity_MeanIntensity_Actin1amask),
                 ##integrated
                 ActinIntensityMeanint=mean(Cell_Intensity_IntegratedIntensity_Actin1amask,trim=0.2),
                 ActinIntensityMedint=median(Cell_Intensity_IntegratedIntensity_Actin1amask),
                 
                 ## corrected intensity min intensity subtracted
                 #median
                 ActinIntensityMeanmed_Nor=mean(Cell_Intensity_MedianIntensity_Actin1amask_Nor,trim=0.2),
                 ActinIntensityMedmed_Nor=median(Cell_Intensity_MedianIntensity_Actin1amask_Nor),
                 ##mean
                 ActinIntensityMeanmean_Nor=mean(Cell_Intensity_MeanIntensity_Actin1amask_Nor,trim=0.2),
                 ActinIntensityMedmean_Nor=median(Cell_Intensity_MeanIntensity_Actin1amask_Nor),
                 ##integrated
                 ActinIntensityMeanint_Nor=mean(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor,trim=0.2),
                 ActinIntensityMedint_Nor=median(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor),
                 
                 ## corrected intensity per batch correction
                 #median
                 ActinIntensityMeanmed_Nor_corr=mean(Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr,trim=0.2),
                 ActinIntensityMedmed_Nor_corr=median(Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr),
                 ##mean
                 ActinIntensityMeanmean_Nor_corr=mean(Cell_Intensity_MeanIntensity_Actin1amask_Nor_corr,trim=0.2),
                 ActinIntensityMedmean_Nor_corr=median(Cell_Intensity_MeanIntensity_Actin1amask_Nor_corr),
                 ##integrated
                 ActinIntensityMeanint_Nor_corr=mean(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor_corr,trim=0.2),
                 ActinIntensityMedint_Nor_corr=median(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor_corr))

##calculate the same for negative control

perobindall.corn<-perobindall.cor[perobindall.cor$Image_Metadata_array==10,]
library(plyr) 
icamintim.neg<-ddply(perobindall.corn,"ImageNumber", summarise,  
                 FeatureIdx=unique(FeatureIdx),
                 Array=unique(Image_Metadata_array),
                 #Icam
                 ## noncorrected intensity 
                 #median
                 IcamIntensityMeanmed=mean(Cell_Intensity_MedianIntensity_Icam1amask,trim=0.2),
                 IcamIntensityMedmed=median(Cell_Intensity_MedianIntensity_Icam1amask),
                 ##mean
                 IcamIntensityMeanmean=mean(Cell_Intensity_MeanIntensity_Icam1amask,trim=0.2),
                 IcamIntensityMedmean=median(Cell_Intensity_MeanIntensity_Icam1amask),
                 ##integrated
                 IcamIntensityMeanint=mean(Cell_Intensity_IntegratedIntensity_Icam1amask,trim=0.2),
                 IcamIntensityMedint=median(Cell_Intensity_IntegratedIntensity_Icam1amask),
                 
                 ## corrected intensity min intensity subtracted
                 #median
                 IcamIntensityMeanmed_Nor=mean(Cell_Intensity_MedianIntensity_Icam1amask_Nor,trim=0.2),
                 IcamIntensityMedmed_Nor=median(Cell_Intensity_MedianIntensity_Icam1amask_Nor),
                 ##mean
                 IcamIntensityMeanmean_Nor=mean(Cell_Intensity_MeanIntensity_Icam1amask_Nor,trim=0.2),
                 IcamIntensityMedmean_Nor=median(Cell_Intensity_MeanIntensity_Icam1amask_Nor),
                 ##integrated
                 IcamIntensityMeanint_Nor=mean(Cell_Intensity_IntegratedIntensity_Icam1amask_Nor,trim=0.2),
                 IcamIntensityMedint_Nor=median(Cell_Intensity_IntegratedIntensity_Icam1amask_Nor),
                                  
                 ## corrected intensity per batch correction
                 #median
                 IcamIntensityMeanmed_Nor_corr=mean(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr,trim=0.2),
                 IcamIntensityMedmed_Nor_corr=median(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr),
                 ##mean
                 IcamIntensityMeanmean_Nor_corr=mean(Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr,trim=0.2),
                 IcamIntensityMedmean_Nor_corr=median(Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr),
                 ##integrated
                 IcamIntensityMeanint_Nor_corr=mean(Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr,trim=0.2),
                 IcamIntensityMedint_Nor_corr=median(Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr),
                 
                 #Actin
                 ## noncorrected intensity 
                 #median
                 ActinIntensityMeanmed=mean(Cell_Intensity_MedianIntensity_Actin1amask,trim=0.2),
                 ActinIntensityMedmed=median(Cell_Intensity_MedianIntensity_Actin1amask),
                 ##mean
                 ActinIntensityMeanmean=mean(Cell_Intensity_MeanIntensity_Actin1amask,trim=0.2),
                 ActinIntensityMedmean=median(Cell_Intensity_MeanIntensity_Actin1amask),
                 ##integrated
                 ActinIntensityMeanint=mean(Cell_Intensity_IntegratedIntensity_Actin1amask,trim=0.2),
                 ActinIntensityMedint=median(Cell_Intensity_IntegratedIntensity_Actin1amask),
                 
                 ## corrected intensity min intensity subtracted
                 #median
                 ActinIntensityMeanmed_Nor=mean(Cell_Intensity_MedianIntensity_Actin1amask_Nor,trim=0.2),
                 ActinIntensityMedmed_Nor=median(Cell_Intensity_MedianIntensity_Actin1amask_Nor),
                 ##mean
                 ActinIntensityMeanmean_Nor=mean(Cell_Intensity_MeanIntensity_Actin1amask_Nor,trim=0.2),
                 ActinIntensityMedmean_Nor=median(Cell_Intensity_MeanIntensity_Actin1amask_Nor),
                 ##integrated
                 ActinIntensityMeanint_Nor=mean(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor,trim=0.2),
                 ActinIntensityMedint_Nor=median(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor),
                 
                 ## corrected intensity per batch correction
                 #median
                 ActinIntensityMeanmed_Nor_corr=mean(Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr,trim=0.2),
                 ActinIntensityMedmed_Nor_corr=median(Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr),
                 ##mean
                 ActinIntensityMeanmean_Nor_corr=mean(Cell_Intensity_MeanIntensity_Actin1amask_Nor_corr,trim=0.2),
                 ActinIntensityMedmean_Nor_corr=median(Cell_Intensity_MeanIntensity_Actin1amask_Nor_corr),
                 ##integrated
                 ActinIntensityMeanint_Nor_corr=mean(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor_corr,trim=0.2),
                 ActinIntensityMedint_Nor_corr=median(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor_corr))


# #make pairs plot  
# library(GGally)
# ggpairs(icamintim,columns=4:27)
##collapsing to feature number                 

icamintft<-ddply(icamintim,"FeatureIdx",numcolwise(function(x) mean(x, trimm=0.2)))

plot(icamintft[order(icamintft$IcamIntensityMeanmed),"IcamIntensityMeanmed"])
plot(icamintft[order(icamintft$IcamIntensityMeanmed_Nor),"IcamIntensityMeanmed_Nor"])
plot(icamintft[order(icamintft$IcamIntensityMeanmed_Nor_corr),"IcamIntensityMeanmed_Nor_corr"])
plot(icamintft$IcamIntensityMeanmed_Nor_corr,icamintft$IcamIntensityMeanmed)
plot(icamintft$IcamIntensityMeanmed_Nor_corr,icamintft$IcamIntensityMeanmed_Nor)
plot(icamintft$IcamIntensityMeanmed,icamintft$IcamIntensityMeanmed_Nor)

# running statistical test to find surfaces from hits that are statistically different.
for (i in unique(icamintft$FeatureIdx)){ 
  temp <- icamintim[icamintim$FeatureIdx==i,]
  pva<-wilcox.test(temp[,"IcamIntensityMeanmed_Nor_corr"],icamintim.neg[,"IcamIntensityMeanmed_Nor_corr"] )
  icamintft[i,"p.value"]<-as.numeric(pva["p.value"])
}
##select only surfaces that passess p value + correction
## BH correction let pass only 2 high Icam intensity surfaces 18 in total
icamintft$p.value.adj<-p.adjust(icamintft$p.value, method= "none")
icamintft.pv<-icamintft[icamintft$p.value.adj<0.05,]
rankration<-icamintft.pv[order(icamintft.pv$IcamIntensityMeanmed_Nor_corr,
                               icamintft.pv$p.value.adj),c("FeatureIdx","IcamIntensityMeanmed_Nor_corr")]
# count total number of passed surfaces
nrow(rankration)
plot(rankration$IcamIntensityMeanmed_Nor_corr)

bottomicam<-head(rankration, n=20)
topicam<-tail(rankration, n=20)

save(bottomicam,topicam,file="hit's based on icam intensities.RDATA")

topicamm<-merge(topicam,icamintim, by="FeatureIdx", sort=F)
bottomicamm<-merge(bottomicam,icamintim, by="FeatureIdx", sort=F)

forplotcoll<-rbind(bottomicamm,topicamm)

forplottransf <- transform(forplotcoll[,c("FeatureIdx", "IcamIntensityMeanmed_Nor_corr.y",
                                          "ActinIntensityMeanmed_Nor_corr") ], FeatureIdx = factor(FeatureIdx, 
     levels = unique(as.character(forplotcoll$FeatureIdx))))
library(reshape2)
forplottransfmelt<-melt(forplottransf, measure.vars =c("IcamIntensityMeanmed_Nor_corr.y",
                                                       "ActinIntensityMeanmed_Nor_corr"))
library(ggplot2)
ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))+
  geom_boxplot()+theme(legend.position="none")+ylim(0,15)
##Repeat analysis on normolized Icam intensity to actin intensity



icamintim$IcamIntensityMeanmed_Nor_corr_norm<-icamintim$IcamIntensityMeanmed_Nor_corr/icamintim$ActinIntensityMeanmed_Nor_corr
icamintim.neg$IcamIntensityMeanmed_Nor_corr_norm<-icamintim.neg$IcamIntensityMeanmed_Nor_corr/icamintim.neg$ActinIntensityMeanmed_Nor_corr
summary(icamintim$IcamIntensityMeanmed_Nor_corr_norm)

icamintft<-ddply(icamintim,"FeatureIdx",numcolwise(function(x) mean(x, trimm=0.2)))
# running statistical test to find surfaces from hits that are statistically different.
for (i in unique(icamintft$FeatureIdx)){ 
  temp <- icamintim[icamintim$FeatureIdx==i,]
  pva<-wilcox.test(temp[,"IcamIntensityMeanmed_Nor_corr_norm"],icamintim.neg[,"IcamIntensityMeanmed_Nor_corr_norm"] )
  icamintft[i,"p.value"]<-as.numeric(pva["p.value"])
}
##select only surfaces that passess p value + correction
## BH correction let pass only 2 high Icam intensity surfaces 18 in total
icamintft$p.value.adj<-p.adjust(icamintft$p.value, method= "none")
icamintft.pv<-icamintft[icamintft$p.value.adj<0.05,]
rankration<-icamintft.pv[order(icamintft.pv$IcamIntensityMeanmed_Nor_corr_norm,
                               icamintft.pv$p.value.adj),c("FeatureIdx","IcamIntensityMeanmed_Nor_corr_norm")]
# count total number of passed surfaces
nrow(rankration)
plot(rankration$IcamIntensityMeanmed_Nor_corr_norm)

bottomicam<-head(rankration, n=20)
topicam<-tail(rankration, n=20)
save(bottomicam,topicam,file="hit's based on normolized icam intensities.RDATA")
topicamm<-merge(topicam,icamintim, by="FeatureIdx", sort=F)
bottomicamm<-merge(bottomicam,icamintim, by="FeatureIdx", sort=F)

forplotcoll<-rbind(bottomicamm,topicamm)

forplottransf <- transform(forplotcoll[,c("FeatureIdx", "IcamIntensityMeanmed_Nor_corr",
                                          "ActinIntensityMeanmed_Nor_corr") ], FeatureIdx = factor(FeatureIdx, 
                                                                                                   levels = unique(as.character(forplotcoll$FeatureIdx))))
library(reshape2)
forplottransfmelt<-melt(forplottransf, measure.vars =c("IcamIntensityMeanmed_Nor_corr",
                                                       "ActinIntensityMeanmed_Nor_corr"))
library(ggplot2)
ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))+
  geom_boxplot()+theme(legend.position="none")+ylim(0,15)

##Integrated Intensity

# running statistical test to find surfaces from hits that are statistically different.
for (i in unique(icamintft$FeatureIdx)){ 
  temp <- icamintim[icamintim$FeatureIdx==i,]
  pva<-wilcox.test(temp[,"IcamIntensityMeanint_Nor_corr"],icamintim.neg[,"IcamIntensityMeanint_Nor_corr"] )
  icamintft[i,"p.value"]<-as.numeric(pva["p.value"])
}
##select only surfaces that passess p value + correction
## BH correction let pass only 2 high Icam intensity surfaces 18 in total
icamintft$p.value.adj<-p.adjust(icamintft$p.value, method= "BH")
icamintft.pv<-icamintft[icamintft$p.value.adj<0.05,]
rankration<-icamintft.pv[order(icamintft.pv$IcamIntensityMeanint_Nor_corr,
                               icamintft.pv$p.value.adj),c("FeatureIdx","IcamIntensityMeanint_Nor_corr")]
# count total number of passed surfaces
nrow(rankration)
plot(rankration$IcamIntensityMeanint_Nor_corr)

bottomicam<-head(rankration, n=20)
topicam<-tail(rankration, n=20)
save(bottomicam,topicam,file="hit's based on icam integrated intensities.RDATA")
topicamm<-merge(topicam,icamintim, by="FeatureIdx", sort=F)
bottomicamm<-merge(bottomicam,icamintim, by="FeatureIdx", sort=F)

forplotcoll<-rbind(bottomicamm,topicamm)

forplottransf <- transform(forplotcoll[,c("FeatureIdx", "IcamIntensityMeanint_Nor_corr.y",
                                          "ActinIntensityMeanint_Nor_corr") ], FeatureIdx = factor(FeatureIdx, 
                                                                                                   levels = unique(as.character(forplotcoll$FeatureIdx))))
library(reshape2)
forplottransfmelt<-melt(forplottransf, measure.vars =c("IcamIntensityMeanint_Nor_corr.y",
                                                       "ActinIntensityMeanint_Nor_corr"))
library(ggplot2)
ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))+
  geom_boxplot()+theme(legend.position="none")
