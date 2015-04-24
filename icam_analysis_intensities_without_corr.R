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
                 ActinIntensityMedint_Nor=median(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor))
                 

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
                     ActinIntensityMedint_Nor=median(Cell_Intensity_IntegratedIntensity_Actin1amask_Nor))
                     


# #make pairs plot  
# library(GGally)
# ggpairs(icamintim,columns=4:27)
##collapsing to feature number                 

icamintft<-ddply(icamintim,"FeatureIdx",numcolwise(function(x) mean(x, trimm=0.2)))

plot(icamintft[order(icamintft$IcamIntensityMeanmed),"IcamIntensityMeanmed"])
plot(icamintft[order(icamintft$IcamIntensityMeanmed_Nor),"IcamIntensityMeanmed_Nor"])
plot(icamintft[order(icamintft$IcamIntensityMeanmed_Nor),"IcamIntensityMeanmed_Nor"])
plot(icamintft$IcamIntensityMeanmed_Nor,icamintft$IcamIntensityMeanmed)
plot(icamintft$IcamIntensityMeanmed_Nor,icamintft$IcamIntensityMeanmed_Nor)
plot(icamintft$IcamIntensityMeanmed,icamintft$IcamIntensityMeanmed_Nor)

# running statistical test to find surfaces from hits that are statistically different.
for (i in unique(icamintft$FeatureIdx)){ 
  temp <- icamintim[icamintim$FeatureIdx==i,]
  pva<-wilcox.test(temp[,"IcamIntensityMeanmed_Nor"],icamintim.neg[,"IcamIntensityMeanmed_Nor"] )
  icamintft[i,"p.value"]<-as.numeric(pva["p.value"])
}
##select only surfaces that passess p value + correction
## BH correction let pass only 2 high Icam intensity surfaces 18 in total
icamintft$p.value.adj<-p.adjust(icamintft$p.value, method= "none")
icamintft.pv<-icamintft[icamintft$p.value.adj<0.05,]
rankration<-icamintft.pv[order(icamintft.pv$IcamIntensityMeanmed_Nor),c("FeatureIdx","IcamIntensityMeanmed_Nor")]
# count total number of passed surfaces
nrow(rankration)
plot(rankration$IcamIntensityMeanmed_Nor)

if(nrow(rankration)>=40) number.surf=20 else number.surf=round(nrow(rankration)/2)

bottomicam<-head(rankration, n=number.surf)
topicam<-tail(rankration, n=number.surf)

topicamm<-merge(topicam,icamintim, by="FeatureIdx", sort=F)
topicamm$Hit<-"Top"
bottomicamm<-merge(bottomicam,icamintim, by="FeatureIdx", sort=F)
bottomicamm$Hit<-"Bottom"
forplotcoll<-rbind(bottomicamm,topicamm)
save(forplotcoll,file="hit's based on icam intensities_without_correct.RDATA")
#save(forplotcoll,file="hit's based on icam intensities_with_correct.RDATA")

forplottransf <- transform(forplotcoll[,c("FeatureIdx", "IcamIntensityMeanmed_Nor.y",
                                          "ActinIntensityMeanmed_Nor") ], FeatureIdx = factor(FeatureIdx, 
                                                                                                   levels = unique(as.character(forplotcoll$FeatureIdx))))
library(reshape2)
forplottransfmelt<-melt(forplottransf, measure.vars =c("IcamIntensityMeanmed_Nor.y",
                                                       "ActinIntensityMeanmed_Nor"))
# ##multiply by 255
# forplottransfmelt$value<-forplottransfmelt$value*255

library(ggplot2)
ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))+
  geom_boxplot()+theme(legend.position="none")+ylim(0.01,0.045)#+geom_jitter()
#############################################################################
#########################################################################

##Repeat analysis on normolized Icam intensity to actin intensity



icamintim$IcamIntensityMeanmed_Nor_norm<-icamintim$IcamIntensityMeanmed_Nor/icamintim$ActinIntensityMeanmed_Nor
icamintim.neg$IcamIntensityMeanmed_Nor_norm<-icamintim.neg$IcamIntensityMeanmed_Nor/icamintim.neg$ActinIntensityMeanmed_Nor
summary(icamintim$IcamIntensityMeanmed_Nor_norm)

icamintft<-ddply(icamintim,"FeatureIdx",numcolwise(function(x) mean(x, trimm=0.2)))
# running statistical test to find surfaces from hits that are statistically different.
for (i in unique(icamintft$FeatureIdx)){ 
  temp <- icamintim[icamintim$FeatureIdx==i,]
  pva<-wilcox.test(temp[,"IcamIntensityMeanmed_Nor_norm"],icamintim.neg[,"IcamIntensityMeanmed_Nor_norm"] )
  icamintft[i,"p.value"]<-as.numeric(pva["p.value"])
}
##select only surfaces that passess p value + correction
## BH correction let pass only 2 high Icam intensity surfaces 18 in total
icamintft$p.value.adj<-p.adjust(icamintft$p.value, method= "BH")
icamintft.pv<-icamintft[icamintft$p.value.adj<0.05,]
rankration<-icamintft.pv[order(icamintft.pv$IcamIntensityMeanmed_Nor_norm,
                               icamintft.pv$p.value.adj),c("FeatureIdx","IcamIntensityMeanmed_Nor_norm")]
# count total number of passed surfaces
nrow(rankration)
plot(rankration$IcamIntensityMeanmed_Nor_norm)

bottomicam<-head(rankration, n=20)
topicam<-tail(rankration, n=20)
save(bottomicam,topicam,file="hit's based on normolized icam intensities.RDATA")
topicamm<-merge(topicam,icamintim, by="FeatureIdx", sort=F)
topicamm$Hit<-"Top"
bottomicamm<-merge(bottomicam,icamintim, by="FeatureIdx", sort=F)
bottomicamm$Hit<-"Bottom"
forplotcoll<-rbind(bottomicamm,topicamm)

forplottransf <- transform(forplotcoll[,c("FeatureIdx", "IcamIntensityMeanmed_Nor",
                                          "ActinIntensityMeanmed_Nor") ], FeatureIdx = factor(FeatureIdx, 
                                                                                                   levels = unique(as.character(forplotcoll$FeatureIdx))))
library(reshape2)
forplottransfmelt<-melt(forplottransf, measure.vars =c("IcamIntensityMeanmed_Nor",
                                                       "ActinIntensityMeanmed_Nor"))
library(ggplot2)
ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))+
  geom_boxplot()+theme(legend.position="none")+ylim(0,0.05)

##Integrated Intensity

# running statistical test to find surfaces from hits that are statistically different.
for (i in unique(icamintft$FeatureIdx)){ 
  temp <- icamintim[icamintim$FeatureIdx==i,]
  pva<-wilcox.test(temp[,"IcamIntensityMeanint_Nor"],icamintim.neg[,"IcamIntensityMeanint_Nor"] )
  icamintft[i,"p.value"]<-as.numeric(pva["p.value"])
}
##select only surfaces that passess p value + correction
## BH correction let pass only 2 high Icam intensity surfaces 18 in total
icamintft$p.value.adj<-p.adjust(icamintft$p.value, method= "BH")
icamintft.pv<-icamintft[icamintft$p.value.adj<0.05,]
rankration<-icamintft.pv[order(icamintft.pv$IcamIntensityMeanint_Nor,
                               icamintft.pv$p.value.adj),c("FeatureIdx","IcamIntensityMeanint_Nor")]
# count total number of passed surfaces
nrow(rankration)
plot(rankration$IcamIntensityMeanint_Nor)

bottomicam<-head(rankration, n=20)
topicam<-tail(rankration, n=20)

topicamm<-merge(topicam,icamintim, by="FeatureIdx", sort=F)
topicamm$Hit<-"Top"
bottomicamm<-merge(bottomicam,icamintim, by="FeatureIdx", sort=F)
bottomicamm$Hit<-"Bottom"
forplotcoll<-rbind(bottomicamm,topicamm)

save(forplotcoll,file="hit's based on icam integrated intensities.RDATA")



forplottransf <- transform(forplotcoll[,c("FeatureIdx", "IcamIntensityMeanint_Nor.y",
                                          "ActinIntensityMeanint_Nor") ], FeatureIdx = factor(FeatureIdx, 
                                                                                                   levels = unique(as.character(forplotcoll$FeatureIdx))))
library(reshape2)
forplottransfmelt<-melt(forplottransf, measure.vars =c("IcamIntensityMeanint_Nor.y",
                                                       "ActinIntensityMeanint_Nor"))
library(ggplot2)
ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))+
  geom_boxplot()+theme(legend.position="none")+ylim(0,600)

###############ribbon ggplot
##another aproach without outliers
library(ggplot2)
library(reshape2)
icam_intensities.f.stats<-ddply(icamintim,"FeatureIdx", summarise, 
                                      ymin=boxplot.stats(IcamIntensityMedmed)$stats[1],
                                      lower=boxplot.stats(IcamIntensityMedmed)$stats[2],
                                      middle=boxplot.stats(IcamIntensityMedmed)$stats[3],
                                      upper=boxplot.stats(IcamIntensityMedmed)$stats[4],
                                      ymax=boxplot.stats(IcamIntensityMedmed)$stats[5])
##ordering
icam_intensities.f.stats<-icam_intensities.f.stats[order(icam_intensities.f.stats$middle),]
icam_intensities.f.stats$FeatureIdx<-factor(icam_intensities.f.stats$FeatureIdx,
                  levels=as.character(icam_intensities.f.stats$FeatureIdx))
# create plot with all statistics 
ggplot(icam_intensities.f.stats, aes(FeatureIdx, 
                                            lower=lower, upper=upper, middle=middle, ymin=ymin, ymax=ymax)) + 
  geom_boxplot(stat="identity", fill = "#E69F00",colour = "#0072B2")+
  geom_point(data=icam_intensities.f.stats, aes(x=c(1:2177),y=middle),
             colour="red",shape="-")+ylab("Median ICAM Intensity per repeat")

# create plot with all statistics ribbon

ggplot(icam_intensities.f.stats, aes(x=FeatureIdx))+
  geom_ribbon(data=icam_intensities.f.stats,
              aes(group = 1,ymin=ymin, ymax=ymax,fill="blue",colour="blue"))+
  geom_ribbon(data=icam_intensities.f.stats,
              aes(group = 1,ymin=lower, ymax=upper,fill="yellow",colour="yellow"))+
  geom_line(data=icam_intensities.f.stats,
            aes(group = 1,y=middle,colour="red"))+
  scale_fill_manual(values=c("blue"="blue","yellow"="yellow","red"="red"))+
  scale_colour_manual(values=c("blue"="blue","yellow"="yellow","red"="red"))+
  ylab("Median ICAM Intensity per repeat")

