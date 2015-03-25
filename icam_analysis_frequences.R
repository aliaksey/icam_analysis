rm(list=ls())
load("Icam_cells_after_correction.RDATA")

#find treshhold factor
library(ggplot2)
library(reshape2)
#based on mean
top_int_data.mean<-perobindall.cor[perobindall.cor$Image_Metadata_array<9,"Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr"]
pos_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==9,"Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr"]
neg_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==10,"Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr"]
lower.limit <- 1
upper.limit <- 50
neg.density <- density(neg_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
pos.density <- density(pos_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
density.difference <- neg.density$y - pos.density$y
intersection.point.mean <- neg.density$x[which(diff(density.difference > 0) != 0) + 1]

combined<-perobindall.cor[perobindall.cor$Image_Metadata_array>8,c("Image_Metadata_array",
                  "Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr")]
ggplot(combined, aes(Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr, fill =as.factor(Image_Metadata_array))) +
  geom_density(alpha = 0.2) + 
  geom_vline(xintercept = intersection.point.mean, color = "red")
#based on median
top_int_data.median<-perobindall.cor[perobindall.cor$Image_Metadata_array<9,"Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr"]
pos_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==9,"Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr"]
neg_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==10,"Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr"]
lower.limit <- 1
upper.limit <- 50
neg.density <- density(neg_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
pos.density <- density(pos_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
density.difference <- neg.density$y - pos.density$y
intersection.point <- neg.density$x[which(diff(density.difference > 0) != 0) + 1]
intersection.point.median<-max(intersection.point)

combined<-perobindall.cor[perobindall.cor$Image_Metadata_array>8,c("Image_Metadata_array",
                                                                   "Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr")]
ggplot(combined, aes(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr, fill =as.factor(Image_Metadata_array))) +
  geom_density(alpha = 0.2) + 
  geom_vline(xintercept = intersection.point.median, color = "red")
#based on integrated
top_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array<9,"Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr"]
pos_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==9,"Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr"]
neg_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==10,"Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr"]
lower.limit <- 1
upper.limit <- 500000
neg.density <- density(neg_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
pos.density <- density(pos_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
density.difference <- neg.density$y - pos.density$y
intersection.point <- neg.density$x[which(diff(density.difference > 0) != 0) + 1]

combined<-perobindall.cor[perobindall.cor$Image_Metadata_array>8,c("Image_Metadata_array",
                                                                   "Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr")]
ggplot(combined, aes(Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr, fill =as.factor(Image_Metadata_array))) +
  geom_density(alpha = 0.2) + 
  geom_vline(xintercept = intersection.point, color = "red")
#based on mean non corrected
top_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array<9,"Cell_Intensity_MeanIntensity_Icam1amask"]
pos_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==9,"Cell_Intensity_MeanIntensity_Icam1amask"]
neg_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==10,"Cell_Intensity_MeanIntensity_Icam1amask"]
lower.limit <- 0.3
upper.limit <- 0.5
neg.density <- density(neg_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
pos.density <- density(pos_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
density.difference <- neg.density$y - pos.density$y
intersection.point <- neg.density$x[which(diff(density.difference > 0) != 0) + 1]

combined<-perobindall.cor[perobindall.cor$Image_Metadata_array>8,c("Image_Metadata_array",
                                                                   "Cell_Intensity_MeanIntensity_Icam1amask")]
ggplot(combined, aes(Cell_Intensity_MeanIntensity_Icam1amask, fill =as.factor(Image_Metadata_array))) +
  geom_density(alpha = 0.2) + 
  geom_vline(xintercept = intersection.point, color = "red")
#based on mean non corrected only min subtracted
top_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array<9,"Cell_Intensity_MeanIntensity_Icam1amask_Nor"]
pos_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==9,"Cell_Intensity_MeanIntensity_Icam1amask_Nor"]
neg_int_data<-perobindall.cor[perobindall.cor$Image_Metadata_array==10,"Cell_Intensity_MeanIntensity_Icam1amask_Nor"]
lower.limit <- 0.02
upper.limit <- 0.1
neg.density <- density(neg_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
pos.density <- density(pos_int_data, from = lower.limit, to = upper.limit, n = max(length(neg_int_data),
                                                                                   length(pos_int_data)))
density.difference <- neg.density$y - pos.density$y
intersection.point <- neg.density$x[which(diff(density.difference > 0) != 0) + 1]

combined<-perobindall.cor[perobindall.cor$Image_Metadata_array>8,c("Image_Metadata_array",
                                                                   "Cell_Intensity_MeanIntensity_Icam1amask_Nor")]
ggplot(combined, aes(Cell_Intensity_MeanIntensity_Icam1amask_Nor, fill =as.factor(Image_Metadata_array))) +
  geom_density(alpha = 0.2) + 
  geom_vline(xintercept = intersection.point, color = "red")
###########################################
###finding treshhold based on Topochip data
int_data.topochip<-perobindall.cor[perobindall.cor$Image_Metadata_array<12,
      c("Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr","Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr",
        "Image_Metadata_array")]
library(reshape2)
int_data.topochip.melt<-melt(int_data.topochip,measure.vars =c("Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr",
                                "Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr"))
#mean                                                              
ggplot(int_data.topochip.melt, aes(as.factor(variable),y=value,fill=as.factor(Image_Metadata_array))) + 
  geom_boxplot()+ scale_y_log10()+ geom_hline(yintercept = intersection.point.mean, color = "red")+
  ggtitle("Mean Icam expression")

ggplot(int_data.topochip.melt[int_data.topochip.melt$value>intersection.point.mean,], aes(as.factor(variable),y=value,fill=as.factor(Image_Metadata_array))) + 
  geom_boxplot()+ scale_y_log10()+ geom_hline(yintercept = intersection.point.mean, color = "red")+
  ggtitle("Mean Icam expression after treshhold")
#median
ggplot(int_data.topochip.melt, aes(as.factor(variable),y=value,fill=as.factor(Image_Metadata_array))) + 
  geom_boxplot()+ scale_y_log10()+ geom_hline(yintercept = intersection.point.median, color = "red")+
  ggtitle("Median Icam expression")

ggplot(int_data.topochip.melt[int_data.topochip.melt$value>intersection.point.median,], aes(as.factor(variable),y=value,fill=as.factor(Image_Metadata_array))) + 
  geom_boxplot()+ scale_y_log10()+ geom_hline(yintercept = intersection.point.median, color = "red")+
  ggtitle("Median Icam expression after treshhold")



###############################################################
#################################################################
##ICAM analysis

#count number of postitive cells, with ICAm expression higher then treshhold

top_int_data.median<-perobindall.cor[perobindall.cor$Image_Metadata_array<9,
      c("Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr","ImageNumber","FeatureIdx")]

neg_int_data.median<-perobindall.cor[perobindall.cor$Image_Metadata_array==10,
                                     c("Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr","ImageNumber","FeatureIdx")]



library(plyr)
ratiorankcorr<-ddply(top_int_data.median,"ImageNumber", summarise, 
                     IcamPositive=sum(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr > intersection.point.median),
                     Total=sum(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr > 0),
                     TrMeanIcam=mean(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr, trimm=0.2),
                     FeatureIdx=unique(FeatureIdx))

ratiorankcorrn<-ddply(neg_int_data.median,"ImageNumber", summarise, 
                      IcamPositive=sum(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr > intersection.point.median),
                      Total=sum(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr > 0),
                      TrMeanIcam=mean(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr, trimm=0.2),
                      FeatureIdx=unique(FeatureIdx))
ratiorankcorr$Ratio<-ratiorankcorr$IcamPositive/ratiorankcorr$Total
ratiorankcorrn$Ratio<-ratiorankcorrn$IcamPositive/ratiorankcorrn$Total

##show number of cells per repeat
median(ratiorankcorr$Total)


ggplot(ratiorankcorr,aes(x=Ratio,y=TrMeanIcam,colour=as.factor(FeatureIdx)))+geom_point(show_guide = FALSE)

#collapse to feaature nuber
ratiorankcorrfff<-ddply(na.omit(ratiorankcorr),"FeatureIdx", summarise, 
                        IcamPos=sum(IcamPositive),
                        Totalc=sum(Total),
                        RatioTrMean=mean(Ratio,trimm=0.2),
                        RatioMedian=median(Ratio),
                        IntensTrMean=mean(TrMeanIcam,trimm=0.2),
                        IntensTrMad=mad(TrMeanIcam),
                        RatioSd=sd(Ratio),
                        RatioMad=mad(Ratio))

ratiorankcorrfff$RatioAbs<-ratiorankcorrfff$IcamPos/ratiorankcorrfff$Totalc
ratiorankcorrfffsm<-ratiorankcorrfff[order(ratiorankcorrfff$RatioMedian),]

ggplot(ratiorankcorrfffsm,aes(x=Featureidx<-c(1:2177),y=RatioMedian,colour=as.factor(FeatureIdx)))+
  geom_point(show_guide = FALSE)+ geom_hline(yintercept = mean(ratiorankcorrn$Ratio), color = "red")+
  ggtitle("Rattio of ICAm positive cells per topounit/calculated from Median")

ratiorankcorrfffs<-ratiorankcorrfff[order(ratiorankcorrfff$RatioTrMean),]

scurve<-ggplot(ratiorankcorrfffs,aes(x=Featureidx<-c(1:2177),y=RatioTrMean,colour=as.factor(FeatureIdx)))+
  geom_point(show_guide = FALSE)+ geom_hline(yintercept = mean(ratiorankcorrn$Ratio), color = "red")+
  geom_hline(yintercept = mean(ratiorankcorr$Ratio), color = "blue")+
  ggtitle("Rattio of ICAm positive cells per topounit/calculated from Mean")
scurve
##scurve +error
library(reshape2)
ratiorankcorr.s<-merge(ratiorankcorrfffs,ratiorankcorr,sort=F)
ratiorankcorr_scurve.tr <- transform(ratiorankcorr.s[,c("FeatureIdx", "Ratio") ], FeatureIdx = factor(FeatureIdx,
levels = unique(as.character(ratiorankcorr.s$FeatureIdx))))
ratiorankcorr_scurve<-melt(ratiorankcorr_scurve.tr, measure.vars ="Ratio")

ggplot(ratiorankcorr_scurve,aes(FeatureIdx,value))+
  geom_point()+ylim(0,0.2)+
geom_smooth(aes(group=1), method="lm")
  

plot_top <- ggplot(ratiorankcorrfffs, aes(RatioTrMean,fill="yellow")) + 
  geom_density(alpha=.5) + 
  theme(legend.position = "none") 
fordensplot<-rbind(data.frame(Intensity=ratiorankcorrn$Ratio,Sample="Negative"),
                   data.frame(Intensity=ratiorankcorr$Ratio,Sample="Topogra"))
plot_right <- ggplot(ratiorankcorrfffs, aes(RatioTrMean,fill="yellow")) + 
  geom_density(alpha=.5) + 
  coord_flip() + 
  theme(legend.position = "none")
empt <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )
library(gridExtra)
grid.arrange(plot_top, empt, scurve, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

##calculate if mean icam significantly different between topography  and negative control

wilcox.test(ratiorankcorr$Ratio,ratiorankcorrn$Ratio, correct=F)
##median of topo=0, not suprisingly that u-test is very significant
##t-test is more applicable
t.test(ratiorankcorr$Ratio,ratiorankcorrn$Ratio)

##thes same plot but with absolute ratio
ratiorankcorrfffsa<-ratiorankcorrfff[order(ratiorankcorrfff$RatioAbs),]

ggplot(ratiorankcorrfffsa,aes(x=Featureidx<-c(1:2177),y=RatioAbs,colour=as.factor(FeatureIdx)))+
  geom_point(show_guide = FALSE)+ geom_hline(yintercept = mean(ratiorankcorrn$Ratio), color = "red")+
  geom_hline(yintercept = mean(ratiorankcorr$Ratio), color = "blue")+
  ggtitle("Rattio of ICAm positive cells per topounit/calculated from Mean")
##correlation between absolute Ratio and mean ratio
ggplot(ratiorankcorrfff, aes(x=RatioTrMean,y=RatioAbs,colour=as.factor(FeatureIdx)))+
  geom_smooth(method = "lm", colour = "grey50", fill = "grey50") + 
  geom_point(show_guide = FALSE) 
  
  
cor(ratiorankcorrfff$RatioTrMean,ratiorankcorrfff$RatioAbs)

##plots with intensities
ggplot(ratiorankcorrfffs,aes(x=Featureidx<-c(1:2177),y=IntensTrMean,colour=as.factor(FeatureIdx)))+
  geom_point(show_guide = FALSE)+ geom_hline(yintercept = mean(ratiorankcorrn$TrMeanIcam), color = "red")+
  ggtitle("ICAM Intensities  per topounit/sorted by persentage of positive")

ratiorankcorrfffsi<-ratiorankcorrfff[order(ratiorankcorrfff$IntensTrMean),]

ggplot(ratiorankcorrfffsi,aes(x=Featureidx<-c(1:2177),y=IntensTrMean,colour=as.factor(FeatureIdx)))+
  geom_point(show_guide = FALSE)+ geom_hline(yintercept = mean(ratiorankcorrn$TrMeanIcam), color = "red")+
  geom_hline(yintercept = mean(ratiorankcorr$TrMeanIcam), color = "blue")+
  ggtitle("ICAM Intensities  per topounit/sorted by mean intensity")
ggplot(ratiorankcorrfffsi,aes(IntensTrMean,fill="yellow"))+geom_density(alpha=.5)

###so the same with actin to show that Icam upregulation is not random



##Actin analysis


topa_inta_data.median<-perobindall.cor[perobindall.cor$Image_Metadata_array<9,
                                     c("Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr","ImageNumber","FeatureIdx")]
neg_inta_data.median<-perobindall.cor[perobindall.cor$Image_Metadata_array==10,
                                     c("Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr","ImageNumber","FeatureIdx")]
library(plyr)
ratiorankcorra<-ddply(topa_inta_data.median,"ImageNumber", summarise, 
                     TrMeanActin=mean(Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr, trimm=0.2),
                     FeatureIdx=unique(FeatureIdx))

ratiorankcorran<-ddply(neg_inta_data.median,"ImageNumber", summarise, 
                      TrMeanActin=mean(Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr, trimm=0.2),
                      FeatureIdx=unique(FeatureIdx))
#collapse to feaature nuber
ratiorankcorrafff<-ddply(na.omit(ratiorankcorra),"FeatureIdx", summarise, 
                        IntensTrMean=mean(TrMeanActin),
                        IntensTrMad=mad(TrMeanActin))
ratiorankcorrafffsi<-ratiorankcorrafff[order(ratiorankcorrafff$IntensTrMean),]
ggplot(ratiorankcorrafffsi,aes(x=Featureidx<-c(1:2177),y=IntensTrMean,colour=as.factor(FeatureIdx)))+
  geom_point(show_guide = FALSE)+ geom_hline(yintercept = mean(ratiorankcorran$TrMeanActin), color = "red")+
  geom_hline(yintercept = mean(ratiorankcorra$TrMeanActin), color = "blue")+
  ggtitle("Actin Intensities  per topounit/sorted by mean intensity")
ggplot(ratiorankcorrafffsi,aes(IntensTrMean,fill="yellow"))+geom_density(alpha=.5)
  
#############################################

#calculating statistics

for (i in 1:length(unique(ratiorankcorrfff[,"FeatureIdx"]))){ 
  temp <- ratiorankcorr[ratiorankcorr$FeatureIdx==ratiorankcorrfff[i,"FeatureIdx"],]
  #temp2<-perobindallfnt[perobindallfnt$FeatureIdx==collobfeat[i,"FeatureIdx"],]
  obs<-c(sum(temp[,"IcamPositive"]),(sum(temp[,"Total"])-sum(temp[,"IcamPositive"])))
  exp<-c(sum(neg_int_data.median$Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr >intersection.point.median)/nrow(neg_int_data.median),
         sum(neg_int_data.median$Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr<intersection.point.median)/nrow(neg_int_data.median))
  pva<-chisq.test(obs,p=exp)
  #wilcox.test(temp[,"IcamIntensityMt"],collob[,"IcamIntensityMt"] )
  #pvak<-kruskal.test(Cell_Intensity_MedianIntensity_Icam1amask~ImageNumber,data=temp2)
  #pvaint<-wilcox.test(temp[,"IcamINtegratedInt"],collob[,"IcamINtegratedInt"] )
  ratiorankcorrfff[i,"p.value"]<-as.numeric(pva["p.value"])
  #collobfeat[i,"int.p.value"]<-as.numeric(pvaint["p.value"])
  #collobfeat[i,"kruskal.p.value"]<-as.numeric(pvak["p.value"])
  #wilcox.test(perobindallfnt[perobindallfnt$FeatureIdx==1559,
  #                   "Cell_Intensity_MedianIntensity_Icam1amask"], 
  #    perobindallfnt$Cell_Intensity_MedianIntensity_Icam1amask)
}
sorted.pvalue<-ratiorankcorrfff[order(ratiorankcorrfff$p.value),]
ggplot(sorted.pvalue,aes(Featureidx<-c(1:2177),p.value))+geom_point()+geom_hline(yintercept=0.05,colour="red")+
  geom_hline(yintercept=0.05/2177,colour="blue")
ggplot(sorted.pvalue,aes(Featureidx<-c(1:2177),p.value))+geom_point()+geom_hline(yintercept=0.05,colour="red")+
  geom_hline(yintercept=0.05/2177,colour="blue")+ylim(0,0.05/2000)
length(ratiorankcorrfff[ratiorankcorrfff$p.value<(0.05/2177),"p.value"])
##select only surfaces that passess p value + correction
ratiorankcorrfff$p.value.adj<-p.adjust(ratiorankcorrfff$p.value, method= "BH")
ratiorankcorrfff.pv<-ratiorankcorrfff[ratiorankcorrfff$p.value.adj<0.05,]
rankration<-ratiorankcorrfff.pv[order(ratiorankcorrfff.pv$RatioTrMean,
                                      ratiorankcorrfff.pv$p.value.adj),c("FeatureIdx","RatioTrMean")]
# count total number of passed surfaces
nrow(ratiorankcorrfff.pv)
plot(rankration$RatioTrMean)

bottomicam<-head(rankration, n=20)
topicam<-tail(rankration, n=20)
########## plot hits
# plotting intensities for actin and icam

dataicamav<-ddply(perobindall.cor[perobindall.cor$Image_Metadata_array<9,c("FeatureIdx", "ImageNumber",
                                        "Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr", 
                "Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr")],"ImageNumber", summarise,
                FeatureIdx=unique(FeatureIdx),
                Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr=mean(Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr,
                                                                         trimm=0.2),
                Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr=mean(Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr,
                                                                         trimm=0.2))

topicamm<-merge(topicam,dataicamav, by="FeatureIdx", sort=F)
bottomicamm<-merge(bottomicam,dataicamav, by="FeatureIdx", sort=F)
forplotcoll<-rbind(bottomicamm,topicamm)

forplottransf <- transform(forplotcoll[,c("FeatureIdx", "Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr",
                    "Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr") ], FeatureIdx = factor(FeatureIdx, 
                                  levels = unique(as.character(forplotcoll$FeatureIdx))))
library(reshape2)
forplottransfmelt<-melt(forplottransf, measure.vars =c("Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr",
                                                       "Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr"))
#forplottransfmelt<-melt(forplottransf, measure.vars ="IcamIntensityMt")
ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))+
geom_boxplot()+theme(legend.position="none")+ylim(-0.5,15)

#####################plot frequences per selected surface

topicamm.fr<-merge(topicam,ratiorankcorr, by="FeatureIdx", sort=F)
topicamm.fr$Hit<-"top"
bottomicamm.fr<-merge(bottomicam,ratiorankcorr, by="FeatureIdx", sort=F)
bottomicamm.fr$Hit<-"bottom"
forplotcoll.fr<-rbind(bottomicamm.fr,topicamm.fr)
#plot all evrything without dots by mean
forplottransf.fr <- transform(forplotcoll.fr[,c("FeatureIdx", "Ratio","Hit")], FeatureIdx = factor(FeatureIdx, 
  levels = unique(as.character(forplotcoll.fr$FeatureIdx))))

forplottransfmelt.fr<-melt(forplottransf.fr, measure.vars ="Ratio")
#forplottransfmelt<-melt(forplottransf, measure.vars ="IcamIntensityMt")
library(doBy)
icam_freq_stat <- summaryBy(Ratio ~ FeatureIdx, forplotcoll.fr,order=F, FUN = c(mean, sd))

SD.quartile <- function(x){
  mean.temp<-mean(x)
  sd1.temp<-sd(x)
  n.temp<-length(x)
  CI.int<-0.995
  error.temp <- qnorm(CI.int)*sd1.temp/sqrt(n.temp)
  out <- c(mean.temp-error.temp,mean.temp,mean.temp+error.temp)
  names(out) <- c("ymin","y","ymax")
  return(out) 
}

#plot dots + confidence interval
ggplot(forplottransfmelt.fr, aes(FeatureIdx, value,group=Hit))+
  stat_summary(fun.data = SD.quartile, geom = "pointrange",color = "Blue")+
  geom_jitter( size = 2)  
  
#   geom_jitter(position = position_jitter(w = 0.1, h = 0.1), size = 1.5) +
#   layer(data = icam_freq_stat, mapping = aes(x = FeatureIdx, y = Ratio.mean, ymin = Ratio.mean - Ratio.sd,
#                                 ymax = Ratio.mean + Ratio.sd), geom = "errorbar", size = 0.9, color = "Blue",
#         width = 0.3) + layer(data = icam_freq_stat, mapping = aes(x = FeatureIdx, y = Ratio.mean),
#                              geom = "point", size = 8, color = "Blue", shape = "+")+
#   theme_bw() + xlab("") + ylab("Range of scores")


ggplot(forplottransfmelt.fr, aes(FeatureIdx, y = value, fill=Hit))+
geom_boxplot()+theme(legend.position="none")+geom_jitter()

# ggplot(forplottransfmelt.fr[forplottransfmelt.fr$FeatureIdx%in%c(532,1359,998,1652,263,628,572,201),], aes(value))+
# stat_density(aes(ymax = ..density..,  ymin = -..density..),
#                  fill = "grey50", colour = "grey50",
#                  geom = "ribbon", position = "identity")+facet_grid(. ~ FeatureIdx) +coord_flip()



##############################here plot for frequences thi is more relevant

##check do we have any surfaces available from the list
#topwellsurf<-read.csv(file ="F:/Alex/ICAM_validation/script/topowellfeat.csv")
mathits<-c(11,42,46,109,229,689,494,1050,1147,2114,1673,
           2150,443,864,1908,2075852,1106,330,1101, 917,991,1901,398,1108,2108,1018)
frhits<-c(11,1147,1050,494,1018,229,1673,2114,2150,864,443,389,46,42)
test<-c(1559,201,255,74,983,456)
#intersect(as.numeric(topwellsurf$FeatureIdx),c(bottomicam$FeatureIdx, topicam$FeatureIdx))
intersect(mathits,c(bottomicam$FeatureIdx, topicam$FeatureIdx))
intersect(frhits,c(bottomicam$FeatureIdx, topicam$FeatureIdx))
intersect(test,c(bottomicam$FeatureIdx, topicam$FeatureIdx))
#######
##save results
save(bottomicam,topicam,file="Hitscalculated from frequences.RDATA")
