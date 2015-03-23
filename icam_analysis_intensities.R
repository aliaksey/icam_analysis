
#hist(perobiint$MnCorrFctr)

#perobiint$ICAM_Intens_Median_Nor<-(perobiint$Cell_Intensity_MedianIntensity_Icam1amask-
#    perobiint$Image_Intensity_MinIntensity_Icam1Crop)

# #hist(perobiint$CorrFctr)
# boxplot(perobiint$ICAM_Intens_Mean_Nor~perobiint$Image_Metadata_array)
# #boxplot(perobiint$ICAM_Intens_Median_Nor~perobiint$Image_Metadata_array)
# 
# #performing kruskal wallis test
# perobiinttk<-perobiint[perobiint$Image_Metadata_array<9,]
# kruskal.test(perobiinttk$Cell_Intensity_MeanIntensity_Icam1amask~perobiinttk$Image_Metadata_array)
# kruskal.test(perobiinttk$Cell_Intensity_MedianIntensity_Icam1amask~perobiinttk$Image_Metadata_array)
# kruskal.test(perobiinttk$ICAM_Intens_Mean_Nor~perobiinttk$Image_Metadata_array)
# kruskal.test(perobiinttk$ICAM_Intens_Median_Nor~perobiinttk$Image_Metadata_array)
# 

# 
# plot(perobiint$ICAM_Intens_Mean_Nor,perobiint$ICAM_Intens_Median_Nor)
# 
# perobindallfnt$Cell_Intensity_MeanIntensity_Icam1amask<-perobindallfnt$Cell_Intensity_MeanIntensity_Icam1amask*255
# perobindallfnt$Cell_Intensity_MeanIntensity_Icam1amask<-perobindallfnt$Cell_Intensity_MeanIntensity_Icam1amask-
#   min(perobindallfnt$Cell_Intensity_MeanIntensity_Icam1amask)
# 
# perobindallfnp$Cell_Intensity_MeanIntensity_Icam1amask<-perobindallfnp$Cell_Intensity_MeanIntensity_Icam1amask*255
# perobindallfnp$Cell_Intensity_MeanIntensity_Icam1amask<-perobindallfnp$Cell_Intensity_MeanIntensity_Icam1amask-
#   min(perobindallfnp$Cell_Intensity_MeanIntensity_Icam1amask)
# 
# perobindallfnn$Cell_Intensity_MeanIntensity_Icam1amask<-perobindallfnn$Cell_Intensity_MeanIntensity_Icam1amask*255
# perobindallfnn$Cell_Intensity_MeanIntensity_Icam1amask<-perobindallfnn$Cell_Intensity_MeanIntensity_Icam1amask-
#   min(perobindallfnn$Cell_Intensity_MeanIntensity_Icam1amask)
# 
# library(ggplot2)
# #perobindallfn[perobindallf$Cell_Intensity_MeanIntensity_Icam1amask<0,"Cell_Intensity_MeanIntensity_Icam1amask"]=0
# #chekking different intensity descriptors:
# perobindall$Cell_Intensity_IntegratedIntensity_Actin1amask
# 
# plottest<-perobindall[,c("Cell_Intensity_MeanIntensity_Icam1amask",
#                          "Cell_AreaShape_Area",
#                          "Cell_Intensity_IntegratedIntensity_Icam1amask")]
# plottest$Calculated<-plottest$Cell_Intensity_IntegratedIntensity_Icam1amask/plottest$Cell_AreaShape_Area
# plot(plottest$Cell_Int)
# ggplot(plottest,aes(x=Cell_Intensity_MeanIntensity_Icam1amask,
#                     y=Calculated)) + geom_point() + geom_density2d()
# 
# ggplot(perobindallt,aes(x=Cell_Intensity_IntegratedIntensity_Icam1amask,
#                         y=Cell_Intensity_IntegratedIntensity_Actin1amask)) + geom_point() + geom_density2d()
# 
# ggplot(perobindallt,aes(x=Cell_Intensity_MeanIntensity_Icam1amask,
#                         y=Cell_Intensity_MeanIntensity_Actin1amask)) + geom_point() + geom_density2d()
# 
# ggplot(perobindallt,aes(x=Cell_Intensity_MeanIntensity_Icam1amask,
#                         y=Cell_Intensity_MeanIntensity_Actin1amask)) + geom_point() + 
#   geom_density2d() +ylim(0.25,0.45)+xlim(0.25,0.5)
# 
# ggplot(perobindallt,aes(x=Cell_Intensity_MeanIntensity_Icam1amask,
#                         y=Cell_Intensity_MeanIntensity_Icam1amask)) + geom_point() + 
#   geom_density2d() +ylim(0.25,0.45)+xlim(0.25,0.5)
# 
# 
# 
# 
# ggplot(perobindallt,aes(x=Cell_Intensity_MeanIntensity_Icam1amask,
#                         y=Cell_Intensity_MeanIntensity_Icam1amask)) + geom_point() + geom_density2d()
# 
# perobindallt$
#   
_________________________________________________________________________________
# for individual cells
#forplots<-merge(rank_tests, perobindallf[,c("ImageNumber","InternalIdx","Image_Metadata_array","Cell_Intensity_MeanIntensity_Icam1amask",
"Cell_Intensity_MeanIntensity_Actin1amask")], by="InternalIdx", sort=F)
#forplotscr<-merge(rank_testscrop, perobindallf[,c("ImageNumber","InternalIdx","Image_Metadata_array","Cell_Intensity_MeanIntensity_Icam1amask",
"Cell_Intensity_MeanIntensity_Actin1amask")], by="InternalIdx", sort=F)
#forplotsin<-merge(rank_testsinter, perobindallf[,c("ImageNumber","InternalIdx","Image_Metadata_array","Cell_Intensity_MeanIntensity_Icam1amask",
"Cell_Intensity_MeanIntensity_Actin1amask")], by="InternalIdx", sort=F)
##############################################
#############################################
###############################################
#calculating hits based on intensities
_______________________________________________________________________________
library(plyr) #at the beginning I used o.3 as trimmed mean coef 
icamintim<-ddply(perobindall.cort,"ImageNumber", summarise,  
                 FeatureIdx=mean(FeatureIdx),
                 IcamIntensityMeanmed=mean(ICAM_Intens_Median_corr,trim=0.2),
                 IcamIntensityMedmed=median(ICAM_Intens_Median_corr),
                 ActinIntensityMeanmed=mean(Actin_Intens_Median_corr,trim=0.2),
                 ActinIntensityMedmed=median(Actin_Intens_Median_corr),
                 IcamIntensityMeanmean=mean(Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr,trim=0.2),
                 IcamIntensityMedmmean=median(Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr),
                 ActinIntensityMeanmean=mean(Cell_Intensity_MeanIntensity_Actin1amask_Nor_corr,trim=0.2),
                 ActinIntensityMedmean=median(Cell_Intensity_MeanIntensity_Actin1amask_Nor_corr),
                 Array=mean(Image_Metadata_array))
# icamintft<-na.omit(ddply(icamintim,"FeatureIdx", summarise, 
#                              IcamIntensityMtf=mean(),

))
icamintft<-aggregate(icamintim[,c("FeatureIdx", "IcamIntensityMeanmed", 
                                  "IcamIntensityMedmed", "ActinIntensityMeanmed",
                                  "ActinIntensityMedmed", "IcamIntensityMeanmean",
                                  "IcamIntensityMedmmean", "ActinIntensityMeanmean",
                                  "ActinIntensityMedmean")], by=list(icamintim$FeatureIdx),
                     FUN=mean, na.action = na.omit)

plot(icamintft[order(icamintft$IcamIntensityMeanmean),"IcamIntensityMeanmean"])

pairs(icamintft[,3:10])
pairs(icamintft[,3:10])
pairs(icamintft[,c(3,4,7,8)])
pairs(icamintft[,c(5,6,9,10)])
pairs(icamintft[,c(5,6,3,4)])

# running statistical test to find surfaces from hits that are statistically different.
for (i in 1:length(unique(icamintft[,"FeatureIdx"]))){ 
  temp <- icamintim[icamintim$FeatureIdx==icamintft[i,"FeatureIdx"],]
  pva<-wilcox.test(temp[,"IcamIntensityMeanmean"],icamintim[,"IcamIntensityMeanmean"] )
  #pvak<-kruskal.test(Cell_Intensity_MeanIntensity_Icam1amask~ImageNumber,data=temp2)
  icamintft[i,"p.value"]<-as.numeric(pva["p.value"])
}
_________________________________________________________________________________
#ordering data based on median and p value 
icamintftp<-icamintft[icamintft$p.value<0.05,]
#on mean intensity
top.hit<-icamintftp[order(-icamintftp$IcamIntensityMeanmean),][1:20,]
bottom.hit<-icamintftp[order(icamintftp$IcamIntensityMeanmean),][1:20,]
top.hit.m<-merge(top.hit,icamintim, by= "FeatureIdx", sort=F, all=F)
bottom.hit.m<-merge(bottom.hit,icamintim, by= "FeatureIdx", sort=F, all=F)
forplotcoll<-rbind(bottom.hit.m[,c("FeatureIdx","IcamIntensityMeanmean.y","ActinIntensityMeanmean.y")], 
                   top.hit.m[length(top.hit.m[,1]):1,c("FeatureIdx","IcamIntensityMeanmean.y",
                                                       "ActinIntensityMeanmean.y")])
#plot all evrything without dots by mean
forplottransf <- transform(forplotcoll, FeatureIdx = factor(FeatureIdx, 
                                                            levels = as.character(forplotcoll$FeatureIdx)))
forplottransfmelt<-melt(forplottransf, measure.vars =c("IcamIntensityMeanmean.y",
                                                       "ActinIntensityMeanmean.y"))
#forplottransfmelt<-melt(forplottransf, measure.vars ="IcamIntensityMt")
p91<- ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))
p91+geom_boxplot()+theme(legend.position="none")+ylim(-0.5,10)

# #on integrated intensity
# collobfeatfpint<-collobfeat[collobfeat$int.p.value<0.05,]
# 
# top.hit.int<-collobfeatfpint[order(-collobfeatfpint$IcamINtegratedIntMean,-collobfeatfpint$CellNumber),][1:20,]
# bottom.hit.int<-collobfeatfpint[order(collobfeatfpint$IcamINtegratedIntMean,-collobfeatfpint$CellNumber),][1:20,]
# top.hit.m.int<-merge(top.hit.int,collob, by= "FeatureIdx", sort=F, all=F)
# bottom.hit.m.int<-merge(bottom.hit.int,collob, by= "FeatureIdx", sort=F, all=F)
# forplotcollint<-rbind(bottom.hit.m.int[,c("FeatureIdx","IcamINtegratedInt")], 
#                       top.hit.m.int[length(top.hit.m.int[,1]):1,c("FeatureIdx","IcamINtegratedInt")])
# _______________________________________________________________
library(ggplot2)
require (plyr)
library(reshape2)

# constructing ranking based on this findingsmedians
rankmedian<-collobfeat[order(collobfeat$IcamIntensityMtf),c("FeatureIdx","IcamIntensityMtf", "CellNumber")]
ranktrans <- transform(rankmedian, FeatureIdx = factor(FeatureIdx, 
                                                       levels = as.character(rankmedian$FeatureIdx)))
ranktransmelt<-melt(ranktrans, measure.vars = "IcamIntensityMtf")
plot(collob$CellNumber~collob$IcamIntensityMt)
median(collob$CellNumber)

plot(rankmedian$IcamIntensityMtf, ylim=range(2:20), log='y',xlab="Rank Index", 
     ylab="Trimmed Mean Icam Intensity")
abline(h=mean(perobindallfnn$Cell_Intensity_MeanIntensity_Icam1amask, trim=0.2))
abline(h=mean(perobindallfnp$Cell_Intensity_MeanIntensity_Icam1amask, trim=0.2))
abline(h=mean(perobindallfnt$Cell_Intensity_MeanIntensity_Icam1amask, trim=0.3), col="red")
#for median per unit
#forplotcoll<-merge(rank_tests, collob,by="InternalIdx", sort=F)
#for median per unit cropped rank
#forplotcollcr<-merge(rank_testscrop, collob,by="InternalIdx", sort=F)
#forplotcollin<-merge(rank_testsinter, collob,by="InternalIdx", sort=F) 
___________________________________________________________________-
  # plotting alll surface medians
  _____________________________________________________________________________________________
# plotting controls
topot<-as.data.frame(perobindallfnt$Cell_Intensity_MeanIntensity_Icam1amask)
colnames(topot)<-"IcamTrMeanIntensity"
topot$sample <- 'Topochip'
negat<-as.data.frame(perobindallfnn$Cell_Intensity_MeanIntensity_Icam1amask)
colnames(negat)<-"IcamTrMeanIntensity"
negat$sample <- 'Negative Control'
posit<-as.data.frame(perobindallfnp$Cell_Intensity_MeanIntensity_Icam1amask)
colnames(posit)<-"IcamTrMeanIntensity"
posit$sample <- 'Positive Control'
intencontnp<-rbind(negat, topot, posit)
intencontnp2<-rbind(negat, posit)

p11<-ggplot(intencontnp2, aes(sample, y = IcamTrMeanIntensity, fill = sample)) 
p11+geom_boxplot()

p22<-ggplot(intencontnp2, aes(IcamTrMeanIntensity, fill = sample)) 
p22+geom_density(alpha = 0.2)
p22+geom_histogram(alpha = 0.7, aes(y = ..density..), position = 'identity')

___________________________________________________________________________________
# calculating hits from ratio
#count number of postitive cells
ratiorankw<-perobindallfnt[perobindallfnt$Cell_Intensity_MeanIntensity_Icam1amask > 10,
                           "FeatureIdx"]
as<-table(ratiorankw)
names(which(table(ratiorankw) == max(table(ratiorankw))))

library(plyr)
icam_ratiotf<-ddply(perobindallfnt,"ImageNumber", summarise, 
                    IcamPositive=sum(Cell_Intensity_MeanIntensity_Icam1amask > 10),
                    Total=sum(Cell_Intensity_MeanIntensity_Icam1amask > 0),
                    TrMeanInt=mean(Cell_Intensity_MeanIntensity_Icam1amask, trimm=0.3),
                    MedianIcam=median(Cell_Intensity_MeanIntensity_Icam1amask),
                    FeatureIdx=mean(FeatureIdx))
#plot(icam_ratiotf$MedianIcam~icam_ratio$TrMeanInt)

icam_ratiotf$Ratio<-icam_ratiotf$IcamPositive/icam_ratiotf$Total

#collapse to feaature nuber
icam_ratiotfff<-ddply(na.omit(icam_ratiotf),"FeatureIdx", summarise, 
                      RatioTrMean=mean(Ratio),
                      RatioMedian=median(Ratio),
                      IntensTrMean=mean(TrMeanInt),
                      IntensTrMad=mad(TrMeanInt),
                      RatioSd=sd(Ratio),
                      RatioMad=mad(Ratio))

rankratio<-icam_ratiotfff[order(icam_ratiotfff$RatioTrMean),c("FeatureIdx","RatioTrMean")]
# count total number of 
head(rankratio)
tail(rankratio)
######## calculating chi squred statistics
hit<-perobindallfnt[perobindallfnt$FeatureIdx==434,"Cell_Intensity_MeanIntensity_Icam1amask"]
negg<-perobindallfnn$Cell_Intensity_MeanIntensity_Icam1amask
tr<-10
obs<-c(sum(hit>tr),sum(hit<tr))
obsratio<-c(sum(hit>tr)/length(hit),sum(hit<tr)/length(hit))
obsratio
obs
exp<-c(sum(negg>tr)/length(negg),sum(negg<tr)/length(negg))
exp
#performing statistics
#chisquared
#chisq.test(ct)
chisq.test(obs,p=exp)



#plot(icam_ratiof$RatioMean~icam_ratiof$IntensTrMean)
#plot(icam_ratiof$IntensTrMad~icam_ratiof$IntensTrMean)


library(ggplot2)
require (plyr)
library(reshape2)
#transforming and melting data for croped
#forplottransf <- transform(forplotscr, InternalIdx = factor(InternalIdx, levels = as.character(forplotscr$InternalIdx)))
#forplottransfmelt<-melt(forplottransf, measure.vars = c("Cell_Intensity_MeanIntensity_Icam1amask","Cell_Intensity_MeanIntensity_Actin1amask"))
______________________________________________________________________
#plot all evrything without dots by mean
forplottransf <- transform(forplotcoll, FeatureIdx = factor(FeatureIdx, 
                                                            levels = as.character(forplotcoll$FeatureIdx)))
#forplottransfmelt<-melt(forplottransf, measure.vars =c("IcamIntensityMt","ActinIntensityMt"))
forplottransfmelt<-melt(forplottransf, measure.vars ="IcamIntensityMt")
p91<- ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))
p91+geom_boxplot()+theme(legend.position="none")+ylim(-0.5,10)

#with integrated intensity
forplottransfint<- transform(forplotcollint, FeatureIdx = factor(FeatureIdx, 
                                                                 levels = as.character(forplotcollint$FeatureIdx)))
forplottransfmeltint<-melt(forplottransfint, measure.vars = "IcamINtegratedInt") #,"Cell_Intensity_MeanIntensity_Actin1amask"
p9991<- ggplot(forplottransfmeltint, aes(FeatureIdx, y = value, fill=variable))
p9991+geom_boxplot()+theme(legend.position="none")
#theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1)) 
