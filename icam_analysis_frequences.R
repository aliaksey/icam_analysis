load("Icam_cells_after_correction.RDATA")

______________________________________________________________________________________
#finding hits based on frequencyies
#find treshhold factor
plot(perobiintc$ICAM_Intens_Mean_corr,perobiintc$Cell_AreaShape_Area, )
colicam=as.factor(perobiintc[perobiintc$Image_Metadata_array>8,"Image_Metadata_array"])
plot(perobiintc[perobiintc$Image_Metadata_array>8,"ICAM_Intens_Mean_corr"],
     perobiintc[perobiintc$Image_Metadata_array>8,"Cell_AreaShape_Area"], 
     col=colicam)
library(ggplot2)
library(reshape2)
forggplot<-perobiintc[perobiintc$Image_Metadata_array>8,c("ICAM_Intens_Mean_corr",
                                                          "Image_Metadata_array")]
p1<-ggplot(forggplot, aes(ICAM_Intens_Mean_corr, fill = as.factor(forggplot$Image_Metadata_array))) 
p1+geom_density(alpha = 0.2)

p2<-ggplot(perobiintc, aes(ICAM_Intens_Mean_corr, fill = as.factor(perobiintc$Image_Metadata_array))) 
p2+geom_density(alpha = 0.2)


hist(perobiintc[perobiintc$Image_Metadata_array>8,"ICAM_Intens_Mean_corr"])
dd<-density(perobiintc[perobiintc$Image_Metadata_array>8,"ICAM_Intens_Mean_corr"])
plot(dd)

boxplot(perobiintc[perobiintc$ICAM_Intens_Mean_corr > 10,
                   "ICAM_Intens_Mean_corr"]~
          perobiintc[perobiintc$ICAM_Intens_Mean_corr > 10, "Image_Metadata_array"])
###############################################################
#################################################################
#count number of postitive cells

perobiintct<-perobiintc[perobiintc$Image_Metadata_array<9,]
perobiintcp<-perobiintc[perobiintc$Image_Metadata_array==9,]
perobiintcn<-perobiintc[perobiintc$Image_Metadata_array==10,]
#making some plots
plot(perobiintcn$Actin_Intens_Mean_corr,
     perobiintcn$ICAM_Intens_Mean_corr)

library(plyr)
ratiorankcorr<-ddply(perobiintct,"ImageNumber", summarise, 
                     IcamPositive=sum(ICAM_Intens_Mean_corr > 10),
                     Total=sum(ICAM_Intens_Mean_corr > 0),
                     TrMeanIcam=mean(ICAM_Intens_Mean_corr, trimm=0.3),
                     FeatureIdx=mean(FeatureIdx))

ratiorankcorrn<-ddply(perobiintcn,"ImageNumber", summarise, 
                      IcamPositive=sum(ICAM_Intens_Mean_corr > 10),
                      Total=sum(ICAM_Intens_Mean_corr > 0),
                      TrMeanIcam=mean(ICAM_Intens_Mean_corr, trimm=0.3),
                      FeatureIdx=mean(FeatureIdx))
#plot(ratiorankcorr$IcamPositive~ratiorankcorr$FeatureIdx)

ratiorankcorr$Ratio<-ratiorankcorr$IcamPositive/ratiorankcorr$Total

plot(ratiorankcorr$Ratio,ratiorankcorr$TrMeanIcam)

#collapse to feaature nuber
ratiorankcorrfff<-ddply(na.omit(ratiorankcorr),"FeatureIdx", summarise, 
                        RatioTrMean=mean(Ratio),
                        RatioMedian=median(Ratio),
                        IntensTrMean=mean(TrMeanIcam),
                        IntensTrMad=mad(TrMeanIcam),
                        RatioSd=sd(Ratio),
                        RatioMad=mad(Ratio))

ratiorankcorrfffs<-ratiorankcorrfff[order(ratiorankcorrfff$RatioTrMean),]

plot(ratiorankcorrfffs$RatioTrMean)
#calculating statistics
for (i in 1:length(unique(ratiorankcorrfff[,"FeatureIdx"]))){ 
  temp <- ratiorankcorr[ratiorankcorr$FeatureIdx==ratiorankcorrfff[i,"FeatureIdx"],]
  #temp2<-perobindallfnt[perobindallfnt$FeatureIdx==collobfeat[i,"FeatureIdx"],]
  obs<-c(sum(temp[,"IcamPositive"]),(sum(temp[,"Total"])-sum(temp[,"IcamPositive"])))
  exp<-c(sum(perobiintcn$ICAM_Intens_Mean_corr >10)/length(perobiintcn$ICAM_Intens_Mean_corr),
         sum(perobiintcn$ICAM_Intens_Mean_corr<10)/length(perobiintcn$ICAM_Intens_Mean_corr))
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

rankration<-ratiorankcorrfff[order(ratiorankcorrfff$RatioTrMean,
                                   ratiorankcorrfff$RatioSd),c("FeatureIdx","RatioTrMean")]
# count total number of 
#checking do we have any surfaces prepared
topicam<-rankration[c(2157:2177),"FeatureIdx"]
bottomicam<-rankration[c(1:20),"FeatureIdx"]
# plotting intensities for actin and icam
library(ggplot2)
dataicamav<-aggregate(perobiintct[,c("FeatureIdx", "ImageNumber", "Actin_Intens_Mean_corr", 
                                     "ICAM_Intens_Mean_corr")], by=list(perobiintct$ImageNumber),
                      FUN=function(x) mean(x, trim=0.2, na.action = na.omit))

topicamm<-dataicamav[dataicamav$FeatureIdx%in%topicam,]
bottomicamm<-dataicamav[dataicamav$FeatureIdx%in%bottomicam,]
forplotcoll<-rbind(topicamm,bottomicamm)
forplotcoll<-forplotcoll[order(forplotcoll$ICAM_Intens_Mean_corr),]

plot(forplotcoll$Actin_Intens_Mean_corr,forplotcoll$ICAM_Intens_Mean_corr)
cor.test(forplotcoll$Actin_Intens_Mean_corr,forplotcoll$ICAM_Intens_Mean_corr)

#plot all evrything without dots by mean
forplottransf <- transform(forplotcoll[,c("FeatureIdx", "Actin_Intens_Mean_corr","ICAM_Intens_Mean_corr") ], FeatureIdx = factor(FeatureIdx, 
                                                                                                                                 levels = as.character(forplotcoll$FeatureIdx)))
library(reshape2)
forplottransfmelt<-melt(forplottransf, measure.vars =c("ICAM_Intens_Mean_corr",
                                                       "Actin_Intens_Mean_corr"))
#forplottransfmelt<-melt(forplottransf, measure.vars ="IcamIntensityMt")
p91<- ggplot(forplottransfmelt, aes(FeatureIdx, y = value, fill=variable))
p91+geom_boxplot()+theme(legend.position="none")+ylim(-0.5,10)




topwellsurf<-read.csv(file ="F:/Alex/ICAM_validation/script/topowellfeat.csv")
mathits<-c(11,42,46,109,229,689,494,1050,1147,2114,1673,
           2150,443,864,1908,2075852,1106,330,1101, 917,991,1901,398,1108,2108,1018)
frhits<-c(11,1147,1050,494,1018,229,1673,2114,2150,864,443,389,46,42)

intersect(as.numeric(topwellsurf$FeatureIdx),c(bottomicam, topicam))
#######


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
icamintim<-ddply(perobiintct,"ImageNumber", summarise,  
                 FeatureIdx=mean(FeatureIdx),
                 IcamIntensityMeanmed=mean(ICAM_Intens_Median_corr,trim=0.2),
                 IcamIntensityMedmed=median(ICAM_Intens_Median_corr),
                 ActinIntensityMeanmed=mean(Actin_Intens_Median_corr,trim=0.2),
                 ActinIntensityMedmed=median(Actin_Intens_Median_corr),
                 IcamIntensityMeanmean=mean(ICAM_Intens_Mean_corr,trim=0.2),
                 IcamIntensityMedmmean=median(ICAM_Intens_Mean_corr),
                 ActinIntensityMeanmean=mean(Actin_Intens_Mean_corr,trim=0.2),
                 ActinIntensityMedmean=median(Actin_Intens_Mean_corr),
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
