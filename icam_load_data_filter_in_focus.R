rm(list=ls())
load("data/icamimageandobjectindex.RData")
#selecting only variables that dealing with focus

perimindallf<-perimindall[,c("ImageNumber","FeatureIdx","InternalIdx","Image_Metadata_array",
              colnames(perimindall)[grepl("Focus",colnames(perimindall))],
              colnames(perimindall)[grepl("focus",colnames(perimindall))],
              colnames(perimindall)[grepl("LogLogSlope",colnames(perimindall))])]

#plotting log logslopes of all images
library (ggplot2)

perobindallfplot <- as.data.frame(perimindallf)
#show were is blamk topo is
highlight.topo <- 2177
perobindallfplot$highlight <- ifelse(perobindallfplot$FeatureIdx == highlight.topo, "highlight", "normal")
textdf <- perobindallfplot[perobindallfplot$FeatureIdx == highlight.topo, ]
mycolours <- c("highlight" = "red", "normal" = "grey50")

p11<-ggplot(data = perobindallfplot, aes(x = Image_ImageQuality_PowerLogLogSlope_Actin1Crop,
                                         y = Image_ImageQuality_PowerLogLogSlope_Dapi1Crop)) 
p11+ geom_point(alpha = 0.4) +geom_density2d()+xlab("PowerLogLogSlope_Phaloidin")+
  ylab("PowerLogLogSlope_DAPI")+ geom_segment(aes(x = -1.25, y = -3, xend = -1.25, yend = 0, colour="red"))+
  geom_segment(aes(y = -1.3, x = -3, yend = -1.3, xend = 0, colour="red"))

focim<-perimindallf[perimindallf$Image_ImageQuality_PowerLogLogSlope_Actin1Crop<(-1.25)&
                      perimindallf$Image_ImageQuality_PowerLogLogSlope_Dapi1Crop<(-1.3),"ImageNumber"]
length(focim)
#plotting results
perobindallfplotinfocus<-perobindallfplot[perobindallfplot$ImageNumber%in%focim,]
ggplot(data = perobindallfplotinfocus, aes(x = Image_ImageQuality_PowerLogLogSlope_Actin1Crop,
                                         y = Image_ImageQuality_PowerLogLogSlope_Dapi1Crop))+ 
geom_point(alpha = 0.4) +geom_density2d()+xlab("PowerLogLogSlope_Phaloidin")+
  ylab("PowerLogLogSlope_DAPI")

# finding images with over exposion (in some cases max and min intensity in ICAm is eqaual
#that image is not suitable for further nanalysis

perimindallftt<-perimindall[perimindall$Image_Metadata_array>8,]
perimindalltt<-perimindall[perimindall$Image_Metadata_array<9,]
p12<-ggplot(data = perimindalltt, aes(x = Image_ImageQuality_MinIntensity_Icam1Crop,
                                      y = Image_ImageQuality_MedianIntensity_Icam1Crop)) 
p12+ geom_point(alpha = 0.4) +geom_density2d()

corrimt<-perimindall[perimindall$Image_ImageQuality_MinIntensity_Icam1Crop<0.28&
                       perimindall$Image_ImageQuality_MedianIntensity_Icam1Crop<0.325&
                       perimindall$Image_Metadata_array<9,"ImageNumber"]
corrimc<-perimindall[perimindall$Image_Metadata_array>8,"ImageNumber"]

corrim<-c(corrimt,corrimc)
length(unique(corrim))

save(focim, file="focim.RData")
save(corrim, file="corrim.RData")