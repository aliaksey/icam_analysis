  rm(list=ls())
  load("D:/projects/icam_project/icam_analysis/data/icamimageandobjectindex.RData")
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
  
  ggplot(data = perobindallfplot, aes(x = Image_ImageQuality_PowerLogLogSlope_Actin1Crop,
                                           y = Image_ImageQuality_PowerLogLogSlope_Dapi1Crop))+
  geom_point(alpha = 0.4) +geom_density2d()+xlab("PowerLogLogSlope_Phaloidin")+
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
  #__________________________________________________________________________________________________________________
  
  # finding images with over exposion (in some cases max and min intensity in ICAm is close
  #that image is not suitable for further nanalysis
  
  ##selecting only images that passed focus filter
  
  perimindall.infoc<-perimindall[perimindall$ImageNumber%in%focim,]
  
  perimindalltt<-perimindall.infoc
  
  imageicamintensity<-colnames(perimindalltt)[grepl("Icam1Crop",colnames(perimindalltt))&
                          grepl("ImageQuality",colnames(perimindalltt))&
                            grepl("Intensity",colnames(perimindalltt))]
  library(GGally)
  for.pairs.plot<-perimindalltt[,c(imageicamintensity,"Image_Metadata_array")]
  colnames(for.pairs.plot)<-gsub("Image_ImageQuality_","",colnames(for.pairs.plot))
  colnames(for.pairs.plot)<-gsub("_Icam1Crop","",colnames(for.pairs.plot))
  for.pairs.plot$Image_Metadata_array<-as.factor(for.pairs.plot$Image_Metadata_array)
  library(gplots)
  palette(rev(rich.colors(10)))
  ggpairs(for.pairs.plot,columns=1:7,colour="Image_Metadata_array")
  palette("default")
  
  ##without controls
  perimindalltt.t<-perimindall.infoc[perimindall.infoc$Image_Metadata_array>8,]
  ##make some boxplot calculations
  
  imageicamintensity<-colnames(perimindalltt)[grepl("Icam1Crop",colnames(perimindalltt))&
                                                grepl("ImageQuality",colnames(perimindalltt))&
                                                grepl("Intensity",colnames(perimindalltt))]
  for.pairs.plot.t<-perimindalltt.t[,c(imageicamintensity,"Image_Metadata_array")]
  colnames(for.pairs.plot.t)<-gsub("Image_ImageQuality_","",colnames(for.pairs.plot.t))
  colnames(for.pairs.plot.t)<-gsub("_Icam1Crop","",colnames(for.pairs.plot.t))
  for.pairs.plot.t$Image_Metadata_array<-as.factor(for.pairs.plot.t$Image_Metadata_array)
  
  ####boxplot of intensities
  
  library(ggplot2)
  library(plyr)
  library(reshape)
  perimindalltt.t.melt<-melt(for.pairs.plot.t,measure.vars = 1:7)
  
  ggplot(perimindalltt.t.melt, aes(fill=Image_Metadata_array, y=value, x=variable))+
    geom_boxplot()+
    facet_grid(variable~., scales="free")+
    theme(axis.text.x=element_blank())
  
  boxplot(perimindall$Image_Crop_AreaRetainedAfterCropping_Icam1Crop~perimindall$Image_Metadata_array)
  ##matrixplot
  library(GGally)
  library(gplots)
  palette(rev(rich.colors(2)))
  ggpairs(for.pairs.plot.t,columns=1:7,colour="Image_Metadata_array")
  palette("default")
  
  ##only controls
  perimindalltt.t<-perimindall.infoc[perimindall.infoc$Image_Metadata_array<9,]
  ##make some boxplot calculations
  
  imageicamintensity<-colnames(perimindalltt)[grepl("Icam1Crop",colnames(perimindalltt))&
                                                grepl("ImageQuality",colnames(perimindalltt))&
                                                grepl("Intensity",colnames(perimindalltt))]
  for.pairs.plot.t<-perimindalltt.t[,c(imageicamintensity,"Image_Metadata_array")]
  colnames(for.pairs.plot.t)<-gsub("Image_ImageQuality_","",colnames(for.pairs.plot.t))
  colnames(for.pairs.plot.t)<-gsub("_Icam1Crop","",colnames(for.pairs.plot.t))
  for.pairs.plot.t$Image_Metadata_array<-as.factor(for.pairs.plot.t$Image_Metadata_array)
  
  ####boxplot of intensities
  
  library(ggplot2)
  library(plyr)
  library(reshape)
  perimindalltt.t.melt<-melt(for.pairs.plot.t,measure.vars = 1:7)
  
  ggplot(perimindalltt.t.melt, aes(fill=Image_Metadata_array, y=value, x=variable))+
    geom_boxplot()+
    facet_grid(variable~., scales="free")+
    theme(axis.text.x=element_blank())
  
  boxplot(perimindall$Image_Crop_AreaRetainedAfterCropping_Icam1Crop~perimindall$Image_Metadata_array)
  ##matrixplot
  library(GGally)
  library(gplots)
  palette(rev(rich.colors(8)))
  ggpairs(for.pairs.plot.t,columns=1:7,colour="Image_Metadata_array")
  palette("default")
### eliminating images based on minintensity value
corrimt<-perimindall.infoc[perimindall.infoc$Image_ImageQuality_MinIntensity_Icam1Crop<0.275&
                               perimindall.infoc$Image_ImageQuality_MedianIntensity_Icam1Crop<0.3&
                               perimindall.infoc$Image_Metadata_array<9,"ImageNumber"]
  corrim<-c(corrimt,perimindall.infoc[perimindall.infoc$Image_Metadata_array>8,"ImageNumber"])
  
#make plots with new filtering

  perimindall.cor<-perimindall.infoc[perimindall.infoc$ImageNumber%in%corrim,]
  perimindalltt.c<-perimindall.cor
  imageicamintensity<-colnames(perimindalltt)[grepl("Icam1Crop",colnames(perimindalltt))&
                                                grepl("ImageQuality",colnames(perimindalltt))&
                                                grepl("Intensity",colnames(perimindalltt))]
  library(GGally)
  for.pairs.plot.c<-perimindalltt.c[,c(imageicamintensity,"Image_Metadata_array")]
  colnames(for.pairs.plot.c)<-gsub("Image_ImageQuality_","",colnames(for.pairs.plot.c))
  colnames(for.pairs.plot.c)<-gsub("_Icam1Crop","",colnames(for.pairs.plot.c))
  for.pairs.plot.c$Image_Metadata_array<-as.factor(for.pairs.plot.c$Image_Metadata_array)
  library(gplots)
  palette(rev(rich.colors(10)))
  ggpairs(for.pairs.plot.c,columns=1:7,colour="Image_Metadata_array")
  palette("default")
  
  perimindalltt.c.melt<-melt(for.pairs.plot.c,measure.vars = 1:7)
  
  ggplot(perimindalltt.c.melt, aes(fill=Image_Metadata_array, y=value, x=variable))+
    geom_boxplot()+
    facet_grid(variable~., scales="free")+
    theme(axis.text.x=element_blank())
##calculate statistics for removed images:
  print("removed after focus and intensity filter")
  table(perimindall[!perimindall$ImageNumber%in%perimindall.cor$ImageNumber,"Image_Metadata_array"])
  print("removed after focus  filter")
  table(perimindall[!perimindall$ImageNumber%in%perimindall.infoc$ImageNumber,"Image_Metadata_array"])
  
    
p12<-ggplot(data = perimindall.cor, aes(x = Image_ImageQuality_MinIntensity_Icam1Crop,
                                      y = Image_ImageQuality_MedianIntensity_Icam1Crop)) 
p12+ geom_point(alpha = 0.4) +geom_density2d()

length(unique(corrim))==length(corrim)
# filter images in objects
  perobindall.cor<-perobindall[perobindall$ImageNumber%in%corrim,]
  perobindall.infoc<-perobindall[perobindall$ImageNumber%in%focim,]
  
save(perimindall.infoc, perobindall.infoc, file="focim.RData")
save(perimindall.cor,perobindall.cor, file="corrim.RData")