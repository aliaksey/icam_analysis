rm(list=ls())
load("cell_objects_after_area_perimeter.RDATA")
load("corrim.RData")

##selecting subpopulation of surfaces
perobindallt<-na.omit(perobindall.filter[perobindall.filter$Image_Metadata_array<9,])
perobindallp<-na.omit(perobindall.filter[perobindall.filter$Image_Metadata_array==9,])
perobindalln<-na.omit(perobindall.filter[perobindall.filter$Image_Metadata_array==10,])

# perobindallinfms<-perobindall.filter
# perimindall<-perimindall.cor
# 

#normolize cell intensity based on min intensity

#boxplot(perobiint$Image_Intensity_MinIntensity_Icam1Crop~perobiint$Image_Metadata_array)
#boxplot(perobiint$Cell_Intensity_MeanIntensity_Icam1amask~perobiint$Image_Metadata_array)
#boxplot(perobiint$Cell_Intensity_MedianIntensity_Icam1amask~perobiint$Image_Metadata_array)

#normalize data based on min intensities in the image 
perobiint<-na.omit(merge(perimindall.cor[,c("InternalIdx","FeatureIdx","UnitIdx","Col",
                  "Row","ImageNumber","Image_Metadata_array",colnames(perimindall.cor)[grepl("Image_Intensity_",
                                                                    colnames(perimindall.cor))|grepl(
                                                                      "Image_ImageQuality_",
                                                                    colnames(perimindall.cor))])],
                         perobindall.filter[,c("ImageNumber","Cell_AreaShape_Area","ObjectNumber", colnames( perobindall.filter)[grepl("_Intensity_",
                         colnames( perobindall.filter))])], by="ImageNumber", all.x=T))
##make sevearl checkups
library(ggplot2)
ggplot(perobiint[perobiint$Image_Metadata_array<9,],
       aes(x=Cell_Intensity_MedianIntensity_Icam1Crop,y=Image_Intensity_MinIntensity_Icam1Crop,
           colour=Image_Metadata_array))+geom_point()
ggplot(perobiint[perobiint$Image_Metadata_array<9,],
       aes(x=Cell_Intensity_MedianIntensity_Icam1amask,y=Image_Intensity_MinIntensity_Icam1Crop,
           colour=Image_Metadata_array))+geom_point()


#subtracting min image intensity from cell intensities 
#correct Icam intensities
perobiint$Cell_Intensity_MeanIntensity_Icam1amask_Nor<-(perobiint$Cell_Intensity_MeanIntensity_Icam1amask-
                                   perobiint$Image_Intensity_MinIntensity_Icam1Crop)
perobiint$Cell_Intensity_MedianIntensity_Icam1amask_Nor<-(perobiint$Cell_Intensity_MedianIntensity_Icam1amask-
                                     perobiint$Image_Intensity_MinIntensity_Icam1Crop)
perobiint$Cell_Intensity_IntegratedIntensity_Icam1amask_Nor<-(perobiint$Cell_Intensity_IntegratedIntensity_Icam1amask-
                                        (perobiint$Image_Intensity_MinIntensity_Icam1Crop*perobiint$Cell_AreaShape_Area))
#correct Actin intensiites
perobiint$Cell_Intensity_MeanIntensity_Actin1amask_Nor<-(perobiint$Cell_Intensity_MeanIntensity_Actin1amask-
                                    perobiint$Image_ImageQuality_MinIntensity_Actin1Crop)
perobiint$Cell_Intensity_MedianIntensity_Actin1amask_Nor<-(perobiint$Cell_Intensity_MedianIntensity_Actin1amask-
                                      perobiint$Image_ImageQuality_MinIntensity_Actin1Crop)
perobiint$Cell_Intensity_IntegratedIntensity_Actin1amask_Nor<-(perobiint$Cell_Intensity_IntegratedIntensity_Actin1amask-
                              (perobiint$Image_ImageQuality_MinIntensity_Actin1Crop*perobiint$Cell_AreaShape_Area))

##visualize changes
#Icam
#Median Intensity
ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_MedianIntensity_Icam1amask,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_MedianIntensity_Icam1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
#Mean intenssity
ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_MeanIntensity_Icam1amask,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_MeanIntensity_Icam1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
#Integrated Intensity
ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_IntegratedIntensity_Icam1amask,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_IntegratedIntensity_Icam1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
#Actin
#median intensity
ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_MedianIntensity_Actin1amask,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_MedianIntensity_Actin1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
##mean intensity
ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_MeanIntensity_Actin1amask,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_MeanIntensity_Actin1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
##integrated intensity
ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_IntegratedIntensity_Actin1amask,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiint,aes(factor(Image_Metadata_array), Cell_Intensity_IntegratedIntensity_Actin1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()



##medians for all chips should be the same
#calculating correction factor for each topochip individially
perobiintc<-c()
for(i in 1:10){
  if (i!=9){
    
    temp<-perobiint[perobiint$Image_Metadata_array==i,]
    # for ICAM correction
    #mean
    temp$corrfct<-(mean(c(median(perobiint[perobiint$Image_Metadata_array==1,"Cell_Intensity_MeanIntensity_Icam1amask_Nor"]),
                        median(perobiint[perobiint$Image_Metadata_array==2,"Cell_Intensity_MeanIntensity_Icam1amask_Nor"]),
                        median(perobiint[perobiint$Image_Metadata_array==3,"Cell_Intensity_MeanIntensity_Icam1amask_Nor"]))))/
      median(temp$Cell_Intensity_MeanIntensity_Icam1amask_Nor)
    
    temp$Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr<-temp$Cell_Intensity_MeanIntensity_Icam1amask_Nor*temp$corrfct*255
    
    #median
    temp$corrfctmd<- (mean(c(median(perobiint[perobiint$Image_Metadata_array==1,"Cell_Intensity_MedianIntensity_Icam1amask_Nor"]),
                             median(perobiint[perobiint$Image_Metadata_array==2,"Cell_Intensity_MedianIntensity_Icam1amask_Nor"]),
                             median(perobiint[perobiint$Image_Metadata_array==3,"Cell_Intensity_MedianIntensity_Icam1amask_Nor"]))))/
      median(temp$Cell_Intensity_MedianIntensity_Icam1amask_Nor)
    
    temp$Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr<-temp$Cell_Intensity_MedianIntensity_Icam1amask_Nor*temp$corrfctmd*255
    
    #Integrated
    temp$corrfctint<- (mean(c(median(perobiint[perobiint$Image_Metadata_array==1,"Cell_Intensity_IntegratedIntensity_Icam1amask_Nor"]),
                              median(perobiint[perobiint$Image_Metadata_array==2,"Cell_Intensity_IntegratedIntensity_Icam1amask_Nor"]),
                              median(perobiint[perobiint$Image_Metadata_array==3,"Cell_Intensity_IntegratedIntensity_Icam1amask_Nor"]))))/
      median(temp$Cell_Intensity_IntegratedIntensity_Icam1amask_Nor)
    
    temp$Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr<-temp$Cell_Intensity_IntegratedIntensity_Icam1amask_Nor*temp$corrfctint*255
    
    # for Actin correction
    #mean
    temp$corrfcta<-(mean(c(median(perobiint[perobiint$Image_Metadata_array==1,"Cell_Intensity_MeanIntensity_Actin1amask_Nor"]),
                          median(perobiint[perobiint$Image_Metadata_array==2,"Cell_Intensity_MeanIntensity_Actin1amask_Nor"]),
                          median(perobiint[perobiint$Image_Metadata_array==3,"Cell_Intensity_MeanIntensity_Actin1amask_Nor"]))))/
      median(temp$Cell_Intensity_MeanIntensity_Actin1amask_Nor)
    
    temp$Cell_Intensity_MeanIntensity_Actin1amask_Nor_corr<-temp$Cell_Intensity_MeanIntensity_Actin1amask_Nor*temp$corrfcta*255
    
    #median
    temp$corrfctmda<- (mean(c(median(perobiint[perobiint$Image_Metadata_array==1,"Cell_Intensity_MedianIntensity_Actin1amask_Nor"]),
                             median(perobiint[perobiint$Image_Metadata_array==2,"Cell_Intensity_MedianIntensity_Actin1amask_Nor"]),
                             median(perobiint[perobiint$Image_Metadata_array==3,"Cell_Intensity_MedianIntensity_Actin1amask_Nor"]))))/
      median(temp$Cell_Intensity_MedianIntensity_Actin1amask_Nor)
    
    temp$Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr<-temp$Cell_Intensity_MedianIntensity_Actin1amask_Nor*temp$corrfctmda*255
    
    #Integrated
    temp$corrfctinta<- (mean(c(median(perobiint[perobiint$Image_Metadata_array==1,"Cell_Intensity_IntegratedIntensity_Actin1amask_Nor"]),
                              median(perobiint[perobiint$Image_Metadata_array==2,"Cell_Intensity_IntegratedIntensity_Actin1amask_Nor"]),
                              median(perobiint[perobiint$Image_Metadata_array==3,"Cell_Intensity_IntegratedIntensity_Actin1amask_Nor"]))))/
      median(temp$Cell_Intensity_IntegratedIntensity_Actin1amask_Nor)
    
    temp$Cell_Intensity_IntegratedIntensity_Actin1amask_Nor_corr<-temp$Cell_Intensity_IntegratedIntensity_Actin1amask_Nor*temp$corrfctinta*255
    
    
  }else{
    temp<-perobiint[perobiint$Image_Metadata_array==i,]
    temp$corrfct<- 1
    temp$corrfctmd<- 1
    temp$corrfctint<- 1
    temp$corrfcta<- 1
    temp$corrfctmda<- 1
    temp$corrfctinta<- 1
    temp$Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr<-temp$Cell_Intensity_MeanIntensity_Icam1amask_Nor*temp$corrfct*255
    temp$Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr<-temp$Cell_Intensity_MedianIntensity_Icam1amask_Nor*temp$corrfctmd*255
    temp$Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr<-temp$Cell_Intensity_IntegratedIntensity_Icam1amask_Nor*temp$corrfctint*255
    
    temp$Cell_Intensity_MeanIntensity_Actin1amask_Nor_corr<-temp$Cell_Intensity_MeanIntensity_Actin1amask_Nor*temp$corrfcta*255
    temp$Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr<-temp$Cell_Intensity_MedianIntensity_Actin1amask_Nor*temp$corrfctmda*255
    temp$Cell_Intensity_IntegratedIntensity_Actin1amask_Nor_corr<-temp$Cell_Intensity_IntegratedIntensity_Actin1amask_Nor*temp$corrfctinta*255
    
  }
  perobiintc<-rbind(perobiintc,temp)
}
##visualize changes
#Icam
#Median Intensity
ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_MedianIntensity_Icam1amask_Nor_corr,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_MedianIntensity_Icam1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
#Mean intenssity
ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_MeanIntensity_Icam1amask_Nor_corr,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_MeanIntensity_Icam1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
#Integrated Intensity
ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_IntegratedIntensity_Icam1amask_Nor_corr,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_IntegratedIntensity_Icam1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
#Actin
#median intensity
ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_MedianIntensity_Actin1amask_Nor_corr,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_MedianIntensity_Actin1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
##mean intensity
ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_MeanIntensity_Actin1amask_Nor_corr,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_MeanIntensity_Actin1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()
##integrated intensity
ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_IntegratedIntensity_Actin1amask_Nor_corr,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

ggplot(perobiintc,aes(factor(Image_Metadata_array), Cell_Intensity_IntegratedIntensity_Actin1amask_Nor,
                     fill=factor(Image_Metadata_array)))+geom_boxplot()

#save results
library(data.table)

perobindall.cor.temp1<-na.omit(perobindall.filter[,c("Cell_Intensity_IntegratedIntensity_Actin1amask",
    "Cell_Intensity_MedianIntensity_Actin1amask","Cell_Intensity_MeanIntensity_Actin1amask",
           "Cell_Intensity_IntegratedIntensity_Icam1amask",
        "Cell_Intensity_MedianIntensity_Icam1amask","Cell_Intensity_MeanIntensity_Icam1amask",      
    colnames( perobindall.filter)[!grepl("_Intensity_", colnames( perobindall.filter))])])
perobindall.cor.temp1$ID<-paste(perobindall.cor.temp1$ImageNumber,perobindall.cor.temp1$ObjectNumber,sep="_")

perobindall.cor.temp2<-perobiintc[,c("ImageNumber","ObjectNumber",colnames(perobiintc)[grepl("_Nor_", 
                    colnames(perobiintc))])]
perobindall.cor.temp2$ID<-paste(perobindall.cor.temp2$ImageNumber,perobindall.cor.temp2$ObjectNumber,sep="_")
perobindall.cor.temp2<-perobindall.cor.temp2[,!colnames(perobindall.cor.temp2)%in%c("ImageNumber","ObjectNumber")]

#perform some checks
length(unique(perobindall.cor.temp1$ID))
length(unique(perobindall.cor.temp2$ID))
intersect(colnames(perobindall.cor.temp1),colnames(perobindall.cor.temp2))

perobindall.cor<-merge(perobindall.cor.temp1, perobindall.cor.temp2, by="ID")
##fixing colnames
#colnames(perobindall.cor)<-gsub("")

save(perobindall.cor,file="Icam_cells_after_correction.RDATA")
