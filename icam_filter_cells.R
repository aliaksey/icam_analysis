rm(list=ls())
load("corrim.RData")
#plotting cell area/ perimeter
library(ggplot2)
ggplot(perobindall.cor,aes(x=Cell_AreaShape_Area,y=Cell_AreaShape_Perimeter,colour=Image_Metadata_array))+geom_point()
##  plotting only topochip data
ggplot(perobindall.cor[perobindall.cor$Image_Metadata_array<9,],
       aes(x=Cell_AreaShape_Area,y=Cell_AreaShape_Perimeter,colour=Image_Metadata_array))+geom_point()

ggplot(perobindall.cor[perobindall.cor$Image_Metadata_array<9,],
       aes(x=Nuclei_AreaShape_Area,y=Nuclei_AreaShape_Perimeter,colour=Image_Metadata_array))+geom_point()
###  filter area perimeter per feature

cell.area<-na.omit(perobindall.cor[perobindall.cor$Image_Metadata_array<9,c("ImageNumber", "FeatureIdx", "ObjectNumber",
                                                                    "Cell_AreaShape_Area","Cell_AreaShape_Perimeter",
                                                                    "Nuclei_AreaShape_Area","Nuclei_AreaShape_Perimeter",
                                                                    "Cell_AreaShape_Solidity","Cell_AreaShape_Compactness")])
# library(GGally)
# # library(gplots)
# # palette(rev(rich.colors(2177)))
# ggpairs(cell.area,columns=4:7)
# palette("default")
##
cell.area.f<-c()
# library(parallel)
# cl <- makeCluster(4)
# registerDoParallel(cl)

for(i in unique(cell.area[,"FeatureIdx"])){
  temp2<-cell.area[cell.area$FeatureIdx==i,]
  ###filter based on cell area   
  lbnda<-as.numeric(quantile(temp2[,"Cell_AreaShape_Area"], probs = 0.25))
  ubnda<-as.numeric(quantile(temp2[,"Cell_AreaShape_Area"], probs = 0.75))
  iuda<-ubnda-lbnda
  rsltac<-temp2[temp2[,"Cell_AreaShape_Area"]<(ubnda+1.5*iuda)&
                 temp2[,"Cell_AreaShape_Area"]>(lbnda-1.5*iuda),]
  ###filter based on cell perimeter
  lbndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Perimeter"], probs = 0.25))
  ubndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Perimeter"], probs = 0.75))
  iudp<-ubndp-lbndp
  rsltpc<-temp2[temp2[,"Cell_AreaShape_Perimeter"]<(ubndp+1.5*iudp)&
                 temp2[,"Cell_AreaShape_Perimeter"]>(lbndp-1.5*iudp),]
  
  ###filter based on cell compactness
  lbndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Compactness"], probs = 0.25))
  ubndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Compactness"], probs = 0.75))
  iudp<-ubndp-lbndp
  rsltcc<-temp2[temp2[,"Cell_AreaShape_Compactness"]<(ubndp+1.5*iudp)&
                  temp2[,"Cell_AreaShape_Compactness"]>(lbndp-1.5*iudp),]
  
  ###filter based on cell solidity
  lbndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Solidity"], probs = 0.25))
  ubndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Solidity"], probs = 0.75))
  iudp<-ubndp-lbndp
  rsltsc<-temp2[temp2[,"Cell_AreaShape_Solidity"]<(ubndp+1.5*iudp)&
                  temp2[,"Cell_AreaShape_Solidity"]>(lbndp-1.5*iudp),]
  
  ###filter based on nuclei area   
  lbnda<-as.numeric(quantile(temp2[,"Nuclei_AreaShape_Area"], probs = 0.25))
  ubnda<-as.numeric(quantile(temp2[,"Nuclei_AreaShape_Area"], probs = 0.75))
  iuda<-ubnda-lbnda
  rsltan<-temp2[temp2[,"Nuclei_AreaShape_Area"]<(ubnda+1.5*iuda)&
                  temp2[,"Nuclei_AreaShape_Area"]>(lbnda-1.5*iuda),]
  ###filter based on nuclei perimeter
  lbndp<-as.numeric(quantile(temp2[,"Nuclei_AreaShape_Perimeter"], probs = 0.25))
  ubndp<-as.numeric(quantile(temp2[,"Nuclei_AreaShape_Perimeter"], probs = 0.75))
  iudp<-ubndp-lbndp
  rsltpn<-temp2[temp2[,"Nuclei_AreaShape_Perimeter"]<(ubndp+1.5*iudp)&
                  temp2[,"Nuclei_AreaShape_Perimeter"]>(lbndp-1.5*iudp),]
  
  
  arpr.ftr<-temp2[row.names(temp2) %in% row.names(rsltac)&
                    row.names(temp2) %in% row.names(rsltpc)&
                    row.names(temp2) %in% row.names(rsltan)&
                    row.names(temp2) %in% row.names(rsltsc)&
                    row.names(temp2) %in% row.names(rsltcc)&
                    row.names(temp2) %in% row.names(rsltpn),]
  cell.area.f<-rbind(cell.area.f,arpr.ftr) 
}
# stopCluster(cl)
# rm(cl)
cell.area.f<-as.data.frame(cell.area.f)
cell.area.f$ID<-paste(cell.area.f$ImageNumber,cell.area.f$ObjectNumber,sep="_")
##joining filtering results with a data
perobindall.cor$ID<-paste(perobindall.cor$ImageNumber,perobindall.cor$ObjectNumber,sep="_")
length(perobindall.cor$ID)

perobindall.filter<-perobindall.cor[perobindall.cor$ID%in%cell.area.f$ID,]

ggplot(perobindall.filter[perobindall.filter$Image_Metadata_array<9,],
       aes(x=Cell_AreaShape_Area,y=Cell_AreaShape_Perimeter,colour=Image_Metadata_array))+geom_point()

ggplot(perobindall.filter[perobindall.filter$Image_Metadata_array<9,],
       aes(x=Nuclei_AreaShape_Area,y=Nuclei_AreaShape_Perimeter,colour=Image_Metadata_array))+geom_point()

##checking with other parameters compactness eccentricity
ggplot(perobindall.filter[perobindall.filter$Image_Metadata_array<9,],
       aes(x=Cell_AreaShape_Solidity,y=Nuclei_AreaShape_Compactness,colour=Image_Metadata_array))+geom_point()

ggplot(perobindall.cor[perobindall.cor$Image_Metadata_array<9,],
       aes(x=Cell_AreaShape_Solidity,y=Nuclei_AreaShape_Compactness,colour=Image_Metadata_array))+geom_point()

##checking with other parameters formfactor solidity
ggplot(perobindall.filter[perobindall.filter$Image_Metadata_array<9,],
       aes(x=Cell_AreaShape_FormFactor,y=Nuclei_AreaShape_Eccentricity,colour=Image_Metadata_array))+geom_point()

ggplot(perobindall.cor[perobindall.cor$Image_Metadata_array<9,],
       aes(x=Cell_AreaShape_FormFactor,y=Nuclei_AreaShape_Eccentricity,colour=Image_Metadata_array))+geom_point()
#apply these filters for controls

cell.area.c<-na.omit(perobindall.cor[perobindall.cor$Image_Metadata_array>8,c("ImageNumber","Col","Row", "Image_Metadata_array" , "ObjectNumber",
                                                                            "Cell_AreaShape_Area","Cell_AreaShape_Perimeter",
                                                                            "Nuclei_AreaShape_Area","Nuclei_AreaShape_Perimeter",
                                                                            "Cell_AreaShape_Solidity","Cell_AreaShape_Compactness")])
cell.area.c$CONID<-paste(cell.area.c$Image_Metadata_array,
                         cell.area.c$Col,sep="_")

cell.area.f.c<-c()
# library(parallel)
# cl <- makeCluster(4)
# registerDoParallel(cl)
for(i in unique(cell.area.c[,"CONID" ])){
  temp2<-cell.area.c[cell.area.c$CONID==i,]
  ###filter based on cell area   
  lbnda<-as.numeric(quantile(temp2[,"Cell_AreaShape_Area"], probs = 0.25))
  ubnda<-as.numeric(quantile(temp2[,"Cell_AreaShape_Area"], probs = 0.75))
  iuda<-ubnda-lbnda
  rsltac<-temp2[temp2[,"Cell_AreaShape_Area"]<(ubnda+1.5*iuda)&
                  temp2[,"Cell_AreaShape_Area"]>(lbnda-1.5*iuda),]
  ###filter based on cell perimeter
  lbndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Perimeter"], probs = 0.25))
  ubndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Perimeter"], probs = 0.75))
  iudp<-ubndp-lbndp
  rsltpc<-temp2[temp2[,"Cell_AreaShape_Perimeter"]<(ubndp+1.5*iudp)&
                  temp2[,"Cell_AreaShape_Perimeter"]>(lbndp-1.5*iudp),]
  ###filter based on cell compactness
  lbndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Compactness"], probs = 0.25))
  ubndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Compactness"], probs = 0.75))
  iudp<-ubndp-lbndp
  rsltcc<-temp2[temp2[,"Cell_AreaShape_Compactness"]<(ubndp+1.5*iudp)&
                  temp2[,"Cell_AreaShape_Compactness"]>(lbndp-1.5*iudp),]
  
  ###filter based on cell solidity
  lbndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Solidity"], probs = 0.25))
  ubndp<-as.numeric(quantile(temp2[,"Cell_AreaShape_Solidity"], probs = 0.75))
  iudp<-ubndp-lbndp
  rsltsc<-temp2[temp2[,"Cell_AreaShape_Solidity"]<(ubndp+1.5*iudp)&
                  temp2[,"Cell_AreaShape_Solidity"]>(lbndp-1.5*iudp),]
  ###filter based on nuclei area   
  lbnda<-as.numeric(quantile(temp2[,"Nuclei_AreaShape_Area"], probs = 0.25))
  ubnda<-as.numeric(quantile(temp2[,"Nuclei_AreaShape_Area"], probs = 0.75))
  iuda<-ubnda-lbnda
  rsltan<-temp2[temp2[,"Nuclei_AreaShape_Area"]<(ubnda+1.5*iuda)&
                  temp2[,"Nuclei_AreaShape_Area"]>(lbnda-1.5*iuda),]
  ###filter based on nuclei perimeter
  lbndp<-as.numeric(quantile(temp2[,"Nuclei_AreaShape_Perimeter"], probs = 0.25))
  ubndp<-as.numeric(quantile(temp2[,"Nuclei_AreaShape_Perimeter"], probs = 0.75))
  iudp<-ubndp-lbndp
  rsltpn<-temp2[temp2[,"Nuclei_AreaShape_Perimeter"]<(ubndp+1.5*iudp)&
                  temp2[,"Nuclei_AreaShape_Perimeter"]>(lbndp-1.5*iudp),]
  
  
  arpr.ftr<-temp2[row.names(temp2) %in% row.names(rsltac)&
                    row.names(temp2) %in% row.names(rsltpc)&
                    row.names(temp2) %in% row.names(rsltan)&
                    row.names(temp2) %in% row.names(rsltsc)&
                    row.names(temp2) %in% row.names(rsltcc)&
                    row.names(temp2) %in% row.names(rsltpn),]
  cell.area.f.c<-rbind(cell.area.f.c,arpr.ftr) 
}
# stopCluster(cl)
# rm(cl)
cell.area.f.c<-as.data.frame(cell.area.f.c)
cell.area.f.c$ID<-paste(cell.area.f.c$ImageNumber,cell.area.f.c$ObjectNumber,sep="_")
##joining filtering results with a data

perobindall.filter.c<-perobindall.cor[perobindall.cor$ID%in%cell.area.f.c$ID,]

ggplot(perobindall.filter.c[perobindall.filter.c$Image_Metadata_array>8,],
       aes(x=Cell_AreaShape_Area,y=Cell_AreaShape_Perimeter,colour=Image_Metadata_array))+geom_point()

ggplot(perobindall.filter.c[perobindall.filter.c$Image_Metadata_array>8,],
       aes(x=Nuclei_AreaShape_Area,y=Nuclei_AreaShape_Perimeter,colour=Image_Metadata_array))+geom_point()

ggplot(perobindall.cor[perobindall.cor$Image_Metadata_array>8,],
       aes(x=Cell_AreaShape_Area,y=Cell_AreaShape_Perimeter,colour=Image_Metadata_array))+geom_point()

ggplot(perobindall.cor[perobindall.cor$Image_Metadata_array>8,],
       aes(x=Nuclei_AreaShape_Area,y=Nuclei_AreaShape_Perimeter,colour=Image_Metadata_array))+geom_point()

## joining all data together

perobindall.filter<-rbind(perobindall.filter,perobindall.filter.c)
##remove objects that have area=0
perobindall.filter<-perobindall.filter[perobindall.filter$Cell_AreaShape_Area!=0,]

save(perobindall.filter, file="cell_objects_after_area_perimeter.RDATA")
