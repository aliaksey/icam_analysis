rm(list=ls())
load("cell_objects_after_area_perimeter.RDATA")
library(chemometrics)
set.seed(763298574)

#select only topography

sh.in.mahal<-perobindall.filter[,c("ImageNumber", "FeatureIdx","ObjectNumber","Image_Metadata_array",
        colnames(perobindall.filter)[grepl("AreaShape",colnames(perobindall.filter))&
                                       !grepl("Center",colnames(perobindall.filter))&
                                       !grepl("Euler",colnames(perobindall.filter))])]
sh.in.mahal.t<-sh.in.mahal[sh.in.mahal$Image_Metadata_array<9,]
#finding outlieras in Mahalanobis'distance based on all available shape parameters ~30

sh.in.mahal.t.f<-c()
for(i in unique(sh.in.mahal.t[,"FeatureIdx"])){
  temp<-sh.in.mahal.t[sh.in.mahal.t$FeatureIdx==i,]
  temp2<-temp[,colnames(temp)[grepl("AreaShape",colnames(temp))]]
  ###filter based on cell area   
  mdres<-Moutlier(temp2, quantile = 0.99, plot = T)
  rsltmd<-temp2[mdres$rd < mdres$cutoff,]
  arpr.ftr<-temp2[row.names(temp2) %in% row.names(rsltmd),]
  sh.in.mahal.t.f<-rbind(sh.in.mahal.t.f,arpr.ftr) 
}
sh.in.mahal.t.f<-as.data.frame(sh.in.mahal.t.f)
nrow(sh.in.mahal.t.f)
