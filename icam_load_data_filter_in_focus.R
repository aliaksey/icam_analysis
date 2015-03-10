trval<-0.29 #set treshhold value for posotive cells
load("data/icamimageandobjectindex.RData")
#selecting only variables that dealing with focus

perimindallf<-perimindall[,c("ImageNumber","FeatureIdx","InternalIdx","Image_Metadata_array",
              colnames(perimindall)[grepl("Focus",colnames(perimindall))],
              colnames(perimindall)[grepl("focus",colnames(perimindall))],
              colnames(perimindall)[grepl("LogLogSlope",colnames(perimindall))])]
LogLogSlope     

perimindallftt<-perimindall[perimindall$Image_Metadata_array>8,]


#finding out of focus images
library (ggplot2)

highlight.topo <- 2177
perobindallfplot <- as.data.frame(perimindallf)

perobindallfplot$highlight <- ifelse(perobindallfplot$FeatureIdx == highlight.topo, "highlight", "normal")
textdf <- perobindallfplot[perobindallfplot$FeatureIdx == highlight.topo, ]
mycolours <- c("highlight" = "red", "normal" = "grey50")

p11<-ggplot(data = perobindallfplot, aes(x = Image_Math_focusa,
                                         y = Image_Math_focusd)) 
p11+ geom_point(alpha = 0.4) +geom_density2d()+xlab("log_actin_local_focus_meas")+
  ylab("log_dapi_local_focus_meas")+ geom_segment(aes(x = -4.2, y = -5, xend = -4.2, yend = 2, colour="red"))+
  geom_segment(aes(x = -6, y = -2.5, xend = 1, yend = -2.5, colour="red"))

# +geom_density2d()+
#   ylim(0,0.1)+xlim(0,0.0025)

library(hexbin)
p11+stat_binhex()
p11+geom_point(alpha = 0.4)
p11 + geom_point() + geom_density2d()
# p11+ geom_point(size = 3, aes(colour = highlight)) +
#   scale_color_manual("Status", values = mycolours) +
#   geom_text(data = textdf, aes(x = Index1, y = Index2, label = "blank")) +
#   theme(legend.position = "none") +
#   theme()

focim<-perimindallf[perimindallf$Image_Math_focusa>-4.2&perimindallf$Image_Math_focusd>-2.5,"ImageNumber"]
# finding images without over exposion 
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