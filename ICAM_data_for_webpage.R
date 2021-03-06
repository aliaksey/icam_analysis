
imagetowebbulk<-function(name,ff,bfidx) for(i in 1:length(bfidx)) ff(paste(name,bfidx[i]),bfidx[i])
featureidx2webpage("1859 best",1859)
featureidx2webpage("1559 best",1559)
featureidx2webpage("572 best",572)
featureidx2webpage("263 best",263)
featureidx2webpage("380 best2",380)
featureidx2webpage("84 best",84)
featureidx2webpage("201 high intensity",201)

featureidx2webpage("1652 some zerow",1652)
featureidx2webpage("integrated intensity hits",forplotcoll)

featureidx2webpage<-function(name, fidx)
{
  load("cell_objects_after_area_perimeter.RDATA")
  load("D:/projects/icam_project/icam_analysis/data/icamimageandobjectindex.RData")
  #   load("focim.RData")
  #   load("corrim.RData")
  htmlfilepath<-paste('D:/projects/icam_project/image_verification/',
                      name,".html", sep = "") 
  ## creating html file
  htmlfile<-file.path(htmlfilepath)
  #specification of html ############################################
  formatspec_head<-'<!DOCTYPE html> \r\n <html> \r\n <body>' 
  formatspec_end<-'</html>' 
  formatspec_h1<-paste('<h1>', name, '</h1>', sep = "")
  formatspec_h1_1<-'<h1> <font color="blue"> ICAm median Intensity </h1>' 
  formatspec_h1_2<-'<h1> Surfaces </h1>' 
  formatspec_h1_3<-'<h1> <font color="blue"> @ </h1>' 
  formatspec_h1_b<-'<h1> <font color="default"> # </h1>' 
  formatspec_h1_n<-'<h1> <font color="red"> 666 </h1>' 
  formatspec_h2<-paste("<h2> Feature number ", fidx, "</h2>\r\n", sep="")   # specify feature number
  formatspec_links<-'<img src='
  formatspec_linke<-' width="300" align="middle"></body> \r\n \r\n' 
  formatspec_links<-'<img src='
  formatspec_h22<-'<h2> Clusters '
  formatspec_h23<-'</h2>\r\n'
  formatspec_h24<-'<h2> Surface Nr '
  formatspec_h25<-'</h2>\r\n'
  formatspec_filt<-'<h2> Filtered cells </h2>\r\n'
  formatspec_left<-'<h2> Remaining cells </h2>\r\n'
  formatspec_left<-'<h2> Remaining cells </h2>\r\n'
  formatspec_linke<-' width="300" align="middle"></body> \r\n \r\n' 
  
  formatspec_<-'<hr>' ########
  ## writing data in html file=======
  cat(formatspec_head,file = htmlfile, append=T)
  cat(formatspec_h1,file = htmlfile, append=T)
  cat(formatspec_,file = htmlfile, append=T)
  #=====
  ## write all
  cat(formatspec_h1_2,file = htmlfile, append=T)
  perimindallinf<-perimindall[perimindall$ImageNumber%in%perobindall.filter$ImageNumber,]
  temp2<-perimindallinf[with(perimindallinf, FeatureIdx %in% fidx$FeatureIdx),]
  temp<-temp2[temp2$Image_Metadata_array<9,]
  ##iterate evrything by class
  for  (cc in unique(fidx$Hit)){
    cat(paste(formatspec_h22,cc,formatspec_h23,sep=""),file = htmlfile, append=T)
    cluster.temp<-fidx[fidx$Hit==cc,]
    #cluster.temp <- cluster.temp.t[sample(1:nrow(cluster.temp.t), 20, replace=T),] 
    tempp<-temp[temp$FeatureIdx %in% cluster.temp$FeatureIdx,]
    for  (ij in unique(cluster.temp$FeatureIdx)){
      cat(paste(formatspec_h24,ij," ",cc,formatspec_h25,sep=""),file = htmlfile, append=T) 
      temppp<-tempp[tempp$FeatureIdx==ij,]
      for  (i in 1:length(temppp[,1])) 
      {
        array<-paste("/array" ,temppp[i,"Image_Metadata_array"],"/", sep="")
        pth0<-paste("../result3",array, temppp[i,"Image_FileName_Dapi1Orig"],"icam_intensity.png", sep="") 
        pth<-gsub("\\.tif","",pth0)
        pth01<-paste("../result3",array, temppp[i,"Image_FileName_Dapi1Orig"],"outlines.png", sep="") 
        pth1<-gsub("\\.tif","",pth01)
        cat(paste(formatspec_links,pth,formatspec_linke,sep=""),file = htmlfile, append=T);
        cat(paste(formatspec_links,pth1,formatspec_linke,sep=""),file = htmlfile, append=T);
      } 
    }
  }
    cat(formatspec_end,file = htmlfile, append=T)
  } 
  
  