#batch plot, 2016Feb, testing p-value for 
#if We can estimate p-value of observed ¾ blacks from ¼ and ½ black distributions. 
#assuming Joint poisson distribution 

# used the same script based on 2013 Dec 10, adjust for low-counts
#tb$Black[tb$Black<=0] = 0.5
#tb$halfBlack[tb$halfBlack<=0] = 0.5

rm=(list=ls())
#setwd("~/projects/LOH-oxidants2012.osX/analysis")
setwd("~/github/LOH_H2O2_2016/analysis")
debug = 0;

FileList = list.files( path="../data.H2O2-LOH/");  FileList; 
if( debug > 5) {FileList = FileList[1:2]}

#ListofFiles = list()
Poisson_Results =list()
Poisson_pValues = list()
#tbAllfiles = list()
for( infile in FileList) {
  
  # infile='M1-2,20111207.modified.csv' #debug for NA, and low lead concentration
  # infile =  'M2-8,08172011.H2O2.LOH.modified.csv';  #viability plot error
  #infile =  'M2-8,06062011.H2O2onLOH.csv';  #b plot error
  
  fullFileName = paste('../data.H2O2-LOH/',infile, sep='');
  mylabel = infile
  
  
  tb = read.csv(fullFileName, colClasses=c("character",NA, NA, "character", rep("numeric",8 ), NA));
  names(tb) = c("Strain", "OD600", "Dilution","Date","H2O2stock", "White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")
  
  ######## set zeros
  mycolumns = c("White","Black","halfBlack", "quarterBlack","ThreeQBlack", "QQBlack", "Other"); 
  for( i in 1:length(tb[,1])) {
    for ( j in mycolumns) {
      if( is.na(tb[i,j]) ) { tb[i,j]= 0 }
    }
  }
  
  #adjust for low-counts
  tb$Black[tb$Black<=0] = 0.5
  tb$halfBlack[tb$halfBlack<=0] = 0.5
  
  tb$H2O2 = tb$H2O2stock/2
  tb$tot = tb$White + tb$Black + tb$halfBlack + tb$quarterBlack + tb$ThreeQBlack + tb$QQBlack + tb$Other
  tb.ori = tb; 
  tb = tb[ ! is.na(tb$White), ]
  
  
  tb$Dilution = tb$Dilution / tb$Dilution[1]
  
  ###### some manual curations here
  #tb$Dilution[c(19,20,21)] = c(10,10,10)  ##this dilution should be 10 times more
  #tb$Black[25] = NA #This number is not right!!!!
  
  ######## normalize all data
  mycolumns = c("White","Black","halfBlack", "quarterBlack","ThreeQBlack", "QQBlack", "Other","tot"); 
  for ( j in mycolumns) {
    tb[,j] = tb[,j] * tb$Dilution
    
    
  }
  
  ####### find out means
  H2O2 = sort(unique( tb$H2O2))
  #s = H2O2
  tbm = data.frame(cbind(H2O2))
  for ( i in 1:length(H2O2)) {
    c = H2O2[i]
    tmp = tb[ tb$H2O2==c, ]  
    tbm$tot[i] = mean(tmp$tot, na.rm=T)
    tbm$White[i] = mean(tmp$White, na.rm=T)
    tbm$Black[i] = mean(tmp$Black, na.rm=T) 
    tbm$halfBlack[i] = mean(tmp$halfBlack, na.rm=T)  
    tbm$quarterBlack[i] = mean(tmp$quarterBlack, na.rm=T)
    tbm$ThreeQBlack[i] = mean(tmp$ThreeQBlack, na.rm=T)  
    tbm$QQBlack[i] = mean(tmp$QQBlack, na.rm=T);  
  }
  
  
  Black = sort(unique( tb$H2O2))
  #s = H2O2
  tbm = data.frame(cbind(H2O2))
  for ( i in 1:length(H2O2)) {
    c = H2O2[i]
    tmp = tb[ tb$H2O2==c, ]  
    tbm$tot[i] = mean(tmp$tot, na.rm=T)
    tbm$White[i] = mean(tmp$White, na.rm=T)
    tbm$Black[i] = mean(tmp$Black, na.rm=T) 
    tbm$halfBlack[i] = mean(tmp$halfBlack, na.rm=T)  
    tbm$quarterBlack[i] = mean(tmp$quarterBlack, na.rm=T)
    tbm$ThreeQBlack[i] = mean(tmp$ThreeQBlack, na.rm=T)  
    tbm$QQBlack[i] = mean(tmp$QQBlack, na.rm=T);  
    
  }
  
    pois_halfBlack<-dpois(tb$halfBlack, mean(tb$halfBlack), log =F)
   pois_quarterBlack<-dpois(tb$quarterBlack, mean(tb$quarterBlack), log =F)
   pois_ThreeQBlack<-dpois(tb$ThreeQBlack, mean(tb$ThreeQBlack), log =F)
  
  Result_pois<-pois_halfBlack*pois_quarterBlack-pois_ThreeQBlack;
  
  for (i in 1:length(Result_pois)) {
  if ( is.na(Result_pois[i])) { Result_pois[i] = 0 }
  }
  
  test_result<-t.test(Result_pois)$p.value
  Poisson_Results[[length(Poisson_Results)+1]] = Result_pois;
  Poisson_pValues[[length(Poisson_pValues)+1]] = test_result;
  
  } 

 for (i in length(Poisson_pValues)){
 m_pvalue<-sum(Poisson_pValues[[i]])/length(Poisson_pValues)
 m_pvalue
 }


