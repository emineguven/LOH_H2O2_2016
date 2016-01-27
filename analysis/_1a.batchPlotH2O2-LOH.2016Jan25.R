# 20160125, try to add error bars to plots

rm=(list=ls())
setwd("~/github/LOH_H2O2_2016/analysis")
debug = 0;

#all data are in 'data.H2O2-LOH'
FileList = list.files( path="../data.H2O2-LOH/");  FileList; 
if( debug > 5) {FileList = FileList[1:2]}

# infile = FileList[2]
for( infile in FileList) {
# infile='M1-2,20111207.modified.csv' #debug for NA, and low lead concentration
# infile =  'M2-8,08172011.H2O2.LOH.modified.csv';  #viability plot error
# infile =  'M2-8,06062011.H2O2onLOH.csv';  #b plot error
   
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

#H2O2 working stock is 2X 
tb$H2O2 = tb$H2O2stock/2
tb$tot = tb$White + tb$Black + tb$halfBlack + tb$quarterBlack + tb$ThreeQBlack + tb$QQBlack + tb$Other
tb.ori = tb; 
tb = tb[ ! is.na(tb$White), ]

tb$Dilution = tb$Dilution / tb$Dilution[1]

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
  tbm$tot.sd[i] = sd(tmp$tot, na.rm=T)
  tbm$White[i] = mean(tmp$White, na.rm=T)
  tbm$White.sd[i] = sd(tmp$White, na.rm=T)
  tbm$Black[i] = mean(tmp$Black, na.rm=T) 
  tbm$Black.sd[i] = sd(tmp$Black, na.rm=T) 
  tbm$halfBlack[i] = mean(tmp$halfBlack, na.rm=T)  
  tbm$halfBlack.sd[i] = sd(tmp$halfBlack, na.rm=T)  
  tbm$quarterBlack[i] = mean(tmp$quarterBlack, na.rm=T)
  tbm$ThreeQBlack[i] = mean(tmp$ThreeQBlack, na.rm=T)  
  tbm$QQBlack[i] = mean(tmp$QQBlack, na.rm=T);  
}

tbm = tbm[tbm$tot>1, ] #remove plates with zero colonies

###### some manual curations here
#tbm$halfBlack[tbm$halfBlack==0 & tbm$H2O2>0] = NA;

###### calculate fractions
tbf = tbm; 
tbf$s = tbf$tot / max(tbf$tot)
tbf$s.sd = tbm$tot.sd / max(tbm$tot)
for ( j in 3:11) {
  tbf[, j] = tbf[,j] / tbf$tot
}

tbf$Black[tbf$Black<0]=NA;  #remove weird experimental data, such as low-lead concentration effect

pdf(paste("sandbox/",infile, ".batch.sd.pdf", sep=''), width=5, height=5)

### full black plot
plot( tbf$s ~ tbf$H2O2, col="blue", axes=F, xlab='', ylab='', ylim=c(-0.05, 1.2), type='p');
lines( tbf$s ~ tbf$H2O2, col="blue",lty=2)
arrows( tbf$H2O2, (tbf$s - tbf$s.sd), tbf$H2O2, (tbf$s + tbf$s.sd), length=0.1, angle=90,code=3, lty=1, col="blue" );

axis( 4, at=pretty(c(0, 1.2)), col='black')
legend ( max(H2O2)*0.7, 0.9, c("viability","black"), col=c("blue","black"), lty=c(2,1), pch=c(1,16) )
par(new=T)
plot( tbf$Black ~ tbf$H2O2, pch=16, xlab='H2O2',ylab="black", ylim=c(-0.01, max(tbf$Black, na.rm = T)*3))
lines(tbf$Black ~ tbf$H2O2)
bottoms = tbf$Black -tbf$Black.sd
bottoms[bottoms<0]= 0.001 
arrows( tbf$H2O2, bottoms, tbf$H2O2, (tbf$Black + tbf$Black.sd), length=0.1, angle=90,code=3, lty=1, col="gray" );

title(mylabel)

## full black -log plot
#tbf$H2O2[tbf$H2O2==0] = min(H2O2[-1])/10;

#with( tbf, plot( tot ~ log10(H2O2), col="blue"));
#par(new=T)
#with( tbf, plot( Black ~ log10(H2O2), pch=16))

#tbf$H2O2[tbf$H2O2==0] = min(H2O2[-1])/10;
#with( tbf, plot( s ~ log10(H2O2), col="blue", axes=F, xlab="H2O2", ylab='viability'), ylim=c(-0.5, 1.1) );
#xlabels = log10(c(0.001,0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.15))
#axis(1, at = xlabels, labels= 10^xlabels)
#ylabels = c(0, 0.25, 0.5, 0.75, 1.0)
#axis(2, at=ylabels, labels=ylabels);
#par(new=T)
#with( tbf, plot( Black ~ log10(H2O2), pch=16, axes=F, xlab='', ylab=''))
#axis(4, pretty(range(tbf$Black)))
#title(mylabel)

### half black plot
plot( tbf$s ~ tbf$H2O2, col="blue", axes=F, xlab='', ylab='', ylim=c(0, 1.2), type='p', main=mylabel);
lines( tbf$s ~ tbf$H2O2, col="blue",lty=2)
axis( 4, at=pretty(c(0, 1.2)))
legend ( max(H2O2)*0.7, 0.5, c("viability","half-black"), col=c("blue","red"), lty=c(2,1), pch=c(1,16) )
par(new=T)
plot( tbf$halfBlack ~ tbf$H2O2, pch=16, xlab='H2O2',ylab="half-black", col='red', ylim=c(-0.005, max(tbf$halfBlack, na.rm=T)*2 ))
lines(tbf$halfBlack ~ tbf$H2O2, col='red')

bottoms = tbf$halfBlack - tbf$halfBlack.sd
bottoms[bottoms<0]= 0
arrows( tbf$H2O2, bottoms, tbf$H2O2, (tbf$halfBlack + tbf$halfBlack.sd), length=0.1, angle=90,code=3, lty=1, col="gray" );

title(mylabel)


### half/full plot
tbf$half.vs.full = tbf$halfBlack / tbf$Black
plot( tbf$s ~ tbf$H2O2, col="blue", axes=F, xlab='', ylab='', ylim=c(0, 1.2), type='p', main=mylabel);
lines( tbf$s ~ tbf$H2O2, col="blue",lty=2)
axis( 4, at=pretty(c(0, 1.2)))
legend ( max(H2O2)*0.7, 0.5, c("viability","half/full"), col=c("blue","green"), lty=c(2,1), pch=c(1,16) )
par(new=T)
plot( tbf$half.vs.full ~ tbf$H2O2, pch=16, xlab='H2O2',ylab="half/full", col='green')
lines(tbf$half.vs.full ~ tbf$H2O2, col='green')
bottoms = tbf$Black -tbf$Black.sd
bottoms[bottoms<0]= 0.001 
title(mylabel)

dev.off() #end pdf 	
} #infile loop
	
#quit("yes")	
   
   
