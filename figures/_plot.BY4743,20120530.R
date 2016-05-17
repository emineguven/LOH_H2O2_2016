# 2013 Dec 10, use layout and mar() to do layered plots (not overlays)

#YPS163 2011April 13, plot

rm=(list=ls())
setwd("~/github/LOH_H2O2_2012-master/figures")
debug = 0;

FileList = list.files( path="../data.H2O2-LOH/", pattern="BY4743");  FileList; 
filename = 'BY4743,20120530,H2O2LOH.csv';
mylabel = 'BY4743 2013 May 30'

#filename = 'BY4743,20120907,H2O2LOH.csv';
#mylabel = 'BY4743 2013 Sep 7'

fullFileName = paste('../data.H2O2-LOH/',filename, sep='');

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

tb$Black[tb$Black<0]=NA;  #remove weird experimental data, such as low-lead concentration effect

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

####### find out means, sd
H2O2 = sort(unique( tb$H2O2))
#s = H2O2
tbm = data.frame(cbind(H2O2))
for ( i in 1:length(H2O2)) {
  c = H2O2[i]
  tmp = tb[ tb$H2O2==c, ]	
  tbm$tot[i] = mean(tmp$tot, na.rm=T)
  tbm$tot.sd[i] = sd(tmp$tot, na.rm=T)
  tbm$White[i] = mean(tmp$White, na.rm=T)
  tbm$Black[i] = mean(tmp$Black, na.rm=T)
  tbm$Black.sd[i] = sd(tmp$Black, na.rm=T)
  tbm$halfBlack[i] = mean(tmp$halfBlack, na.rm=T)  
  tbm$halfBlack.sd[i] = sd(tmp$halfBlack, na.rm=T)  
  tbm$quarterBlack[i] = mean(tmp$quarterBlack, na.rm=T)
  tbm$ThreeQBlack[i] = mean(tmp$ThreeQBlack, na.rm=T)  
  tbm$QQBlack[i] = mean(tmp$QQBlack, na.rm=T);  
}

tbm = tbm[tbm$tot>2, ] #remove plates with <2 colonies. At least 1 white and 1 black is needed. 

###### some manual curations here
# ... 

###### calculate fractions
tbf = tbm; 
tbf$s = tbf$tot / max(tbf$tot)
for ( j in c("Black","halfBlack", "quarterBlack","ThreeQBlack","QQBlack")) {
  tbf[, j] = tbf[,j] / tbf$tot
}
tbf$tot.sd = tbf$tot.sd / max(tbf$tot); #20131211
for ( j in c("Black.sd","halfBlack.sd")) {
 tbf[, j] = tbf[,j] / tbf$tot; #change 20131211
  tbf[, j] = ifelse(tbf[, j]<=1/700, 1/700, tbf[, j]  ) #at least 0.2% error
}
head(tbf)
plot(tbf$tot.sd ~ tbf$H2O2)


#It is not clear whether color figures would be charged more or not by Aging Cell. 

#### define functions
logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }
logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }
# genome.integrity <- function(b.max, b.min, T, w, t) { 2 * (b.max - b.min) / (1 + (t/T)^w) + (1 - 2 * b.max); }

######## calclulate Cv, Cb
require(nlme)
t= tbf$H2O2; s= tbf$s; b = tbf$Black

model1 = function(v,w) {  1/( 1 + ( t / v )^ w ) }  
fm.s <- gnls( s ~ model1(v,w) , start = list( v = 0.02, w=2)  );
t2 = seq(0,1,by=0.001)
fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t2 );
Cv = fm.s$coefficients[1]

#estimate starting values
b.tmp = sort( b[ b>0 ] );
b.min = sum(b.tmp[1])*1.2
b.tmp = rev(b.tmp); 
b.max = sum(b.tmp[1:4])/4
  
modelBlack = function(v, w) { b.max - (b.max-b.min)/(1 + (t/v)^w) }
modelBlack1 = function(v) { b.max - (b.max-b.min)/(1 + (t/v)^0.4) }
#
fm.b <- gnls( b ~ modelBlack1(v) , start = list( v = 0.008) ); fm.b
fit.b = logistical.black(b.max, b.min, fm.b$coef[1], 2, t2 )

#fm.b <- gnls( b ~ modelBlack(v,w) , start = list( v = 0.0025, w= 0.4  )); fm.b
#fit.b = logistical.black(b.max, b.min, fm.b$coef[1], fm.b$coef[2], t2 )
Cb = fm.b$coef[1]

plot( b ~ t)

####################################################################
#   do the plot
####################################################################

#### specificy the layout
#pdf(paste("_20131210", filename, "pdf", sep="."), width=6,height=9); 
pdf(paste("_20131211", filename, "pdf", sep="."), width=5,height=8); 
#tiff(paste("_20131210", filename, "tif", sep="."), width=240,height=500); 

mat = matrix( seq(1,3), nrow=3, ncol= 1 ); 
layout(mat, heights= c( 1.15, 1, 1.2) );


### viability
par(mar=c(0,5,5,1))
plot( tbf$s ~ tbf$H2O2, col="blue", axes=F, xlab='', ylab='Viability', ylim=c(-0.1, 1.2), type='p', pch=19);
axis( 1, at = pretty(tbf$H2O2), labels=F, tcl=0.3);
axis( 2, at=pretty(tbf$s), tcl=0.2, las=2 )  
lines( fit.s ~ t2, col="blue",lty=1);
#legend ( max(H2O2)*0.5, 0.5, c("viability","black", "half black"), col=c("blue","black"), lty=c(2,1), pch=c(1,16) )
arrows( tbf$H2O2, (tbf$s - tbf$tot.sd), tbf$H2O2, (tbf$s + tbf$tot.sd), length=0.1, angle=90,code=3, lty=2 );
#legend ( max(H2O2)*0.5, 0.5, c("viability"), col=c("blue"), lty=c(1), pch=c(16) )
points( Cv, 0.5, pch=15, col="red", cex=1.2 );
arrows( Cv, 0.5, Cv, -1, lty=2, col="red");
mtext( "Cv",side=1,at=c(Cv*0.85), line=-1, cex=0.8 );
box()
title(mylabel)

############ full blacks
par(mar=c(0,5,0,1))
plot( tbf$Black ~ tbf$H2O2, pch=16, xlab='',ylab="Black", log='y', ylim=c(0.8E-3, 0.15), axes=F)
axis( 2, at = c(0.001,0.01,0.05,0.1), labels=c(0.001,0.01,0.05,0.1), tcl=0.2, las=2);
axis( 1, at = pretty(tbf$H2O2), labels=F, tcl=0.3);
arrows( tbf$H2O2, (tbf$Black - tbf$Black.sd), tbf$H2O2, (tbf$Black + tbf$Black.sd), length=0.1, angle=90,code=3, lty=2 );
lines( fit.b ~ t2, lty=1, col='black');
points ( Cb,  ( b.max/2 + b.min/2), pch=15, col="red", cex=1.2);
arrows( Cb, (b.max/2 + b.min/2), Cb, 1E-8, lty=2, col="red", length=0.1);
mtext( "Cb",side=1,at=c(Cb*1.5), line=-1, cex=0.8 );
box()


####################half black plots

#estimate starting values
b0.5 = tbf$halfBlack
# Do I need to adjust the dilution bias in [0.02,0.03]? No. Plot as is. 
# b0.5[6:8] = b0.5[6:8] + 0.01

b0.5.tmp = sort( b0.5[ b0.5>0 ] );
b0.5.min = sum(b0.5.tmp[1])*1.3
b0.5.tmp = rev(b0.5.tmp); 
b0.5.max = sum(b0.5.tmp[1:3])/3 #account for volatity

modelBlack0.5B = function(v) { b0.5.max - (b0.5.max-b0.5.min)/(1 + (t/v)^2) }
fm.b0.5 <- gnls( b0.5 ~ modelBlack0.5B(v) , start = list( v = 0.008) ); fm.b0.5
fit.b0.5 = logistical.black(b0.5.max, b0.5.min, fm.b0.5$coef[1], 2, t2 )

#logistical.black <- function(b.max, b.min, v, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }
modelBlack0.5 = function(v, w) { b0.5.max - (b0.5.max-b0.5.min)/(1 + (t/v)^w) }
fm.b0.5 <- gnls( b0.5 ~ modelBlack0.5(v,w) , start = list( v = 0.005, w=1)  );   fm.b0.5
fit.b0.5 = logistical.black(b0.5.max, b0.5.min, fm.b0.5$coef[1], fm.b0.5$coef[2], t2 )

#logistical.blackC <- function(b.max, b.min, v, w, u, t ) { ret <- b.max - (b.max - b.min) /( u + ( t / v )^ w );  }
#modelBlack0.5C = function(v, w, u) { b0.5.max - (b0.5.max-b0.5.min)/(u + (t/v)^w) }
#fm.b0.5 <- gnls( b0.5 ~ modelBlack0.5C(v,w,u) , start = list( v = 0.008, w=1.2, u=1.1)  );   fm.b0.5
#fit.b0.5 = logistical.blackC(b0.5.max, b0.5.min, fm.b0.5$coef[1], fm.b0.5$coef[2], fm.b0.5$coef[3], t2 )

#plot(b0.5 ~ t )
#lines(fit.b0.5 ~ t2)
Cb0.5 = fm.b0.5$coef[1]

#par(new=T)
par(mar=c(5,5,0,1))
#plot( tbf$halfBlack ~ tbf$H2O2, pch=16, xlab='H2O2',ylab="half-black", col='red')
plot( b0.5 ~ t, pch=16, xlab='H2O2',ylab="Half-black", col='green', log='y', ylim=c(0.8E-3,0.15), axes=F)
axis( 1, at = pretty(t), labels=T, tcl=0.3);
#axis( 2, at = c(0.001,0.01,0.05,0.1), labels=T, tcl=0.2, las=2);
axis( 2, at = c(0.001,0.01,0.05,0.1), labels=c(0.001,0.01,0.05,0.1), tcl=0.2, las=2);
box()
lines(fit.b0.5 ~ t2, col='green')
arrows( tbf$H2O2, (tbf$halfBlack - tbf$halfBlack.sd), tbf$H2O2, (tbf$halfBlack + tbf$halfBlack.sd),
        length=0.1, angle=90,code=3, lty=2, lwd=1 );
points ( Cb0.5,  ( b0.5.max/2 + b0.5.min/2), pch=15, col="red", cex=1.2);
arrows( Cb0.5, (b0.5.max/2 + b0.5.min/2), Cb0.5, 1E-8, lty=2, col="red", length=0.1);
mtext( "Cb0.5",side=1,at=c(Cb*2.2),line=-1, cex=0.8 );

    
dev.off()