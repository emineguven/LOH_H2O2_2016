# logistical models for Parnell's thesis on H2O2-LOH effect
rm(list=ls())
setwd("~/~/github/LOH_H2O2_2016/model-diagram")

#these are my old function for CLS-H2O2
logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }
logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }

### for Cb, Cv figure
Tv = 0.042
Tb = 0.013
t = seq(0, 0.3, by=0.001)
w = 2;

pdf( 'diagram-CbCv.pdf',width=6,height=5)
#pdf( 'diagram-compareCvCv-TgTc.pdf')

#mat = matrix( seq(1,2, nrow=2, ncol= 1 ); 
#layout(mat);

par(mai=c(1,1,0.5,1))
v = logistical.viability( Tv, w, t)
plot( v ~ t, type='l', col='blue', log='x', xlim=c(0.001,0.5), xlab='H2O2 fraction', ylab='viability', lwd=4)
points( Tv, 0.5, pch=19, col="red", cex=1.2 );
arrows( Tv, 0.5, Tv, -1, lty=2, col="red", lwd=3);
mtext( "Cv",side=1,at=c(Tv) );
mtext( 'black percentage', side=4, line=3 ) #2012 July 28 change

par(new=T)
b.max = 0.04; b.min = 0.002;
b = logistical.black( b.max, b.min, Tb, w, t)
plot( b ~ t, log='x', axes=F, xlab='', ylab='', type='l', lwd=4, ylim=c(b.min, 1.5*b.max))
axis(4, pretty(range(b)))
points( Tb, (b.max+b.min)/2, pch=19, col="green", cex=1.2 );
arrows( Tb, (b.max+b.min)/2, Tb, -1, lty=2, col="green", lwd=3);
mtext( "Cb",side=1,at=c(Tb) );

dev.off()

### for Tg, Tc figures

Tg = 5.5
Tc = 5
t = seq(0, 10, by=0.01)
w = 5;
                        
########
pdf( 'diagram-TgTc.pdf',width=6,height=5)
par(mai=c(1,1,0.5,1))
v = logistical.viability( Tc, w, t)
par(mai=c(1,1,0.5,1))              
plot( v ~ t, type='l', col='blue', log='', xlim=c(01,10), xlab='Days', ylab='viability', lwd=4);           
              
points( Tc, 0.5, pch=19, col="red", cex=1.5 );
arrows( Tc, 0.5, Tc, -1, lty=2, col="red", lwd=3);
mtext( "Tc",side=1,at=c(Tc) );
mtext( "black percentage", side=4, line=3)
par(new=T)
b.max = 0.04; b.min = 0.002;
b = logistical.black( b.max, b.min, Tg, w, t)
plot( b ~ t, log='', axes=F, xlab='', ylab='', type='l', lwd=4, ylim=c(b.min, 1.5*b.max))
axis(4, pretty(range(b)))
points( Tg, (b.max+b.min)/2, pch=19, col="green", cex=1.2, lwd=4 );
arrows( Tg, (b.max+b.min)/2, Tg, -1, lty=2, col="green", lwd=3);
mtext( "Tg",side=1,at=c(Tg) );
dev.off()

pdf( 'diagram-CbTg.pdf',width=6,height=5)
Tv = 0.042
Tb = 0.013
t = seq(0, 0.3, by=0.001)
w = 2;

par(mai=c(1,1,0.5,1))
v = logistical.viability( Tv, w, t)
plot( v ~ t, type='l', col='blue', xaxt = "n", log='x', xlim=c(0.001,0.5), xlab='H2O2 fraction OR days', ylab='viability', lwd=4)

par(new=T)          
Tg = 5.5
Tc = 5
w = 5;  

par(new=T)
b.max = 0.04; b.min = 0.002;
b = logistical.black( b.max, b.min, Tb, w, t)
plot( b ~ t, log='x', axes=F, xlab='', ylab='', yaxt= "n",xaxt = "n",type='l', lwd=4, ylim=c(b.min, 1.5*b.max))
axis(4, pretty(range(b)))
points( Tb, (b.max+b.min)/2, pch=19, col="green", cex=1.2 );
arrows( Tb, (b.max+b.min)/2, Tb, -1, lty=2, col="green", lwd=3);
mtext( "Cb",side=1,at=c(Tb) );

Tg = 6.9
Tc = 5
w = 5;
b.max = 0.04; b.min = 0.002;
t = seq(0, 10, by=0.01)
par(new=T)
b = logistical.black( b.max, b.min, Tg, w, t)
plot(b~t, log='', col="magenta", xlab='', ylab='',yaxt="n", xaxt ="n" ,type='l', lwd=4, ylim=c(b.min, 1.5*b.max))
axis(4, pretty(range(b)))
points( Tg, (b.max+b.min)/2, pch=19, col="green", cex=1.2, lwd=4 );
arrows( Tg, (b.max+b.min)/2, Tg, -1, lty=2, col="green", lwd=3);
mtext( "Tg",side=1,at=c(Tg) );
legend("topright",col=c("magenta","black","blue"),lty=1,legend=c("Days","H2O2","Viability"))
dev.off()

