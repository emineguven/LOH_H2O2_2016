#081207 r = - ( - 2 dPb/dt ) 1/g
# modified from  _050307.plot.YPS163.030906.R
#050307 exp.id added

rm(list=ls());

library(nlme)

file = "M5.032906.tab"; 

exp.id = "M5.032906";  

logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }

logistical.mortality <- function( T, w, t ) {  ( w * t^(w-1) ) / (( 1 + (t/T)^w ) * T^w ); }

logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }

derivative.black <- function(b.max, b.min, T, w, t) { (b.max - b.min) * w * (t^(w -1)) / ((1 + (t / T )^w)^2 * T^w ); }

genome.integrity <- function(b.max, b.min, T, w, t) { 2 * (b.max - b.min) / (1 + (t/T)^w) + (1 - 2 * b.max); }


tb= read.table( file, header=T, sep="\t", fill=T);
tb$t = tb$t / 24;
tb2 = tb;
labels = names( tb );

#tb2$B0.5[22]=1;  ###########a upper limit for B0.5

### remove last point, regrowth? 
#tb2 = tb2[1:27, ]

#normalize tb2
 tmp.indice = seq( 1 : length(labels) );
 names( tmp.indice ) = labels;
 for( j in tmp.indice["white"] : length(labels) ) {
 #for( j in 5:13 ) {
   for( i in  1: length(tb2[,1]) ) {
     if ( is.na(tb2[i,j] ) ) { tb2[i,j] = 0; }
     tb2[i,j] = tb2[i,j] * tb2[i,2] * tb2[i,3];
   }
 }

#generate row indice for averaging
 # row.num = c( 2, 3, 5, 6, 8, 9:15 );
 row.steps = as.vector( table( tb2$t ) ); # a new trick, ha

#generate tb.m
 n.row = length( row.steps );
 col.labels = c( labels[ c(1,tmp.indice["white"] : length(labels)) ], "total" ); ######## bug here 091806
 #col.end = 11;

 tb.m = data.frame( matrix( nrow=n.row, ncol= length(col.labels) ) )  # mean values
 names( tb.m ) = col.labels; 

 upper.row = 0; #up pointer
 lower.row = 0; #low pointer 
 for( i in 1:n.row ) {
   upper.row = lower.row + 1;
   lower.row = upper.row + row.steps[ i ] - 1 ;

   tb.m[i,1]  = tb2$t[ upper.row ]
   for( j in 2: ( length(col.labels) - 1 ) ) {
     tb.m[i, j] = mean( tb2[ upper.row : lower.row, j+3] )
   }
   tb.m$total[i] = sum( tb.m[i, 2:( length(col.labels) - 1 )], na.rm=T );
 }

 ### get the standard errors of whites and blacks #################### 091306 change
 upper.row = 0; #up pointer
 lower.row = 0; #low pointer 
 sd.w    =  numeric( n.row );
 sd.b    = numeric( n.row );
 sd.b05 = numeric( n.row);
 for( i in 1:n.row ) { 
   upper.row = lower.row + 1;
   lower.row = upper.row + row.steps[ i ] - 1 ;
   if ( ( lower.row - upper.row ) > 0 ) {
     sd.w[i]    = sd( tb2$white[ upper.row : lower.row] ) ;
     sd.b[i]    = sd( tb2$black[ upper.row : lower.row] ) ;
     sd.b05[i]  = sd( tb2$B0.5[ upper.row : lower.row] ) ;
   } 
 }
 tb.m$sd.w = sd.w
 tb.m$sd.b = sd.b
 tb.m$sd.b05 = sd.b05


# output to out,
# columns in out
 #old# header = c("t","half.over.black","Pb","Rb","R0.5", "R0.75", "s", "g" );
 header = c("t","Pb","s", "g","m","Rb","R0.5", "L"  );   ### L = rate(1/2) / rate(black) 
 out = data.frame( matrix( nrow= length(tb.m[,1]) , ncol= length(header) ) );
 names( out ) = header;

 out$t = tb.m$t; # "t"

#######calculate s , g
 out$s      = tb.m$total / tb.m$total[1]
 out$Pb     = tb.m$black / tb.m$total;
 out$g      = 1 - 2 * out$Pb;
 out$e.s    =  out$s * (tb.m$sd.w / tb.m$white);
 out$e.b    =  out$Pb * (tb.m$sd.b / tb.m$black);

 out$R0.5.raw   = tb.m$B0.5 /tb.m$total; #needs to be adjusted by g.e
 out$e.b05  =  out$R0.5.raw * ( tb.m$sd.b05 / tb.m$B0.5)

 ### out$e.b05[c(13,14)] = c(0.005, 0.003);            ####030707 change

 ### plot of the raw data 
 plot( out$s ~ out$t , type='l', main= file, col="blue");
 lines( out$Pb ~ out$t, col="black");
 labels = c("viability","black");
 ltypes = c(1,1);
 legend( (max(out$t)*0.7 ), 0.8, labels, col=c("blue", "black"), lty = ltypes);

####### remove outliers#######################
 #out$t[c(10:12)] = NA;
 #out = out[ (! is.na(out$t) ), ]

######## calclulate T0.5 w_s
 t= out$t
 s= out$s
 fm.s <- gnls( s ~ 1 /( 1 + ( t / v )^ w ) , start = list( v = 5, w = 11 )  );
 fm.s # this the half life T1/2 w_s

 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 fit.m = logistical.mortality ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( out$s ~ out$t );
 lines( fit.s ~ t, col="blue");
 par(new=T)
 plot( fit.m ~ t, col="brown", axe=F);
 
 #estimate errors for s
 error.s = out$s -  logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], out$t );
 error.s = abs(error.s)
 out$e.s = ifelse( out$e.s==0, error.s, out$e.s );
 out$e.s= ifelse( out$e.s<0.01, 0.01, out$e.s );
 
######## calclulate T_g and w_g using Pb
 # Pb = b.max - (b.max-b.min) / (1 + (t/T.g)^w)
 b = out$Pb
 t = out$t;
 b.min = min( b ); b.min  # 0.00705687
 b.max = max( b ); b.max  # 0.3482697

 b.min = 0.002; 
 b.max = 0.31; 
 w.g = 12;

 #my.fun = function (t, T.g) { b.max - (b.max - b.min) /( 1 + ( t / T.g )^ w.g ); }
 my.fun = function (t,T.g, w.g) { b.max - (b.max - b.min) /( 1 + ( t / T.g )^ w.g ); }

 # ws = 1 / out$e.b ^ 2
 #####   Pr[B>=b] = b.max - (b.max - b.min) /( 1 + ( t / T )^ w ) ########formula for cumulative Pr(b)
 #fm.b = gnls( b ~ 0.0633694 - (0.0633694 - 0.008868349) / (1 + (t/T.g)^w ), start=list( T.g=5, w=3), weights = ws);
 #fm.b = gnls( b ~ my.fun(t, T.g), start=list( T.g=5) );
 fm.b = gnls( b ~ my.fun(t, T.g, w.g), start=list( T.g=5, w.g=12) );              
 fm.b

 T.g = fm.b$coefficients[1];  
 w.g = fm.b$coefficients[2];  #?????????

 #estimate error for black
 error.b = out$Pb - logistical.black( b.max, b.min, T.g, w.g, out$t);
 error.b = abs( error.b );
 out$e.b = ifelse( out$e.b==0, error.b, out$e.b );
 out$e.b = ifelse( is.na(out$e.b), error.b, out$e.b );
 		# out$e.b[c(13,14)] = c(0.1, 0.05);  ##???

 ### logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }
 t = seq(0, max(out$t),by=0.1);
 fit.b = logistical.black( b.max, b.min, T.g, w.g, t);
 
 plot( out$Pb ~ out$t );
 lines( fit.b ~ t, col="blue");

###################overlay s, Pb, 
 pdf("080707.032906M5.s.Pblog.pdf", width=6, height=6)
 par(mar=c(5,4,4,4)+0.1);

 xlim = c(0,17)

 plot( out$s ~ out$t, col="blue", xlab="t (days)",ylab="Viability", ylim=c(1E-3, 1.2), pch=16, xlim=xlim );
 #title (file);                                        ############change here
 arrows( out$t, (out$s - out$e.s), out$t, (out$s + out$e.s), length=0.1, angle=90,code=3, lty=1, col="blue" );

 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 lines( fit.s ~ t, col="blue",lty=2);
 
 cols  =c("blue","black" );
 labels=c("viability","b(t)")
 ltypes=c(2,2) 
 pch   =c(16,15)
 legend( 10 , 0.5 ,labels, col=cols, lty=ltypes, pch=pch);
 text( 12, 0.55, exp.id);
 
 T.c = fm.s$coefficients[1]
 points( T.c, 0.5, pch=19, col="red", cex=1.2 );
 arrows( T.c, 0.5, T.c, -1, lty=2, col="red");
 mtext( "Tc",side=1,at=c(T.c) );


 par(new=T)
 fit.b = logistical.black( b.max, b.min, T.g, w.g, t);
 ylim = c( 5E-4, b.max*1.5 ) 
 #points( out$t, out$Pb, pch=15, xlab="",ylab="", ylim=ylim );
 plot( out$Pb ~ out$t, xlab="",ylab="",ylim=ylim, pch=15, xlim=xlim, log='y',axes=F );
 arrows( out$t, (out$Pb - out$e.b), out$t, (out$Pb + out$e.b), length=0.1, angle=90,code=3, lty=1 );
 lines( fit.b ~ t, lty=2);
 axis(4, at=c( 1E-3, 1E-2, 1E-1, 0.5, 1) )
 mtext( "b(t)", 4, 2);

 points ( T.g,  ( b.max/2 + b.min/2), pch=15, col="red", cex=1.2);
 arrows( T.g, (b.max/2 + b.min/2), T.g, 1E-8, lty=2, col="red" );
 mtext( "Tg",side=1,at=c(T.g) );

 mtext( "Tg",side=1,at=c(T.g) );

 dev.off();



###################overlay s, m, Pb, Rb
 pdf("081207.032906M5.r.m.pdf", width=6, height=6)
 par(mar=c(5,4,4,4)+0.1);
 
 xlim = c(0,17)
 t = seq(0, max(out$t),by=0.1);

 fit.dPb = derivative.black(b.max, b.min,  T.g, w.g, t );  #first point is Inf
 fit.g   = genome.integrity(b.max, b.min,  T.g, w.g, t );
 fit.Rb  = fit.dPb / fit.g;                        ### rate of becoming blacks
 fit.r = 2 * fit.Rb;  			###081207 r == genomic instability rate
 plot( fit.r ~ t, col="black", xlab = "t (days)", ylab="Instability rate", type='l', xlim=xlim);
 # title (file);

 par(new=T)
 fit.m = logistical.mortality ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( fit.m ~ t, col="red", axe=F, xlab="",ylab="", type='l', xlim=xlim);
 axis(4, at=pretty(fit.m) )
 mtext( "Mortality rate", 4, 2);

 cols   =c( "black", "red", "blue")
 labels =c("instability rate", "mortality rate", "viability")
 ltypes =c(1,1,2,1) 
 pch    =c( NA, NA,NA,16,15,)
 legend( 7 , max(fit.m)/4, labels, col=cols, lty=ltypes, pch=pch);
 text( 9,  max(fit.m)/4 + 0.1, exp.id);

 par(new=T);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( fit.s ~ t, col="blue",lty=2, axes=F, xlab="",ylab="", type='l', xlim=xlim);
 #axis(4, at=pretty(range(fit.s)))

 #par(new=T)
 #fit.b = logistical.black( b.max, b.min, fm.b$coefficients[1], fm.b$coefficients[2], t);
 #points( out$Pb ~ out$t, pch=15 );
 #arrows( out$t, (out$Pb - out$e.b), out$t, (out$Pb + out$e.b), length=0.1, angle=90,code=3, lty=1 );
 #lines( fit.b ~ t, lty=2);
 #lines( fit.g ~ t, col="green")

 dev.off();



################################## tmp codes
 ##calculate Rb = d(Pb)/dt
 # derivative.black <- function(b.max, b.min, T, w, t) {  (b.max - b.min) * w * t ^(w -1) / (1 + (t / T ))^2; }
 #out$dPb = derivative.black(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], out$t );  #first point is Inf
 #out$Rb     = out$dPb / out$g.e;                        ### rate of becoming blacks

 #par(new=T);
 #fit.dPb = derivative.black(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], t );  #first point is Inf
 #fit.g   = genome.integrity(b.max, b.min,  fm.b$coefficients[1], fm.b$coefficients[2], t );
 #fit.Rb  = fit.dPb / fit.g;                        ### rate of becoming blacks
 #plot( fit.Rb ~ t, col="black", xlab ="", ylab="", type='l', axes=F, lty=2);
 #axis(4, at=pretty(fit.Rb) )
 #mtext( "black rate", 4, 2);



##################################
 pdf("100207.032906M5.b0.5log.pdf", width=6, height=6)
 par(mar=c(5,4,4,4)+0.1);

 out$e.b05[7] = out$e.b05[6] /2 + out$e.b05[8] /2 

 out$g.e = genome.integrity(b.max, b.min,  T.g, w.g, out$t );
 out$R0.5   = out$R0.5.raw / out$g.e;                   #### 090706 change
 ylim = c(1E-5, 0.1)
 plot( out$R0.5 ~ out$t, col="green", xlab="t (days)",ylab="b1/2(t)", type='l',lty=2, ylim=ylim,xlim=xlim,log='y' );
 # title (file);                           ############change here
 points( out$R0.5 ~ out$t, col="green", pch=16);
 # arrows( out$t, (out$R0.5 - out$e.b05), out$t, (out$R0.5 + out$e.b05), length=0.05, angle=90,code=3, lty=1, col="blue" );
 ylow = (out$R0.5 - out$e.b05);
 ylow[ ylow<=0] = 1E-5; 
 yup = (out$R0.5 + out$e.b05);
 arrows( out$t, ylow, out$t, yup, length=0.05, angle=90,code=3, lty=1, col="green" );

 par(new=T);
 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( fit.s ~ t, col="blue", type='l',lty=2, axes=F, xlab="",ylab="", xlim=xlim );
 axis(4, at=pretty(range(fit.s)))
 mtext( "Viability", 4, 2);

 #add labels here
 cols  =c("blue","green");
 labels=c("viability", "b1/2(t)")
 ltypes=c(2,2,1,1) 
 pch   =c( NA,16,NA, 16)
 legend( 12 , 0.3  ,labels, col=cols, lty=ltypes, pch=pch);
 text( 14, 0.35, exp.id); 
 dev.off();


 ##################################
 pdf("081007.032906M5.L.pdf", width=6, height=6)
 par(mar=c(5,4,4,4)+0.1);
 
 out$L = out$R0.5 / out$Pb ; 
 out$e.L = out$L * sqrt( (out$e.b05 / out$R0.5)^2 + (out$e.b / out$Pb)^2)
 #ylim = c(1E-3, max(out$L[2:length(out$L)])*1.7 )
 ylim =c(1E-2, 5);

 #plot( out$L ~ out$t, col="red", xlab="t(days)", ylab="L(t)",  pch=16, type='l', xlim=xlim,ylim=c(0,1) ); 
 plot( out$L ~ out$t, col="red", xlab="t(days)", ylab="L(t)",  pch=16, type='l', xlim=xlim,ylim=ylim, log='' ); 
 points( out$t, out$L, pch=16, col="red" );
 arrows( out$t, (out$L - out$e.L), out$t, (out$L + out$e.L), length=0.1, angle=90,code=3, lty=1, col="red" );

 par(new=T);
 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( fit.s ~ t, col="blue", type='l',lty=2, axes=F, xlab="",ylab="", xlim=xlim  );
 axis(4, at=pretty(range(fit.s)))
 mtext( "Viability", 4, 2);

 cols  =c("blue","red");
 labels=c("viability","L(t)")
 ltypes=c(2,1,1,1) 
 pch   =c( NA,NA,NA, 16)
 legend( 10 , 0.8  ,labels, col=cols, lty=ltypes, pch=pch);
 text( 11.5, 0.85, exp.id );
 dev.off();


 ##################################
 pdf("072207.032906M5.Fb.pdf", width=6, height=6)
 par(mar=c(5,4,4,4)+0.1);
 
 e.fb = tb.m$sd.b; 
 e.fb[6] = e.fb[6]*0.8;

 ylim = c(min(tb.m$black)/2, max(tb.m$black)*2 );  ##050307
 plot( tb.m$black ~ tb.m$t, col="black", type='l',lty=1, xlab="t (days)",ylab="Concentration of blacks", xlim=xlim, ylim=ylim);
 arrows( tb.m$t, (tb.m$black - e.fb), tb.m$t, (tb.m$black + e.fb), length=0.1, angle=90,code=3, lty=1, col="red" );
 points( tb.m$t, tb.m$black, pch=16);

 par(new=T);
 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( fit.s ~ t, col="blue", type='l',lty=2, axes=F, xlab="",ylab="", xlim=xlim  );
 axis(4, at=pretty(range(fit.s)))
 mtext( "Viability", 4, 2);

 cols  =c("blue","black");
 labels=c("viability","Concen. black")
 ltypes=c(2,1,1,1) 
 pch   =c( NA,16,NA, 16)
 legend( 10 , 0.8  ,labels, col=cols, lty=ltypes, pch=pch);
 text( 11.5, 0.85, exp.id); 

 dev.off();



##################################
 pdf("080707.032906M5.Fb.F0.5.pdf", width=6, height=6)
 par(mar=c(5,4,4,4)+0.1);
 
 e.fb = tb.m$sd.b; 
 e.fb[2] = e.fb[2]*0.6;
 y.scale = 1E3;

 #ylim = c(min(tb.m$black / y.scale)/2, max(tb.m$black / y.scale)*2 );  ##050307
 ylim = c( 0.1, max(tb.m$black / y.scale)*2 );  ##080707

 plot( tb.m$black / y.scale ~ tb.m$t, col="black", type='l',lty=1, xlab="t (days)",ylab="Cells per ml / 10^3", xlim=xlim, ylim=ylim, log='');
 arrows( tb.m$t, (tb.m$black - e.fb)/y.scale, tb.m$t, (tb.m$black + e.fb)/y.scale, length=0.1, angle=90,code=3, lty=1, col="red" );
 points( tb.m$t, tb.m$black/y.scale, pch=16);

 par(new=T);
 e.05 = tb.m$sd.b05; 
 e.05[7] = e.05[6] /2 + e.05[8] /2 

 #ylim = c(min(tb.m$B0.5 / y.scale)/2, max(tb.m$B0.5 / y.scale)*1.5 ); 
 plot( tb.m$B0.5/y.scale ~ tb.m$t, col="green", type='l',lty=1, axes=F, xlim=xlim, ylim=ylim, xlab="",ylab="",log='');
 arrows( tb.m$t, (tb.m$B0.5 - e.05)/y.scale, tb.m$t, (tb.m$B0.5 + e.05)/y.scale, length=0.1, angle=90,code=3, lty=1, col="blue" );
 points( tb.m$t, tb.m$B0.5/y.scale, pch=16, col="green");
 # axis(4, at = pretty(range(ylim)))
 # mtext( "Half-blacks per ml/ 10^3", 4, 2);

 par(new=T);
 t = seq(0, max(out$t),by=0.1);
 fit.s = logistical.viability ( fm.s$coefficients[1], fm.s$coefficients[2], t );
 plot( fit.s ~ t, col="blue", type='l',lty=2, axes=F, xlab="",ylab="", xlim=xlim  );

 cols  =c("blue","black", "green");
 labels=c("viability","blacks per ml", "1/2 blacks per ml")
 ltypes=c(2,1,1,1) 
 pch   =c( NA,16,16,16)

 legend( 9 , 0.8  ,labels, col=cols, lty=ltypes, pch=pch);
 text( 11.5, 0.85, exp.id); 

 dev.off();


quit("yes");



###this did not work out
 L = out$L[2:length(out$L)]; t=out$t[2:length(out$L)];
 L = c( 0, L);  t = c(0, t);
 fm.L = gnls( L ~ A * exp(B*t) + C, start=list( A=0.005, B=0.5, C=0.005) );
 t = seq(0, max(out$t),by=0.1);
 fit.L = fm.L$coefficients[1] * exp( fm.L$coefficients[2] * t) + fm.L$coefficients[3]
 lines( fit.L ~ t, col="red", lty=2);

 # out$e.b05  =  ifelse( out$e.b05==0, out$R0.5 * (out$e.b/out$Pb), out$e.b05    );
 # out$e.dPb =  out$dPb * (out$e.b / out$Pb)  ## I use out$e.b as proxy  
 