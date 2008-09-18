`curvetest` <-
function(XX1, YY, XX2=NULL, ZZ=NULL, 
  kernel=c("Trio",  "Gaussian", "Uniform", "Triweight", "Triangle", "Epanechnikov", "Quartic"), 
  equal.var=TRUE, hh=0.5, plotit=FALSE, conf.level=0.05, nn=100)
{   
    if(missing(kernel)) kernel="Trio"
    pm=pmatch(kernel, c("Trio", "Gaussian", "Uniform", "Triweight", "Triangle", "Epanechnikov", "Quartic")) 
    if(is.na(pm) |pm==1) ww<-function(x) {(1-abs(x)^3)^3*(abs(x)<1)}  else
    if(pm==2)  ww<-function(x) { exp(-x^2/2)/sqrt(2*pi)}  else
    if(pm==3)  ww<-function(x) { (abs(x)<=1)/2}  else
    if(pm==4)  ww<-function(x) { (1-x^2) ^3 *35/32*(abs(x)<=1) } 
    if(pm==5)  ww<-function(x) {(1-abs(x))(abs(x)<=1) }  
    if(pm==6)  ww<-function(x) { (1-x^2)   *3/4*(abs(x)<=1) }  
    if(pm==7)  ww<-function(x) { (1-x^2) ^2 *15/16*(abs(x)<=1) }   
    distance<-function(x,y) {sqrt(sum((x-y)^2))};
    app1<-function(ij, t1,t2,t3){ww((t1[ij[1]]-t2[ij[2]])/t3)}
    app2<-function(ij, t1,t2,t3,t4) {t4[ij[1],ij[2]]*sum( t4[ij[1],]*(t1[ij[1]]-t2)*(t2[ij[2]]-t2))} 
    n1<-length(XX1);
    if(n1!=length(YY)) stop("Lengths  of Response and Predictor  are not equal.")
    xmax<-max(XX1); 
    xmin<-min(XX1)
    if(!is.null(XX2)&!is.null(ZZ)){
        n2<-length(XX2);
        if(n2!=length(ZZ)) stop("Lengths  of Response and Predictor  are not equal.")
        xmax<-max(xmax,XX2); 
        xmin<-min(xmin,XX2)
    }
   myx<-seq(xmin, xmax, length=nn)
   if(is.null(XX2)) {
      ##Compute k0; 
      exnnn1<-expand.grid(1:nn, 1:n1)
      temp<-apply(exnnn1, 1, FUN=function(ij, t1=myx, t2=XX1,t3=hh) app1(ij,t1,t2,t3))
      w1ij<- matrix(temp,nrow=nn,ncol=n1)
      temp<-apply(exnnn1, 1, FUN=function(ij, t1=myx, t2=XX1,t3=hh,t4=w1ij) app2(ij, t1,t2,t3,t4))
      sx<-matrix(temp, nrow=nn,ncol=n1)
      sx<-sweep(sx, 1, rowSums(sx),FUN="/");#function sx in f1_hat=s(x)*YY

      TX<-sweep(sx, 1, sqrt(rowSums(sx^2)),FUN="/")
      k0<-sum( sqrt(sum((TX[1:(nn-1),]-TX[2:nn,])^2)) )
      exn1n1<-expand.grid(1:n1,1:n1)
      temp<-apply(exn1n1,1,FUN=function(ij,t1=XX1,t2=XX1,t3=hh) app1(ij,t1,t2,t3))
      L1<-matrix(temp,nrow=n1,ncol=n1)  #The matrix L in Y_hat=LY.
      L1<-sweep(L1,1, rowSums(L1),FUN="/")
      AA<-diag(rep(1,n1))-L1; delta11<-sum((AA)^2);#delta11<-trace((I-L1)*(I-L1)').
      delta12<- sum((AA%*%AA)^2)   #delta12<-trace(AA*t(AA))^2.
      vv<-delta11^2/delta12        #Estimated degree of freedom
      eYY<-L1%*%YY;                ###Estimated YY<-f1(x) 
      resid1<-YY-eYY;
      sigma.square<- sum(resid1^2)/delta11 ##Estimated sigma^2
      CC <- max(abs( sx%*%YY)/( sqrt(sigma.square)*sqrt(rowSums(sx^2))))  ##Test statistics of H0:f1=0.
      ff <- k0/3.1415926 * (1+(CC*CC)/vv)^(-vv/2)+pt(-CC, vv)*2
      if(ff>1) ff=1
      if(plotit){
        xn=deparse(substitute(XX1))
        yn=deparse(substitute(YY))          
        plot(XX1,YY, xlab=xn, ylab=yn,col=2);
        lines(myx, sx%*%YY,col=2,lwd=3, lty=2)
      }
  }
 else if(!is.null(XX2) & !is.null( ZZ)& equal.var){
      ##Compute k0; 
      exnnn1<-expand.grid(1:nn, 1:n1)
      temp<-apply(exnnn1, 1, FUN=function(ij, t1=myx, t2=XX1,t3=hh) app1(ij,t1,t2,t3))
      w1ij<- matrix(temp,nrow=nn,ncol=n1)
      temp<-apply(exnnn1, 1, FUN=function(ij, t1=myx, t2=XX1,t3=hh,t4=w1ij) app2(ij, t1,t2,t3,t4))
      sx<-matrix(temp, nrow=nn,ncol=n1)
      sx<-sweep(sx, 1, rowSums(sx),FUN="/");#function sx in f1_hat=s(x)*YY

      exnnn2<-expand.grid(1:nn, 1:n2)
      temp<-apply(exnnn2, 1,  FUN=function(ij, t1=myx, t2=XX2,t3=hh) app1(ij,t1,t2,t3))
      w2ij<- matrix(temp,nrow=nn,ncol=n2)
      temp<-apply(exnnn2, 1,FUN=function(ij, t1=myx, t2=XX2,t3=hh,t4=w2ij) app2(ij,t1,t2,t3,t4))
      tx<-matrix(temp,nrow=nn,ncol=n2)
      tx<-sweep(tx, 1, rowSums(tx),FUN="/")   ###function tx in f2_hat=t(x)*ZZ
      stx<-sqrt(rowSums(sx^2)+rowSums(tx^2))   ###sqrt(||sx||^2+||tx||^2)
      T1X <- sweep(sx, 1,stx,FUN="/");  
      T2X <- sweep(tx, 1,stx,FUN="/")
      k0<-0;
      for(ii in 2:nn)
           k0 <- k0+sum(sqrt(distance(T1X[ii,],T1X[ii-1,])^2+distance(T2X[ii,],T2X[ii-1,])^2))
     #######Calculate the degree of freedoms for the two models.
     exn1n1<- expand.grid(1:n1, 1:n1)
     temp<- apply(exn1n1, 1,  FUN=function(ij, t1=XX1,t2=XX1,t3=hh) app1(ij,t1,t2,t3))
     wxx1<-matrix(temp,nrow=n1,ncol=n1)
     temp<-apply(exn1n1, 1, FUN=function(ij, t1=XX1, t2=XX1,t3=hh,t4=wxx1) app2(ij,t1,t2,t3,t4))
     L1<- matrix(temp, nrow=n1,ncol=n1)  #The matrix L in Y_hat=LY.
     L1<-sweep(L1,1, rowSums(L1),FUN="/")
     exn2n2<-expand.grid(1:n2, 1:n2)
     temp<-apply(exn2n2, 1,  FUN=function(ij, t1=XX2, t2=XX2,t3=hh) app1(ij,t1,t2,t3))
     wxx2<-matrix(temp, nrow=n2,ncol=n2)
     temp<-apply(exn2n2, 1,FUN=function(ij,t1=XX2,t2=XX2,t3=hh,t4=wxx2)app2(ij,t1,t2,t3,t4))
     L2<-matrix(temp,  nrow=n2,ncol=n2)  #The matrix L in Y_hat=LY.
     L2<-sweep(L2,1, rowSums(L2),FUN="/")
     AA<-diag(rep(1,n1))-L1; delta11<-sum((AA)^2);#delta11<-trace((I-L1)*(I-L1)').
     delta12<- sum((AA%*%AA)^2)   #delta12<-trace(AA*t(AA))^2.
     v1<-delta11^2/delta12  ##Degree of freedom of model 1, see paper of Cleveland.
     AA<-diag(rep(1,n2))-L2;  
     delta21<-sum((AA)^2);
     delta22<- sum((AA%*%AA)^2)
     v2<-delta21^2/delta22
     vv<-(delta11+delta21)^2/(delta12+delta22)
   ####################################Fit Study
     eYY<-L1%*%YY;
     eZZ<-L2%*%ZZ; ###Estimated YY<-f1(x) and ZZ<-f2(x)
     resid1<-YY-eYY;
     resid2<-ZZ-eZZ
     sigma.square<-(sum(resid1^2)+sum(resid2^2))/(delta11+delta21)##Estimated sigma^2
     CC <- max(abs( sx%*%YY-tx%*%ZZ)/( sqrt(sigma.square)* stx))  ##Test statistics of H0:f1=f2.
     ff <- k0/3.1415926*(1 + (CC*CC)/vv)^(-vv/2)+pt(-CC,vv)*2
     if(ff>1) ff=1
     if(plotit) {        
        xn=paste(deparse(substitute(XX1)), "or",  deparse(substitute(XX2)))
        yn=paste(deparse(substitute(YY)), "or",  deparse(substitute(ZZ)))          
        plot(XX1,YY,xlab=xn,ylab=yn, col=2, xlim=c(min(myx), max(myx)), 
              ylim=c(min(YY,ZZ),max (YY,ZZ)));
        lines(myx, sx%*%YY, lwd=4,col=2, lty=2);
        points(XX2,ZZ, col=4, pch=2);
        lines(myx, tx%*%ZZ,col=4, lwd=4, lty=4)
     }
    }
    if(!is.null(XX2) & !is.null( ZZ)& !equal.var){
          ##Compute k0; 
          exnnn1<-expand.grid(1:nn, 1:n1)
          temp<-apply(exnnn1, 1, FUN=function(ij, t1=myx, t2=XX1,t3=hh) app1(ij,t1,t2,t3))
          w1ij<- matrix(temp,nrow=nn,ncol=n1)
          temp<-apply(exnnn1, 1, FUN=function(ij, t1=myx, t2=XX1,t3=hh,t4=w1ij) app2(ij, t1,t2,t3,t4))
          sx<-matrix(temp, nrow=nn,ncol=n1)
          sx<-sweep(sx, 1, rowSums(sx),FUN="/");#function sx in f1_hat=s(x)*YY
          #######Calculate the degree of freedoms for the two models.
          exnn<-expand.grid(1:n1,1:n1)
          temp<-apply(exnn,1,FUN=function(ij,t1=XX1,t2=XX1,t3=hh) app1(ij,t1,t2,t3))
          wxx1<-matrix(temp ,nrow=n1,ncol=n1)
          temp<-apply(exnn,1,FUN=function(ij,t1=XX1,t2=XX1,t3=hh,t4=wxx1) app2(ij,t1,t2,t3,t4))
          L1<-matrix(temp, nrow=n1,ncol=n1)  #The matrix L in Y_hat=LY.
          L1<-sweep(L1,1, rowSums(L1),FUN="/")
          exnn<-expand.grid(1:n2, 1:n2)
          temp<-apply(exnn, 1, FUN=function(ij, t1=XX2, t2=XX2,t3=hh) app1(ij,t1,t2,t3))
          wxx2<-matrix(temp,nrow=n2,ncol=n2)
          temp<-apply(exnn, 1, FUN=function(ij, t1=XX2, t2=XX2,t3=hh,t4=wxx2) app2(ij,t1,t2,t3,t4))
          L2<- matrix(temp, nrow=n2,ncol=n2)  #The matrix L in Y_hat=LY.
          L2<-sweep(L2,1, rowSums(L2),FUN="/")
          AA<-diag(rep(1,n1))-L1; delta11<-sum((AA)^2);#delta11<-trace((I-L1)*(I-L1)').
          delta12<- sum((AA%*%AA)^2)   #delta12<-trace(AA*t(AA))^2.
          v1<-delta11^2/delta12  ##Degree of freedom of model 1, see paper of Cleveland.
          AA<-diag(rep(1,n2))-L2;  
          delta21<-sum((AA)^2); 
          delta22<- sum((AA%*%AA)^2)
          v2<-delta21^2/delta22
          vv <- (n1+n2)^2*v1*v2/(n2^2*v1 + n1^2* v2)#Total degree of freedom
          ####################################Fit Study
          eYY<-L1%*%YY;eZZ<-L2%*%ZZ; ###Estimated YY<-f1(x) and ZZ<-f2(x)
          resid1<-YY-eYY;resid2<-ZZ-eZZ  ##The residuals.
          esigma1<-sqrt((sum(resid1^2))/delta11 )    #Estimated sigma_1^2
          esigma2<-sqrt(sum(resid2^2) /delta21 )     ##Estimated sigma2^2
          exnn<-expand.grid(1:nn,1:n2)
          temp<-apply(exnn, 1, FUN=function(ij, t1=myx, t2=XX2,t3=hh) app1(ij,t1,t2,t3))
          w2ij<- matrix(temp,    nrow=nn,ncol=n2)
          temp<- apply(exnn,1,FUN=function(ij,t1=myx,t2=XX2,t3=hh,t4=w2ij) app2(ij,t1,t2,t3,t4) )
          tx<-matrix(temp, nrow=nn, ncol=n2)
          tx<-sweep(tx, 1, rowSums(tx),FUN="/")   ###function tx in f2_hat=t(x)*ZZ
          stx<-sqrt(esigma1^2*rowSums(sx^2)+esigma2^2*rowSums(tx^2)) ###sqrt(\sigma_1*||sx||^2+\sigma2*||tx||^2)
          T1X <- esigma1*sweep(sx, 1,stx,FUN="/");
          T2X <- esigma2*sweep(tx, 1,stx,FUN="/")
          k0<-0; 
          for(ii in 2:nn) 
          k0 <- k0+sum(sqrt(distance(T1X[ii,],T1X[ii-1,])^2+distance(T2X[ii,],T2X[ii-1,])^2))
          CC <- max(abs( sx%*%YY-tx%*%ZZ)/stx)  ##Test statistics of H0:f1=f2.
          ff<- k0/3.1415926 *  exp(-CC^2/2)+ pnorm( - CC) * 2
          if(ff>1) ff=1
      if(plotit){
        xn=paste(deparse(substitute(XX1)), "or",  deparse(substitute(XX2)))
        yn=paste(deparse(substitute(YY)), "or",  deparse(substitute(ZZ)))
        plot(XX1,YY,xlab=xn,ylab=yn, col=2,ylim=c(min(YY,ZZ),max(YY,ZZ)));
        lines(myx, sx%*%YY, lwd=4,col=2, lty=2) ;
        points(XX2,ZZ, col=4, pch=2);
        lines(myx, tx%*%ZZ,col=4, lwd=4, lty=4)
        }
    }

     cat("\n \t \t ======================\n")
     cat("\t\t Curve Test  Procedures \n")
     cat("\t \t ====================== " ); 
     cat("\n The p-value to test H0:", if(is.null(XX2)) "f(x)=0 is \t" else "f1(x)=f2(x) is ", ff, "\n" ); 
     cat( "\n With test statistics equals\t\t",CC, "," ,  "\n Estimated degree of freedom is \t",vv,",\n" );
     if(is.null(XX2) ){  
     cat("\n Estimated sigma^2 is\t \t\t",sigma.square,".\n")
     cat(" ============================\n") 
     } else if(!is.null(XX2)&equal.var) {
     cat("\n Equal variances assumed. \n")
     cat(" Estimated common sigma^2 is \t \t", sigma.square,".\n")
     cat(" ====================\n" )
     } else if(!equal.var)  {
     cat("\n Unequal variances assumed.\n")
     cat(" Estimated sigma^2 are\t \t \t", esigma1^2,"  \n")
     cat(" and \t\t\t\t \t" , esigma2^2,".\n")
     cat(" =========================\n")
     }
if(equal.var|is.null(XX2)) 
temp<-list(Esigma2=sigma.square, Est.df=vv,     k0=k0 ,
     fitted=(if(is.null(XX2)) as.vector(eYY) else  
             list(fit1=as.vector(eYY), fit2=as.vector(eZZ))),
     Residual=(if(is.null(XX2)) as.vector(resid1) else 
              list(resid1=as.vector(resid1),
              resid2=as.vector(resid2)))) else  
if(!equal.var)  
temp=list( Esigma2=c(esigma1^2,esigma2^2), Est.df=c(v1,v2),    k0=k0 ,
     fitted=list(fit1=as.vector(eYY), fit2=as.vector(eZZ)),
     Residual=list(resid1=as.vector(resid1),   resid2=as.vector(resid2)));
invisible(temp); 
}

