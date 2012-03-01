\name{getoptimalalpha}
\alias{getoptimalalpha}
\title{
 Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors.
}
\description{
 To calculate the optiomal smoothing alpha if it si not specified. 
}
\usage{
getoptimalalpha(formula, data, plotit = F)
}
 \arguments{
\item{formula}{
 Formula to do the regression.  
 }
  \item{data}{
 A data frame of n rows by 2 columns. The column names should agree with the variable names specified in the formula. 
 }
  \item{plotit}{
  Whether plot will be geneated to show how the different choices of alpha will affect the generalized cross validation values.
 }
}
\details{
  When this routine is invoked, it will fit a sereies regressions specified by the formula on data set. For each one, the generalized cross validation
  will be calculated and the "best" (minimal) GCV will be found with the corresponding alpha returned. 
 }
\value{
A numeric value of alpha value will be returned. }
 
 \references{Zhongfa Zhang, et al: Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors. To appear  }
\author{Zhongfa Zhang, Jiayang Sun}


 
\seealso{
 curvefit, curvetest. 
}
\examples{  
    x1=seq(0,1, length=n1<-50); f1<-function(x){x*(1-x)+sin(2*pi*x)}; 
    y1=f1(x1)+rnorm(n1, 0, 0.2) 
    getoptimalalpha(formula=y~x, data.frame(x=x1, y=y1), plotit = TRUE) 
}
 
