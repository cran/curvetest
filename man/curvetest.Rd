\name{curvetest}
\alias{curvetest.raw}
\alias{curvetest}
\alias{print.curvetest} 
\title{
  Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors.
}
\description{
 Main test routine for formally testing the equality of two curves
  represented by discrete data points over a domain with homoscedacity  or heteroscedacity errors. 
 curvetest is the wrapper of curvetest.raw.
}
\usage{
curvetest(formula, data1 = NULL, data2 = NULL, equal.var = TRUE, alpha = 0.5, bw = NULL, plotit = TRUE, conf.level = 0.05, kernel = c("Triangle", "Gaussian", "Trio", "Uniform", "Triweight", "Epanechnikov", "Quartic"), nn = 100, myx = NULL, bcorrect = "simple",...)
curvetest.raw(fits1, fits2, equal.var, conf.level, plotit)
\method{print}{curvetest}(x,...) 
}
\arguments{
  \item{formula}{ A formula to the data set such as y~x. }
  \item{data1}{A data frame of 2 columns representing the underlying curve 1. The column names must agree with the names in formula.}
  \item{data2}{A data frame for curve 2. If it is NULL, the test is whether curve 1 is statistically equal to 0 over the defining domain.}
  \item{equal.var}{Whether the variances are equal. Default to TRUE. }
  \item{alpha}{Smoothing parameter. }
  \item{bw}{Window bandwidth for both curves.}
  \item{plotit}{WHether plot the fitted curves or not. Default: FALSE.}
  \item{conf.level}{The confidence level to claim the curves are different. Default: 0.05.}
  \item{kernel}{One of the kernel functions to use to fit the curves. Must be one  of "Triangle",  "Gaussian", "Trio","Uniform", "Triweight",  "Epanechnikov", 
    "Quartic". partial match is allowed.}
  \item{nn}{Number of data points in the test domain to calculate the curve values. }
  \item{myx}{x-values  in the test domain to calculate the curve values. If it is specified, nn will be suppressed.}
  \item{bcorrect}{Boundary correction method. Right now, except for 'none', meaning no corrections, the only other opton is 'simple'.}
  \item{\dots}{When plotit is true, plot parameters can be specified such as pch, lty, col etc. }
  \item{fits1,fits2}{ The fitted results from fitting the first or second  curve by curvefit procedure. } 
  \item{x}{Test results from curvetest.} 
 } 
\details{
  The algorithm works by first fitting the curves using local regression method specified by formula  on 
  data1 and/or data 2 with smoothing parameters specified in the function calls. Then it will test on the
   fitted curve 1 and curve 2 to see if they are statistically equal or not.  
}

\value{An R object of class curvetest will be generatd, containing curve fitting and testing results.  
   \item{Statistic,p}{Test statistic and p value of testing whether $f_1(x)==f_2(x)$ or $f_1(x)==0$.}
   \item{eDF}{Estimated degree of freedom of the fitting. }
   \item{equal.var}{The model specification of whether the two variances are equal or not }
   \item{esigma1,esigma2}{Estimated varaince of the fitted curves. }
   \item{k0}{The calculated value from the tube. See detail in paper.} 
   \item{fits1,fits2}{Objects of class 'curvefit' from curvefit ruotine. }
 }
 
\references{Zhongfa Zhang, et al: Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors. To appear}
\author{Zhongfa Zhang, Jiayang Sun}
\seealso{ curvefit, print.curvetest, plot.curvetest}
\examples{
    n1=150; n2=155 ##numbers of data points for the two curves. 
    f1<-f2<-function(x){x*(1-x)+sin(2*pi*x)}; ##True functions.
    x1=seq(0,1, length=n1);
    x2=seq(0, 1, length=n2);
    y1=f1(x1)+rnorm(n1, 0, 0.2)
    y2=f2(x2)+rnorm(n2, 0, 0.2) ###Measured data for the two curves with noises. 
    curvetest(y~x,data.frame(x=x1,y=y1), data.frame(x=x2,y=y2), alpha =  0.7, equal.var=TRUE,plotit=TRUE) 
 }
 
\keyword{curvetest}
\keyword{curvefit}
