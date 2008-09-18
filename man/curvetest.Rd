% --- Source file: man/curvetest.Rd ---
\name{curvetest}
\alias{curvetest}
\title{Test equality of two curves or one curve with a constant.}
\description{Curvetest is used when two smoothing curves with either
   homogeneous or heterogeneous error are equal or not using Tube theory to obtain the tail probability.}
\usage{
    curvetest(XX1, YY, XX2 = NULL, ZZ = NULL, kernel = c("Trio", "Gaussian", "Uniform", "Triweight", "Triangle",
               "Epanechnikov", "Quartic"), equal.var = TRUE, hh = 0.5,plotit = FALSE, conf.level = 0.05, nn = 100)
}
\arguments{
  \item{XX1, YY}{Data that defines the first curve. YY=response, XX1=predictors. Required.}
   \item{XX2, ZZ}{Data that defines the second curve. Optional. If missing,
   the test is to test  curve one to 0.}
   \item{kernel}{The kernel function when smoothing the data.}
  \item{equal.var}{logical value. If TRUE, equal variances are assumed. }
  \item{hh}{Window width for smoothing the data. Can be specified through
 visually  checking the fitting by setting plotit=TRUE.}
  \item{plotit}{Logical value. if TRUE, a scatter plot and smoothing curves
  will be produced. }
  \item{conf.level}{The alpha value.}
  \item{nn}{The number of points equally spaced in the test domain so fitting
  will calculate the values at these points. }
}
\details{
 The package will first try to using smoothing technique to calculate the
 values of the curves. Then using  the Tube formula to calculate the tail
 probability that the two curves are not equal, under null hypothesis that
 they are equal; In one curve case, it calculate the tail probability that the
 curve is not equal to 0, under the null hypothesis that it is 0 in the whole
 test domain.
}
\value{   
  If plotit is true, the scatter plot of data points with smoothing curves will be    
produced, and the test summarization will be printed out. Otherwise, only the latter is printed out. 
}
\references{Test Equality of Curves, Statistic Sinica,2008,  to appear.  }
\author{Jacob Zhang,  Jiayang Sun}
\note{ First try a window width hh equal to one fourth of the whole length of test interval, specify the plotit=TRUE to check the fitting of curves to the
   points. Then choose a suitable hh by visually checking the plot. This may need a few times. 
}
\seealso{locfit, t.test, test}
\examples{
   n1=50; n2=55
   x1=seq(0,1, length=n1); x2=seq(0, 1, length=n2);
   y1=x1*(1-x1)+rnorm(n1, 0, 0.02)
   y2=x2*(1-x2)+rnorm(n2, 0, 0.01)
   curvetest(x1, y1, x2, y2, equal.var=FALSE, plotit=TRUE) 
   #choose hh by checking the fitting of curves.
   curvetest(x1, y1, x2, y2, equal.var=FALSE, hh=0.2,  plotit=TRUE)
}
\keyword{methods} 
