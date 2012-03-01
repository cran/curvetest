\name{curvefit}
\alias{curvefit}
\alias{lines.curvefit}
\alias{lines}
\alias{plot.curvefit}
\alias{plot}
\alias{print.curvefit}
\alias{print}
\alias{distance}
\title{ Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors.
}
\description{
  Fit the smoothing curves. 
}
\usage{
curvefit(formula, data, kernel = "Quartic", alpha = 0.5, bw = NULL, myx, bcorrect = "simple", getit = F)
\method{print}{curvefit}(x,...) 
\method{plot}{curvefit}(x,y=NULL, add = F, get.data = TRUE, ...)
\method{lines}{curvefit}(x,...)
} 
\arguments{
	\item{formula}{ A formula to the data set such as y~x. }
	\item{data}{A data frame of 2 columns representing the underlying curve. 
	The column names must agree with the names in formula.}
	\item{alpha}{Smoothing parameter. Default=0.5.}
	\item{bw}{Window bandwidth for fitting the  curve.}
	\item{kernel}{One of the kernel functions to use to fit the curves. Must be one  of "Triangle",  "Gaussian", "Trio","Uniform", "Triweight",  "Epanechnikov", 
	"Quartic". partial match is allowed.}
	\item{myx}{x-values  in the test domain to calculate the curve values.}
	\item{bcorrect}{Boundary correction method. Right now, except for 'none', meaning no corrections, the only other opton is 'simple'.}
	\item{getit}{unused for this function.}
	\item{add}{logical, Tf true, add the curves to the plot.Otherwise, add fitted lines to the plot. }
	\item{get.data}{logical, not used in this function. }
	\item{x}{The fitted results from fitting the first or second  curve by curvefit procedure. } 
	\item{y}{dummy variable for compatible with parameters in the base definition of plot.}
	\item{\dots}{parameters for plot such as  pch, lty, col etc. }
}
\details{
  For a 2 column data, the curve will be fitted according to formula using local regression method. Boundary corrections can be made. 
  The fitted result will be returned as a 'curvefit' object, that can be plotted and printted by the associated S3 method print and plot. 
 }
\value{ 
An R object of class 'curvefit' will be generated including the fitted values of the curves with original specification of parameters.
}
  \references{Zhongfa Zhang, et al: Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors. To appear  }
\author{Zhongfa Zhang, Jiayang Sun}
 \seealso{ curvefit, print.curvetest, plot.curvetest}
\examples{
     x=seq(0,1, length=n<-150);
     f<-function(x){x*(1-x)+sin(2*pi*x)};
     y=f(x)+rnorm(n, 0, 0.5)
     fit<-curvefit(y~x,data.frame(x=x,y=y), bw=0.4,getit=T) 
     plot(fit)
     lines(fit)
     fit ##print 
 } 
\keyword{curvefit }
\keyword{plot.curvefit}

 
