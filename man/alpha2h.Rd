\name{alpha2h}
\alias{alpha2h}
\title{
   Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors.
}
\description{
For each value $at$ in the defining interval, find a bandwidth  $h$ so that alpha*100 percent of data points specified in $xseq$
should be within the $(x-h, x+h)$ window.  This is a utility function.
}
\usage{
alpha2h(alpha, at, xseq)
}
\arguments{
  \item{alpha}{Smoothing parameter that  for each point in the domain, use a  windown width that should have 
  alpha*100 percent of data points falling in the window.}
  \item{at}{ a point in the x domain.  }
  \item{xseq}{ Sequence of the data points.}
}
 
\value{A numeric value h that will be used as bandwidth in the next curve fitting process.  }
\references{Zhongfa Zhang, et al: Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors. To appear.}
\author{Zhongfa Zhang, Jiayang Sun}
 \seealso{curvetest, curvefit, print.curvetest, plot.curvetest}
\examples{
    x= runif(100)
    (h=alpha2h(0.3, at=0.5, xseq=x))  ##get the window width h around x=.5 so that 30% data points of xseq fall in the area.  
  }
 \keyword{ alpha2h}
 
