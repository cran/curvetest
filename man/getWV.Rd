\name{getWV}
\alias{getWV} 
\title{
  Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors.
}
\description{
  This is a utility function. 
}
\usage{
getWV(x, myx, kernel = kernel, alpha = NULL, bw = NULL, bc = "simple", getit =TRUE)
}
\arguments{
  \item{x}{The x values in the data set that defines the curve. }
  \item{myx}{The x-values in the x-domain that will be used to calculate the curve values.  }
  \item{kernel}{  Kernel functions.  }
  \item{alpha}{  Smoothing parameter.  }
  \item{bw}{Bandwidth. If both alpha and bw are specified, hh will be used instead.  }
  \item{bc}{ Boundary correction method. }
  \item{getit}{Logical, tell the algorithm whether to calculate or to laod the file to speed up calculation. }
}
\details{
  This function will accomplish the bulk load of calculations in the algorithm. If getit=TRUE, the algorithm will run the calculation and save the result in somewhere for future load. Otherwise, it will just load the saved output. 
  This will save a lot of time when do a simulation for a large number of iterations. 
}
\value{
  The function will return the calculated k0, delta,delta2,
   degree of freedon vv, etc that will be used in the curve fit and test step. 
} 
 
   
  
