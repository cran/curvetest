\name{get.weight.function}
\alias{get.weight.function}
\alias{get.data}
\title{ Test Equality of Curves with Homoscedastic or  Heteroscedastic Errors.}
\description{ For each kernel, return a function corresponding to the name specified. This is a utility function.
}
\usage{
  get.weight.function(type)
}
 \arguments{
  \item{type}{ A character string of the name for the kernel type function. }
  
} 

\value{ A kernel function will be returned with attribute 'name' storing the function name for future possible identification of the kernel function.  
} 

\examples{  get.weight.function("Uniform")}
