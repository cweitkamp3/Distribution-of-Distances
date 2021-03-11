\name{DoD-package}
\alias{DoD-package}
\alias{DoD}
\docType{package}
\title{
  Compare the Distributions of Distances
}
\description{
 This package implements the (bootstrapped) Distribution-of-Distances-Test (DoD-Test) proposed in \insertCite{weitkamp2020gromov;textual}{DoD}. 
}
\details{
Given two point clouds (or distance matrices), the DoD-test compares the corresponding distributions of interpoint distances based on a trimmed, empirical Kantorovich distance. It is an asymptotic level-alpha test for the null hypothesis
 
 \strong{H:}  \emph{"The distributions of (interpoint) distances of both inputs are equal"}
   
against the alternative that they are not. The limit distribution required for testing is approximated by the resampling scheme proposed in Section 3 of \insertCite{weitkamp2020gromov;textual}{DoD}.
}
\author{
Christoph Weitkamp

Maintainer: Christoph Weitkamp <cweitka@mathematik.uni-goettingen.de>
}
\references{
\insertAllCited{}
}
\keyword{ package }

\examples{
  \dontrun{
## n= 100
## 
## square = data.frame("x"=numeric(n),"y"=numeric(n))
## square$x = runif(n,0,1)
## square$y = runif(n,0,1)
## 
## disc = data.frame("x"=numeric(n),"y"=numeric(n))
## radius = 0.5
## r = radius*sqrt(runif(n,0,1))
## theta = runif(n,0,1)*2*pi
## disc$x= r*cos(theta)
## disc$y= r*sin(theta)
## 
## DoD.test(square, disc)  
  }
}
