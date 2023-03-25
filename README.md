# Distribution-of-Distances

This packages implements the bootstrapped Distribution of Distances (DoD)-test proposed in [Weitkamp et al. (2022)](https://www.tandfonline.com/doi/abs/10.1080/01621459.2022.2127360).
It is the basis for all simulation presented in this paper.

## Installation

In order to install this package, run

```R
library(devtools)
install_github("cweitkamp3/Distribution-of-Distances")
```
## Usage
The usage of this package is quite intuitive and illustrated in the following. First of all, the DoD-statistic of a given sample can be calculated as follows: 
```R
n = 100

square = data.frame("x"=numeric(n),"y"=numeric(n))
square$x = runif(n,0,1)
square$y = runif(n,0,1)

disc = data.frame("x"=numeric(n),"y"=numeric(n))
radius = 0.5
r = radius*sqrt(runif(n,0,1))
theta = runif(n,0,1)*2*pi
disc$x= r*cos(theta)
disc$y= r*sin(theta)

DoD.stat(square, disc)
```
A Bootstrap sample can be calculated like this:
```R
boot.samp = DoD.boot.samp(square, beta = 0.01, p =2)
```
And the bootstrapped DoD-test between two independent samples can be performed like this:
```R
DoD.test(square, disc)
```
Note that the DoD-test is not symmetric in its input arguments due to the Bootstrap procedure used for estimating the quantiles required (see [Weitkamp et al. (2022)](https://www.tandfonline.com/doi/abs/10.1080/01621459.2022.2127360) for a detailed explanation). 
However, the results usually remain comparable.

## References
C. Weitkamp, K. Proksch, C. Tameling and A. Munk. "Distribution of Distances based Object Matching: Asymptotic Inference". Journal of the American Statistical Association, 1-14, 2022.
