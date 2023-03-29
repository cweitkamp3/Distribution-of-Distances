#' The (bootstrapped) Distribution-of-Distances (DoD) Test
#'
#' Perform the (bootstrapped) DoD-test proposed in \insertCite{DoD2021;textual}{DoD} for two samples of size n and m.
#' @param input.asym either a n x d numeric coordinate matrix (resp. data frame) or a n x n distance matrix. The asymptotic distribution of the test statistic is approximated based on this input.
#' @param input.comp either a m x d numeric coordinate matrix (resp. data frame) or a m x m distance matrix.
#' @param beta the trimming parameter of the DoD-statistic. It can attain values in [0,1/2).
#' @param p  the order of the DoD-statistic. By definition p>= 1.
#' @param boot.m.out.of.n the bootstrap sample size.
#' @param boot.rep the number of bootstrap repetitions. 
#' @param force.coordinates boolean that indicates whether the input is forced to be interpreted as coordinate matrix.
#' @param method the distance measure to be used in case the input is considered as coordinate matrices. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @details 
#' The DoD-Test is an asymptotic level-alpha test for the null hypothesis
#' 
#' \strong{H:}  \emph{"The distributions of (interpoint) distances of both inputs are equal"}
#'   
#' against the alternative that they are not. The test is based on the (trimmed) empirical Kantorovich distance of order p between the distribution of distances of the input
#' samples \code{input.asym} and \code{input.comp} (for more details see \insertCite{DoD2021;textual}{DoD}).
#'  
#' If both inputs are quadratic, numeric matrices (and \code{force.coordinates}==F), they will be interpreted as distance matrices. Otherwise, they are considered as coordinate matrices
#' and the interpoint distances are calculated with \code{dist} and the method specified.
#' 
#' The limit distribution under the null hypothesis, which is required for testing,
#' is approximated by a resampling scheme (see Section 3 of \insertCite{DoD2021;textual}{DoD}) based on the input \code{input.asym}.    
#' @return 
#' The function returns a list containing the following components: \tabular{ll}{
#' \code{scaled.statistic}:  the value of the DoD-statistic (defined in \insertCite{DoD2021;textual}{DoD}) scaled by mn/(m+n). \cr
#' \tab \cr
#' \code{p.value}: the p-value of the comparison.
#'}
#' @examples
#' n <- 100
#' 
#' square <- data.frame("x" = numeric(n),"y" = numeric(n))
#' square$x <- runif(n,0,1)
#' square$y <- runif(n,0,1)
#' 
#' disc <- data.frame("x" = numeric(n),"y" = numeric(n))
#' radius <- 0.5
#' r <- radius*sqrt(runif(n, 0, 1))
#' theta <- runif(n, 0, 1)*2*pi
#' disc$x <- r*cos(theta)
#' disc$y <- r*sin(theta)
#' 
#' DoD.test(square, disc)
#'@references{
#'   \insertAllCited{}
#' }
#' @importFrom Rdpack reprompt
#' @export
DoD.test <- function(input.asym, input.comp, beta =0.01 ,p=2,boot.m.out.of.n =dim(input.asym)[1], boot.rep = 1000,force.coordinates=F, method = "euclidean"){
  dim.quant <- dim(input.asym)
  dim.comp <- dim(input.comp)

  if(beta <0 || beta >= 0.5)
    stop("The parameter beta has to be in [0,1/2).")
  
  if(p < 1)
    stop("The parameter p has to be greater than or equal to 1.")
  
  dist.mat.input = F
  if(dim.quant[1] == dim.quant[2] && dim.comp[1] == dim.comp[2] && force.coordinates == F){
    dist1 <- input.asym
    dist2 <- input.comp
    dist.mat.input <- T
  } else{
    dist1 <- dist(input.asym, method = method)
    dist2 <- dist(input.comp, method = method)
  }

  if(dist.mat.input){
    dist.vec1 <- sort(dist1[lower.tri(dist1)])
    dist.vec2 <- sort(dist2[lower.tri(dist2)])
  } else {
    dist.vec1 <- sort(as.vector(dist1))
    dist.vec2 <- sort(as.vector(dist2))
    dist1 <- as.matrix(dist1)
    dist2 <- as.matrix(dist2)
  }

  boot.samp <- Bootstrap(boot.rep,  beta, dist1, boot.m.out.of.n, dist.vec1, p)
  stat <- (dim.quant[1]*dim.comp[1])/(dim.quant[1]+dim.comp[1])*trimmed_quantile_diff(dist.vec1, dist.vec2, beta, p)
  pval <- 1 - ecdf(boot.samp)(stat)

  return(list("scaled.statistic" = stat, "p.value" = pval))
}


#' The DoD-Statistic
#' 
#' Calculate the (unscaled) DoD-statistic proposed in \insertCite{DoD2021;textual}{DoD} for two samples of size n and m.
#' @param input.1 either a n x d numeric coordinate matrix (resp. data frame) or a n x n distance matrix. 
#' @param input.2 either a m x d numeric coordinate matrix (resp. data frame) or a m x m distance matrix.
#' @param beta the trimming parameter of the DoD-statistic. It can attain values in [0,1/2).
#' @param p  the order of the DoD-statistic. By definition p>= 1.
#' @param force.coordinates boolean that indicates whether the input is forced to be interpreted as coordinate matrix.
#' @param method the distance measure to be used in case the input is considered as coordinate matrices. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @details 
#' The DoD-Statistic is defined as the (trimmed) empirical Kantorovich distance of order p between the distribution of distances of the input
#' samples \code{input1} and \code{input2} (for more details see \insertCite{DoD2021;textual}{DoD}).
#'  
#' If both inputs are quadratic numeric matrices (and \code{force.coordinates}==F), they will be interpreted as distance matrices. Otherwise, they are considered as coordinate matrices
#' and the interpoint distances are calculated with \code{dist} and the method specified.
#'
#' @return 
#' The function returns the value of the DoD-statistic (defined in \insertCite{DoD2021;textual}{DoD}).
#' 
#' @examples
#' n <- 100
#' 
#' square <- data.frame("x" = numeric(n),"y" = numeric(n))
#' square$x <- runif(n, 0, 1)
#' square$y <- runif(n, 0, 1)
#' 
#' disc <- data.frame("x" = numeric(n),"y" = numeric(n))
#' radius <- 0.5
#' r <- radius*sqrt(runif(n,0,1))
#' theta <- runif(n,0,1)*2*pi
#' disc$x <- r*cos(theta)
#' disc$y <- r*sin(theta)
#' 
#' DoD.stat(square, disc)
#'@references{
#'   \insertAllCited{}
#' }
#' @importFrom Rdpack reprompt
#' @export

DoD.stat <- function(input.1, input.2, beta = 0.01 ,p = 2, force.coordinates = F, method = "euclidean"){
  dim.1 = dim(input.1)
  dim.2 = dim(input.2)

  if(beta < 0 || beta >= 0.5)
    stop("The parameter beta has to be in [0,1/2).")

  if(p < 1)
    stop("The parameter p has to be greater than or equal to 1.")


  dist.mat.input = F
  if(dim.1[1] == dim.1[2] && dim.2[1] == dim.2[2] && force.coordinates == F){
    dist1 <- input.1
    dist2 <- input.2
    dist.mat.input <- T
  } else{
    dist1 <- dist(input.1, method = method)
    dist2 <- dist(input.2, method =method)
  }

  if(dist.mat.input){
    dist.vec1 <- sort(dist1[lower.tri(dist1)])
    dist.vec2 <- sort(dist2[lower.tri(dist2)])
  } else {
    dist.vec1 <- sort(as.vector(dist1))
    dist.vec2 <- sort(as.vector(dist2))
  }

  dod.stat <- trimmed_quantile_diff(sort(dist.vec1), sort(dist.vec2), beta, p)

  return(dod.stat)
}


#' Approximation of the limit distribution under the null hypothesis
#' 
#' Calculate a bootstrap approximation of the limit distribution of the DoD-statistic under the null
#' hypothesis based on a sample of size n (see Section 3 of \insertCite{DoD2021;textual}{DoD} for additional information).
#' @param input either a n x d numeric coordinate matrix (resp. data frame) or a n x n distance matrix. 
#' @param beta the trimming parameter of the DoD-statistic. It can attain values in [0,1/2).
#' @param p  the order of the DoD-statistic. By definition p>= 1.
#' @param boot.m.out.of.n the bootstrap sample size.
#' @param boot.rep the number of bootstrap repetitions. 
#' @param force.coordinates boolean that indicates whether the input is forced to be interpreted as coordinate matrix.
#' @param method the distance measure to be used in case the input is considered as coordinate matrices. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @details 
#' In order to apply the DoD-Test the limit distribution of the DoD-Statistic under the null
#' hypothesis
#'  
#' \strong{H:}  \emph{"The distributions of (interpoint) distances of both inputs are equal."}
#'  
#' is required. This function implements the resampling scheme proposed in Section 3 of \insertCite{DoD2021;textual}{DoD} to approximate this limit distribution.  
#'  
#' If the input is a quadratic numeric matrix (and \code{force.coordinates}==F), it will be interpreted as distance matrix. Otherwise, it is considered as coordinate matrix
#' and the interpoint distances are calculated with \code{dist} and the method specified.
#'
#' @return 
#' The function returns a numeric vector that contains the specified number (given by \code{boot.rep}) of m-out-of-n bootstrap realizations (m=\code{boot.m.out.of.n }) of the resampling
#' scheme proposed in Section 3 of \insertCite{DoD2021;textual}{DoD} in order to approximate the limit distribution of the DoD-Statistic under the null hypothesis.
#' 
#' @examples
#' n <- 100
#' 
#' square <- data.frame("x" = numeric(n),"y" = numeric(n))
#' square$x <- runif(n, 0, 1)
#' square$y <- runif(n, 0, 1)
#' 
#' 
#' DoD.boot.samp(square)
#' 
#'@references{
#'   \insertAllCited{}
#' }
#' @importFrom Rdpack reprompt
#' @export

DoD.boot.samp <- function(input, beta = 0.01 ,p = 2, boot.m.out.of.n = dim(input)[1], boot.rep = 1000, force.coordinates=F, method = "euclidean"){
  dim.input = dim(input)

  if(beta < 0 || beta >= 0.5)
    stop("The parameter beta has to be in [0,1/2).")
  
  if(p < 1)
    stop("The parameter p has to be greater than or equal to 1.")
  
  dist.mat.input <- F
  if(dim.input[1] == dim.input[2] && force.coordinates == F){
    dist1 <- input
    dist.mat.input <- T
  } else{
    dist1 <- dist(input, method = method)
  }

  if(dist.mat.input){
    dist.vec1 <- sort(dist1[lower.tri(dist1)])
  } else {
    dist.vec1 <- sort(as.vector(dist1))
    dist1 <- as.matrix(dist1)
  }

  boot.samp <- Bootstrap(boot.rep,  beta, dist1, boot.m.out.of.n, dist.vec1, p)

  return(boot.samp)
}
