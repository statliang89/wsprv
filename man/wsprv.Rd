% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WSP.R
\name{wsprv}
\alias{wsprv}
\alias{weight_sp}
\title{A weighted selection probability is developed to locate individual rare variants associated with multiple phenotypes.}
\usage{
weight_sp(
  x,
  y,
  alpha = 1,
  penalty.factor = NULL,
  standardize = TRUE,
  type.multinomial = c("grouped", "ungrouped"),
  rep = 100,
  rate = 0.05,
  gamma = 0.01
)
}
\arguments{
\item{x}{A \eqn{n \times (m+p)} matrix with \eqn{n} samples, \eqn{m} covariates and \eqn{p} rare variants where \eqn{m} can be zero, i.e., there does not exist covariates.}

\item{y}{A \eqn{n \times Q} phenotype matrix with \eqn{n} samples and \eqn{Q} phenotypes where \eqn{Q>1}.}

\item{alpha}{The mixing parameter of elastic-net, \code{alpha=1} is the lasso, and \code{alpha=0} is the ridge. Default value is 1.}

\item{penalty.factor}{Separate penalty factors factors can be applied to each coefficient. Can be \code{0} for some variables, which implies no shrinkage, and that variable is always included in the model.}

\item{standardize}{Genotype standardization. Default is \code{TRUE}.}

\item{type.multinomial}{A group lasso penalty is used on the multinomial coefficients for a variable when 'grouped'. It ensures the multinomial coefficents are all in or out. Default is 'grouped'.}

\item{rep}{The number of bootstrap replications. We recommend to use 100 or more to compute weighted selection probability. Default value is 100.}

\item{rate}{A tuning parameter represents rate of degree of freedom to the number of rare variants. Default value is 0.05.}

\item{gamma}{The upper \code{gamma} quantile of selection frequencies of individual variants each phenotype to compute the threshold. Default value is 0.01.}
}
\value{
\item{res}{A matrix contains the order of weighted selection probabilities from the largest to the smallest and the corresponding weighted selection probabilities.}
   \item{eta}{eta used.}
   \item{bootstrap.rep}{The number of bootstrap replications used.}
   \item{rate}{The tuning parameter \code{rate} used.}
   \item{gamma}{The upper \code{gamma} quantile of selection frequencies of individual rare variants each phenotype used.}
}
\description{
Recently, rare variant association studies with multiple phenotypes have drawn a lot of attentions because association signals can be boosted when rare variants are related with more than one phenotype. Most of existing statistical methods to identify rare variants associated with multiple phenotypes are based on a group test, where a gene or a genetic region is tested one at a time. However, these methods are not designed to locate individual rare variants within a gene or a genetic region. We propose a weighted selection probability to locate individual rare variants within a group after a multiple-phenotype based group test finds significance.
}
\details{
The penalty function of \code{elastic-net} is defined as \deqn{\lambda(\alpha||\beta||_1+\frac{(1-\alpha)}{2}||\beta||_2^2),} where \eqn{\alpha} is a mixing proportion of ridge and the lasso, and \eqn{\beta} is regression coefficients. This penalty is equivalent to the Lasso penalty if \code{alpha=1}. \cr \cr Let \eqn{\eta} be the degree of freedom and it depends on the tuning parameter \eqn{\lambda}, and \code{rate} is computed as \deqn{rate=\frac{\eta}{p},} Note that \eqn{\eta \leq n} is set up in \code{weight_sp} function. \cr \cr Let \eqn{\delta_{\gamma}} be a threshold of \eqn{SF} and it depends on the upper \eqn{\gamma^{th}} qunatile value of \eqn{SF}. Where \eqn{SF=\left\{SF_{11}(\eta),SF_{21}(\eta),\cdots,SF_{pQ}(\eta) \right\}} is a set that contains selection frequencies of individual rare variants each phenotype.
}
\examples{

# Generate simulation data
 n <- 400
 p <- 100
 q <- 5
 MAF <- 0.01
 geno.prob <- rbind((1-MAF)^2,2*(1-MAF)*MAF,MAF^2)
 x <- matrix(NA,n,p)
 set.seed(1)
 for(i in 1:p) x[,i] <- sample(0:2,n,prob=geno.prob,replace=TRUE)
 beta <- c(rep(3.0,10),rep(0,(p-10)))
 cova <- matrix(0.75,q,q)
 diag(cova) <- 1
 require(mnormt)
 err.mat <- rmnorm(n,rep(0,q),cova)

 y1 <- x \%*\% beta+err.mat[,1]
 y2 <- x \%*\% beta+err.mat[,2]
 y <- cbind(y1,y2,err.mat[,3:5])
 # Weighted selection probabilities for individual rare variants without covariates.
 #If rep=100, time consuming.
 wsp.rv1 <- weight_sp(x,y,rep=5) # continuous phenotypes

 # Weighted selection probabilities for individual rare variants with covariates.
 #If rep=100, time consuming.
 cx <- cbind(rnorm(n),sample(0:1,n,replace=TRUE))
 x <- cbind(cx,x)
 penalty.factor <- c(rep(0,2),rep(1,p))
 colnames(x) <- c('Age','Gender',paste0('V',3:102))
 wsp.rv2 <- weight_sp(x,y,penalty.factor=penalty.factor,rep=5) # continuous phenotypes


}
