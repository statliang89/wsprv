#'A weighted selection probability is developed to locate individual rare variants associated with multiple phenotypes.
#'@name wsprv
#'@description Recently, rare variant association studies with multiple phenotypes have drawn a lot of attentions because association signals can be boosted when rare variants are related with more than one phenotype. Most of existing statistical methods to identify rare variants associated with multiple phenotypes are based on a group test, where a gene or a genetic region is tested one at a time. However, these methods are not designed to locate individual rare variants within a gene or a genetic region. We propose a weighted selection probability to locate individual rare variants within a group after a multiple-phenotype based group test finds significance.
#'@param x A \eqn{n \times (m+p)} matrix with \eqn{n} samples, \eqn{m} covariates and \eqn{p} rare variants where \eqn{m} can be zero, i.e., there does not exist covariates.
#'@param y A \eqn{n \times Q} phenotype matrix with \eqn{n} samples and \eqn{Q} phenotypes where \eqn{Q>1}.
#'@param alpha The mixing parameter of elastic-net, \code{alpha=1} is the lasso, and \code{alpha=0} is the ridge. Default value is 1.
#'@param penalty.factor Separate penalty factors factors can be applied to each coefficient. Can be \code{0} for some variables, which implies no shrinkage, and that variable is always included in the model.
#'@param standardize Genotype standardization. Default is \code{TRUE}.
#'@param type.multinomial A group lasso penalty is used on the multinomial coefficients for a variable when 'grouped'. It ensures the multinomial coefficents are all in or out. Default is 'grouped'.
#'@param rep The number of bootstrap replications. We recommend to use 100 or more to compute weighted selection probability. Default value is 100.
#'@param rate A tuning parameter represents rate of degree of freedom to the number of rare variants. Default value is 0.05.
#'@param gamma The upper \code{gamma} quantile of selection frequencies of individual variants each phenotype to compute the threshold. Default value is 0.01.
#'@details The penalty function of \code{elastic-net} is defined as \deqn{\lambda(\alpha||\beta||_1+\frac{(1-\alpha)}{2}||\beta||_2^2),} where \eqn{\alpha} is a mixing proportion of ridge and the lasso, and \eqn{\beta} is regression coefficients. This penalty is equivalent to the Lasso penalty if \code{alpha=1}. \cr \cr Let \eqn{\eta} be the degree of freedom and it depends on the tuning parameter \eqn{\lambda}, and \code{rate} is computed as \deqn{rate=\frac{\eta}{p},} Note that \eqn{\eta \leq n} is set up in \code{weight_sp} function. \cr \cr Let \eqn{\delta_{\gamma}} be a threshold of \eqn{SF} and it depends on the upper \eqn{\gamma^{th}} qunatile value of \eqn{SF}. Where \eqn{SF=\left\{SF_{11}(\eta),SF_{21}(\eta),\cdots,SF_{pQ}(\eta) \right\}} is a set that contains selection frequencies of individual rare variants each phenotype.
#'@importFrom glmnet glmnet
#'@importFrom mnormt rmnorm
#'@importFrom stats median quantile rnorm
#'@returns
#'    \item{res}{A matrix contains the order of weighted selection probabilities from the largest to the smallest and the corresponding weighted selection probabilities.}
#'    \item{eta}{eta used.}
#'    \item{bootstrap.rep}{The number of bootstrap replications used.}
#'    \item{rate}{The tuning parameter \code{rate} used.}
#'    \item{gamma}{The upper \code{gamma} quantile of selection frequencies of individual rare variants each phenotype used.}
#'@examples
#'
#' # Generate simulation data
#'  n <- 400
#'  p <- 100
#'  q <- 5
#'  MAF <- 0.01
#'  geno.prob <- rbind((1-MAF)^2,2*(1-MAF)*MAF,MAF^2)
#'  x <- matrix(NA,n,p)
#'  set.seed(1)
#'  for(i in 1:p) x[,i] <- sample(0:2,n,prob=geno.prob,replace=TRUE)
#'  beta <- c(rep(3.0,10),rep(0,(p-10)))
#'  cova <- matrix(0.75,q,q)
#'  diag(cova) <- 1
#'  require(mnormt)
#'  err.mat <- rmnorm(n,rep(0,q),cova)
#'
#'  y1 <- x %*% beta+err.mat[,1]
#'  y2 <- x %*% beta+err.mat[,2]
#'  y <- cbind(y1,y2,err.mat[,3:5])
#'  # Weighted selection probabilities for individual rare variants without covariates.
#'  #If rep=100, time consuming.
#'  wsp.rv1 <- weight_sp(x,y,rep=5) # continuous phenotypes
#'
#'  # Weighted selection probabilities for individual rare variants with covariates.
#'  #If rep=100, time consuming.
#'  cx <- cbind(rnorm(n),sample(0:1,n,replace=TRUE))
#'  x <- cbind(cx,x)
#'  penalty.factor <- c(rep(0,2),rep(1,p))
#'  colnames(x) <- c('Age','Gender',paste0('V',3:102))
#'  wsp.rv2 <- weight_sp(x,y,penalty.factor=penalty.factor,rep=5) # continuous phenotypes
#'
#'
#'@export
weight_sp <- function(x,y,
                      alpha=1,
                      penalty.factor = NULL,standardize=TRUE,type.multinomial=c('grouped','ungrouped'),
                      rep=100,rate=0.05,gamma=0.01){

  x <- as.matrix(x)
  n <- as.integer(nrow(x))
  p <- as.integer(ncol(x))
  if(floor(p*rate) >= n){
    seq.df <- n
  }else{
    seq.df <- floor(p*rate)
  }
  if (is.null(colnames(x))){
    var.names <- paste('V',1:ncol(x),sep='')
  }else{
    var.names <- colnames(x)
  }
  if(is.matrix(y)==FALSE) stop('y should be a matrix.')
  if (nrow(y)!=n) stop('x and y have different number of rows.')
  if(is.null(penalty.factor)) penalty.factor <- rep(1, ncol(x))
  if(sum(penalty.factor==0)!=0) seq.df <- seq.df + sum(penalty.factor==0)
  type.multinomial <- match.arg(type.multinomial)

  q <- as.integer(ncol(y))
  boot.n <- n
  boot.rep <- rep
  boot.mat <- matrix(0,p,q)
  boot.app <- rep(0,p)
  for(i in 1:boot.rep){
    set.seed(125*i)
    wt.bt <- sample(n,boot.n,replace=TRUE)
    x.boot <- x[wt.bt,]
    boot.sum <- apply(x.boot,2,sum)
    if(min(boot.sum)==0){
      wh <- which(boot.sum==0)
      wh.n <- which(boot.sum!=0)
      x.boot <- x.boot[,-wh]
    }else{
      wh.n <- 1:ncol(x.boot)
    }
    y.boot <- y[wt.bt,]
    boot.app[wh.n] <- boot.app[wh.n] + 1

    families <- rep('gaussian',q)
    families.check <- apply(y,2,function(x) length(table(x)))
    if(sum(families.check==2)!=0) families[which(families.check==2)] <- 'binomial'
    if(sum(families.check==3|families.check==4|families.check==5)!=0) families[which(families.check==3|families.check==4|families.check==5)] <- 'multinomial'

    penalty.factor.temp <- penalty.factor[wh.n]
    for(j in 1:q){
      if(families[j]=='multinomial'){
        opt.lambda <- df_lambda(x=x.boot,y=y.boot[,j],penalty.factor=penalty.factor.temp,alpha=alpha,family=families[j],
                                seq.df=seq.df,standardize=standardize,type.multinomial=type.multinomial)

        fit <- glmnet(x=x.boot,y=y.boot[,j],penalty.factor=penalty.factor.temp,family=families[j],alpha=alpha,
                      lambda=opt.lambda$lambda,standardize=standardize,type.multinomial=type.multinomial)
        boot.mat[,j][wh.n] <- boot.mat[,j][wh.n]+as.numeric(fit$beta[[1]]!=0)

      }else{
        opt.lambda <- df_lambda(x=x.boot,y=y.boot[,j],penalty.factor=penalty.factor.temp,alpha=alpha,family=families[j],
                                seq.df=seq.df,standardize=standardize)

        fit <- glmnet(x=x.boot,y=y.boot[,j],penalty.factor=penalty.factor.temp,family=families[j],alpha=alpha,
                      lambda=opt.lambda$lambda,standardize=standardize)
        boot.mat[,j][wh.n] <- boot.mat[,j][wh.n]+as.numeric(fit$beta!=0)
      }
    }
  }

  if(!is.numeric(penalty.factor)) stop('penalty.factor should be a numeric vector.')
  if(sum(penalty.factor==0)==0){
    counts.qtl <- quantile(boot.mat,(1-gamma))
    boot.sp.mat <- boot.mat/boot.app
    wts.temp <- NA
    for(k in 1:q){
      wts.temp[k] <- sum(boot.mat[,k][boot.mat[,k] > counts.qtl])
    }
    wts <- wts.temp/sum(wts.temp)
    if(sum(is.na(wts))!=0){
      wts[is.na(wts)] <- 0
    }else{
      wts <- wts
    }
    wsp <- apply(boot.sp.mat%*%diag(wts),1,max)
    wsp <- wsp/max(wsp,na.rm=TRUE)
    wsp.order <- order(wsp,decreasing=TRUE)
    mat <- cbind(wsp.order,wsp[wsp.order])
    colnames(mat) <- c('variable','wsp')
    rownames(mat) <- var.names[wsp.order]
  }else{
    cov.loc <- which(penalty.factor==0)
    boot.mat.temp <- boot.mat[-cov.loc,]
    counts.qtl <- quantile(boot.mat.temp,(1-gamma))
    boot.sp.mat.temp <- boot.mat.temp/boot.app[-cov.loc]
    wts.temp <- NA
    for(k in 1:q){
      wts.temp[k] <- sum(boot.mat.temp[,k][boot.mat.temp[,k] > counts.qtl])
    }
    wts <- wts.temp/sum(wts.temp)
    if(sum(is.na(wts))!=0){
      wts[is.na(wts)] <- 0
    }else{
      wts <- wts
    }
    wsp <- apply(boot.sp.mat.temp%*%diag(wts),1,max)
    wsp <- wsp/max(wsp,na.rm=TRUE)
    wsp.order <- order(wsp,decreasing=TRUE)
    wsp.order.mat <- wsp.order + sum(penalty.factor==0)
    mat <- cbind(wsp.order.mat,wsp[wsp.order])
    colnames(mat) <- c('variable','wsp')
    var.names.mat <- var.names[-cov.loc]
    rownames(mat) <- var.names.mat[wsp.order]
  }
  return(list(res=mat,eta=seq.df,bootstrap.rep=boot.rep,rate=rate,gamma=gamma))
}



