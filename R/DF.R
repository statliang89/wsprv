#'@importFrom stats median quantile
#'@importFrom glmnet glmnet
df_lambda <- function(x, y,
                      alpha,
                      family,type.multinomial,seq.df,
                      penalty.factor=NULL,standardize){

  x <- as.matrix(x)
  n=as.integer(nrow(x))
  p=as.integer(ncol(x))
  if(length(y)!=n) stop('x and y have different number of rows.')

  get_range_glmnet <- function(glmnet_fit,DFDF){
    df.loss <- abs(glmnet_fit$df-DFDF)
    out.df <- glmnet_fit$df[which.min(df.loss)]
    out.lambda <- glmnet_fit$lambda[which.min(df.loss)]
    if(sum(df.loss==0)>0){
     lambda_upper <- glmnet_fit$lambda[median(which(df.loss==0))]
     lambda_lower <- glmnet_fit$lambda[median(which(df.loss==0))]
    }else{
     lambda_minmax <- glmnet_fit$lambda[order(df.loss,decreasing=FALSE)[1:2]]
     lambda_upper <- max(lambda_minmax)
     lambda_lower <- min(lambda_minmax)
    }
    list(df=out.df,lambda=out.lambda,minmax=c(lower=lambda_lower,upper=lambda_upper))
  }

  fit.lambda_range <- glmnet::glmnet(x=x,y=y,alpha=alpha,nlambda=nrow(x)*15,
                                     family=family,dfmax=nrow(x)*3,
                                     type.multinomial=type.multinomial,
                                     penalty.factor=penalty.factor,standardize=standardize)

  lambda.sequence <- fit.lambda_range$lambda
  df.sequence <- fit.lambda_range$df
  out.lambda.seq <- NULL

  for(DFDF in seq.df){
    fit.lambda_minmax <- get_range_glmnet(fit.lambda_range,DFDF)
    avg.lambda.minmax <- mean(fit.lambda_minmax$minmax)
    out.lambda.seq <- c(out.lambda.seq,avg.lambda.minmax)
  }

  tmp_glmnet_fit <- glmnet::glmnet(x=x,y=y,
                           alpha=alpha,lambda=out.lambda.seq,family=family,
                           type.multinomial=type.multinomial,
                           penalty.factor=penalty.factor,
                           standardize=standardize)

  list(DF=seq.df,df=tmp_glmnet_fit$df,lambda=out.lambda.seq)
}


