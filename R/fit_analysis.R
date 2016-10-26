########################################################################
# fit_analysis.R
# License: MIT

# ""
# Jesse Tweedle
# Oct, 2016
########################################################################

fit_glmnet <- function(X,c,alpha=1,nlambda=100,lambda.min.ratio=1e-5,penalty=rep_len(1,dim(X)[2])) {
  fit <- glmnet(X,c,alpha=1,nlambda=nlambda,intercept=FALSE,lower.limits=0,upper.limits=1,lambda.min.ratio=lambda.min.ratio)
  
  # pick last iteration.
  t = min(length(fit$df),nlambda)
  
  y <- coef(fit) # could make this more efficient by calling coef(fit,s="correct lambda here")
  pred <- predict(fit,newx=X,s=c(fit$lambda[t])) # some small lambda.
  #df <- fit$df[t] #/ (R*N+N^2)
  
  coefs <- y %>% summary() %>% tbl_df() %>% 
    filter(j==t) %>% # take the last iteration from glmnet
    mutate(i=i-1,j=1) # remove intercept row (i=i-1), and add a column so I can convert to matrix (j=1)
  
  
  return(list(coefs=coefs,pred=pred))
}

# analysis:

# Check A

# Check G

# Check specificity, sensitivity.