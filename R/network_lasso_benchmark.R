# Function that calls the benchmarking thing for another program.
network_lasso_benchmark <- function(R,N,I,ik,kk,s,beta,summary=FALSE) {
  # source("R/helpers.R") # if I source the directory where this file is, I should also get the helpers. Just need
  # to make sure they don't overwrite over things with the same name. Shouldn't.
  # startup()
  load_libraries()
  
  # these should all be arguments to the function. need to load data elsewhere.
  # I <- args$I[,1]
  # ik <- args$ik
  # kk <- args$kk
  # s <- args$s[,1]
  # beta <- args$beta
  
  upper_bound <- upper_bound(args$A,args$G,region_density,firm_density) 
  
  # Lower bound is a sub-sample of true data A and G.
  # lower_bound <- lower_bound(args$A,args$G,region_sample_frac=0.5,firm_sample_frac=0.1)
  
  x <- upper_bound
  dim(x) <- c((R+N)*N,1)
  nonzero_vars <- x %>% summary() %>% 
    tbl_df() %>% 
    rownames_to_column(var="i_hat") %>% 
    rename(i_original=i) %>%
    select(-j,-x) %>% 
    mutate(i_original=as.numeric(i_original),i_hat=as.numeric(i_hat))
  
  # dim(lower_bound) <- c((R+N)*N,1)
  # x <- cbind(x,lower_bound)
  # rm(lower_bound)
  # penalty <- 1 - 1*x[rowSums(x)>0,2] # next, gonna try to change this to re-weight up firm expenditure relative to region expenditure (e.g., R/N or something)
  # rm(x)
  
  X_mc <- create_X_mc(I,beta,s,upper_bound)
  X_ag <- create_X_ag(I,beta,s,upper_bound)
  X_ind <- create_X_ind(beta=beta,s=s,upper_bound=upper_bound,ik=ik,nonzero_vars=nonzero_vars)
  
  c_mc <- s # RHS for market clearing equations
  c_a <- I #rep_len(1,R) # RHS for rowSums(A) = 1 equations
  c_g <- (1-beta) * s #(1-args$beta)*s #rep_len(1,N) #1-beta # RHS for rowSums(G) = 1-beta equations.
  c_ind <- kk[,1]
  
  # Apply it all together.
  c <- c(c_mc,c_a,c_g,c_ind)
  X <- rbind(X_mc,X_ag,X_ind)
  rm(X_mc,X_ag,X_ind)
  
  # deviance for some reason can go above 100%? I just want it to be as high as possible.
  glmnet.control(devmax = 5) 
  
  nlambda <- 100
  
  # fit; returns coefs, and prediction.
  
  # Try penalty with different for R and N. should be able to use upper_bound for that.
  # penalty <- c(rep_len(1/R,R),rep_len(1/N,N)) %>% to_sdiag() %*% upper_bound
  penalty <- c(I^(0.5),s^(0.025)) %>% to_sdiag() %*% upper_bound
  dim(penalty) <- c(R*N+N^2,1)
  penalty <- penalty[rowSums(penalty)>0,1]
  
  # system.time(fit <- fit_glmnet(X,c,alpha=1,nlambda=nlambda,lambda.min.ratio=1e-12))
  fit <- fit_glmnet(X,c,alpha=1,nlambda=nlambda,lambda.min.ratio=1e-12,penalty=penalty)
  pm <- fit %>% predicted_matrices(R=R,N=N,nonzero_vars=nonzero_vars)
  
  if (summary) {
    prhs <- fit %>% predicted_rhs(R=R,N=N)
    # summary of fit.
    fit_summary(prhs,c_mc,c_a,c_g,c_ind)
  }
  
  # benchmarked things.
  # maybe also s that is internally consistent? nah.
  return(pm) #includes A_hat, G_hat
}