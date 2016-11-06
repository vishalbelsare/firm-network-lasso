########################################################################
# glmnet_z.R
# Function that calculates lasso on industry subset and returns 
# potential non-zero coefficients.
# License: MIT

# ""
# Jesse Tweedle
# Oct, 2016
########################################################################


glmnet_z <- function(z,args) {
# Apply glmnet, return coefficients from specific lambda. 
  
  # This for parallel on windows; it operates in new environment, doesn't have any of package information.
  # May have to change this if new environment doesn't know the package working directory
  dir <- getwd()
  
  # Get helper functions (including load_libraries())
  # This for parallel on windows; it operates in new environment, doesn't have any of package information.
  source(paste0(dir,"/R/helpers.R"))
  load_libraries(inside_parallel=TRUE)
  
  # Try to manage memory from previous calls. 
  gc()
  sprintf("industry: %d", z)
  
  # Save important arguments.
  I <- args$I
  iz <- args$iz
  Er <- args$Er
  En <- args$En
  s <- args$s
  R <- args$R
  N <- args$N
  rmgc(args)
  
  sz <- s[iz[z,]==1] # select firm sizes from industry z
  cz <- sz # The RHS argument that goes into glmnet.
  
  Erz <- Er[,iz[z,]==1]
  Enz <- En[,iz[z,]==1]
  # Number of firms in this industry.
  z_len <- length(sz)
  
  penalty <- rbind(Erz,Enz)
  rm(Er,En)
  dim(penalty) <- c((R+N)*z_len,1)
  penalty <- 1 - penalty

  # random penalty.
  # penalty <- as.integer(rnorm(n=(R+N)*z_len) > 0)
  
  # if (FALSE) { # a test.
  #   dim(penalty) <- c(R+N,z_len)
  #   erx <- penalty[1:R,1:z_len]
  #   erx + Erz # all elements should be = 2
  #   enx <- penalty[(R+1):(R+N),1:z_len]
  #   enx + Enz # all elements should be = 2
  # }  
  
  # The market clearing equations; would like to remove s_i from each row i. Could just calculate the relevant
  # position in X_mc and set them all to 0? Should only be N operations.
  X_mc <- c(I,s) %>% t() %>% list() %>%  rep(z_len) %>% bdiag()
  
  # Estimate elastic net / lasso.
  # Objective: find potential non-zero coefficients.
  # Include 0 penalties for links identified in a previous step.
  fit <- glmnet(X_mc,sz,alpha=0.1,nlambda=100,intercept=FALSE,lower.limits=0,upper.limits=1,penalty.factor=penalty,lambda.min.ratio=0.0001)
  y <- coef(fit)

  # pick one with max number of non-zerosz.
  ncoef <- which(fit$df==max(fit$df))[1]
  yz <- y %>% summary() %>% tbl_df() %>% filter(j==ncoef) %>% select(-j) %>% mutate(i=i-1) 
  
  nnz <- dim(yz)[1]
  yzs <- sparseMatrix(i=yz$i,j=rep_len(1,nnz),x=1,dims=c((R+N)*z_len,1))
  X_mc <- (X_mc %*% (yzs %>% to_sdiag())) %>% summary() %>% tbl_df() %>% filter(x>0)
  X_mc <- X_mc %>% df_to_s(dims=c(z_len,(R+N)*z_len))
  
  X_ag <- Reduce(cbind, .sparseDiagonal(R+N) %>% list() %>% rep(z_len))
  X_ag <- X_ag %*% (yzs %>% to_sdiag()) %>% summary() %>% tbl_df() %>% filter(x>0) %>% df_to_s(dims=c(R+N,(R+N)*z_len))
  
  return(list(y_z=yzs,X_mc=X_mc,X_ag=X_ag))
}
