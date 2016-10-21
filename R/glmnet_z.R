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
  # Try to manage memory from previous calls. 
  gc()
  print(str_c("industry: ", z))
  
  # Save important arguments.
  I <- args$I
  iz <- args$iz
  s <- args$s
  R <- args$R
  N <- args$N
  
  sz <- s[iz[z,]==1] # select firm sizes from industry z
  cz <- sz # The RHS argument that goes into glmnet.
  
  # Number of firms in this industry.
  z_len <- length(sz)
  
  # The market clearing equations
  X_mc <- c(I,s) %>% t() %>% list() %>%  rep(z_len) %>% bdiag()
  
  # Estimate elastic net / lasso.
  # Objective: find potential non-zero coefficients.
  fit <- glmnet(X_mc,sz,alpha=1,nlambda=100,intercept=FALSE,lower.limits=0,upper.limits=1)
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
