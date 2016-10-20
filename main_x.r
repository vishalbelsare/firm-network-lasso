########################################################################
# main.R
# Clean and load relevant files
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
########################################################################

library(dplyr)
library(tibble)
library(Matrix)
library(stringr)
library(tidyr)
library(ggplot2)
library(glmnet)

rm(list=ls(all=TRUE))
gc()

options( warn = 2 )

# Tip from http://jeromyanglim.tumblr.com/post/33418725712/how-to-source-all-r-files-in-a-directory
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir(paste0(getwd(),"/R"))

set.seed(9) 
R <- 20
N <- 2000
Z <- 10 # industries

args <- initialize_fake_s(R,Z,N)

iz <- args$iz
# then split on industry rows? apply lasso? etc.

I <- args$I[,1]
s <- sinit <- args$s

#iz %>% summary() %>% tbl_df() %>% filter(i==1) %>% select(-x,-i)

# z <- 2 # industry 1, to start.
# sz <- s[iz[z,]==1] # select
# cz <- sz
# z_len <- length(sz)
# 
# I_mc <- I %>% t() %>% list() %>%  rep(z_len) %>% bdiag()
# s_mc <- s %>% t() %>% list() %>%  rep(z_len) %>% bdiag()
# X_mc <- cbind(I_mc,s_mc)
# 
# X_mc <- c(I,s) %>% t() %>% list() %>%  rep(z_len) %>% bdiag()
# fit <- glmnet(X_mc,sz,intercept=FALSE,lower.limits=0,upper.limits=1)
# y <- coef(fit)
# # # pick one with max df.
# ncoef <- which(fit$df==max(fit$df))[1]
# yz <- y %>% summary() %>% tbl_df() %>% filter(j==ncoef) %>% select(-j,-x) %>% mutate(i=i-1)
# nnz <- dim(yz)[1]
# yzs <- sparseMatrix(i=yz$i,j=rep_len(1,nnz),x=1,dims=c((R+N)*z_len,1))
#can either return this directly, or process it to save a and g?

# X_mc <- (X_mc %*% (yzs %>% to_sdiag())) %>% summary() %>% tbl_df() %>% filter(x>0)

glmnet_z <- function(z,R,N) {
  gc()
  print(str_c("industry: ", z))
  sz <- s[iz[z,]==1] # select 
  cz <- sz
  
  z_len <- length(sz)
  
  X_mc <- c(I,s) %>% t() %>% list() %>%  rep(z_len) %>% bdiag()
  # I_mc <- I %>% t() %>% list() %>%  rep(z_len) %>% bdiag()
  # s_mc <- s %>% t() %>% list() %>%  rep(z_len) %>% bdiag()
  # X_mc <- cbind(I_mc,s_mc)
  
  fit <- glmnet(X_mc,sz,alpha=0.95,nlambda=100,intercept=FALSE,lower.limits=0,upper.limits=1)
  y <- coef(fit)
  # pick one with max df.
  ncoef <- which(fit$df==max(fit$df))[1]
  yz <- y %>% summary() %>% tbl_df() %>% filter(j==ncoef) %>% select(-j) %>% mutate(i=i-1) 

  nnz <- dim(yz)[1]
  yzs <- sparseMatrix(i=yz$i,j=rep_len(1,nnz),x=1,dims=c((R+N)*z_len,1))
  X_mc <- (X_mc %*% (yzs %>% to_sdiag())) %>% summary() %>% tbl_df() %>% filter(x>0)
  X_mc <- X_mc %>% df_to_s(dims=c(z_len,(R+N)*z_len))

  X_ag <- Reduce(cbind, .sparseDiagonal(R+N) %>% list() %>% rep(z_len))
  # this is A and G together.
#  X_ag %>% str()
  # ok this works
  X_ag <- X_ag %*% (yzs %>% to_sdiag()) %>% summary() %>% tbl_df() %>% filter(x>0) %>% df_to_s(dims=c(R+N,(R+N)*z_len))
  
  return(list(y_z=yzs,X_mc=X_mc,X_ag=X_ag))
}

# number of industries = Z. might have to keep I_mc, etc, outside the function?
res <- lapply(1:Z,FUN=glmnet_z,R=R,N=N)

# Z = length(res)
somey <- unlist(res)[seq(1,3*Z,3)]
# indexed by z. so 1st element is y_1, second is y_2, etc.
# now apply these to X_a and X_g, somehow.
yyy <- Reduce(rbind,somey)

# return y too. then slice the X_mcs from the return and combine them, then use y to create the X_a and X_g?
X_mc <- unlist(res)[seq(2,3*Z,3)] %>% bdiag()
X_ag <- Reduce(cbind, unlist(res)[seq(3,3*Z,3)]) #%>% str() #%>% bdiag()

# so here, either do cbind or bdiag, I think cbind.
# I think c.
# ok, sure, but do this by z.
# X_ag <- Reduce(cbind, .sparseDiagonal(R+N) %>% list() %>% rep(N))
# # this is A and G together.
# 
# # ok this works
# X_ag <- X_ag %*% (yyy %>% to_sdiag()) %>% summary() %>% tbl_df() %>% filter(x>0) %>% df_to_s(dims=c(R+N,(R+N)*N))

cz_list <- 
  lapply(1:Z,FUN=function(z) {
    s[iz[z,]==1]
  })
szs <- Reduce(c,cz_list)
c_mc <- szs
c_a <- rep_len(1,R)
c_g <- 1-args$beta

c <- c(c_mc,c_a,c_g)
X <- rbind(X_mc,X_ag)

###### DO THIS NEXT. APPLY FIRST STEP ZEROS TO THESE MATRICES.
# X_ari <- Reduce(cbind, .sparseDiagonal(R) %>% list() %>% rep(N))
# X_gij <- Reduce(cbind, .sparseDiagonal(N) %>% list() %>% rep(N))
# 
# X <- rbind(X_mc,bdiag(X_ari,X_gij))

# c_ari <- rep_len(1,R)
# c_gij <- 1-args$beta
# 
# c <- c(c_mc,c_ari,c_gij)

# deviance for some reason can go above 100%? I just want it to be as high as possible.
glmnet.control(devmax = 5) 

# have a lot of lambdas. just want to get as close as possible to the data,
# don't need to worry about overfitting. it should exactly fit the "training" data.
#,lower.limits=0,upper.limits=1) <- may or may not need this.
fit <- glmnet(X,c,alpha=1,nlambda=100,intercept=FALSE,lower.limits=0,upper.limits=1)
y <- coef(fit)
# print(fit$df)
# print(fit)

# next, rearrange shit to get back A, keeping in mind the re-organization of s via z.
# last one is to get rid of intercept index.
yx <- y %>% summary() %>% tbl_df() %>% filter(j==max(j)) %>% mutate(i=i-1,j=1)
yx <- yx %>% df_to_s(dims=c((R+N)*N,1))
nnz <- fit$df[length(fit$df)] / (R*N+N^2)
print(str_c("sparsity: ",nnz))

#y---first R are a_{.1}, next N are g_{.1}, etc. so split by R+N.
#yx %>% matrix(nrow=R+N,ncol=N,byrow=FALSE)
dim(yx) <- c(R+N,N)

A <- yx[1:R,1:N]
G <- yx[(R+1):(R+N),1:N]

# doesn't work, because now they split up like a,g,a,g,a,g,etc.
# ok, cool. now it's like (a,g)_1,(a,g)_2, etc. but 1,2,...are organized by industry.

# Rearrange y to get A
# y_sa <- yx %>% filter(i<=R*N) %>% mutate(col=floor((i-1)/R)+1,row=i-(col-1)*R)
# A <- sparseMatrix(i=y_sa$row,j=y_sa$col,x=y_sa$x,dims=c(R,N))
# 
# # Rearrange y to get G
# y_ga <- yx %>% filter(i>R*N) %>% mutate(i=i-R*N) %>% mutate(col=floor((i-1)/N)+1,row=i-(col-1)*N)
# G <- sparseMatrix(i=y_ga$row,j=y_ga$col,x=y_ga$x,dims=c(N,N))

# Do A and G satisfy the equations? #sx <- solve_s(R,N,list(beta=args$beta,ir=args$ir,A=A,G=G))

# Solve for implied s
s_hat <- t(A) %*% (I %>% matrix(nrow=R,ncol=1)) + t(G) %*% (((1-args$beta)*s) %>% matrix(nrow=N,ncol=1)) 
df <- tibble(s_hat=s_hat[,1],s=szs)
lm(s_hat %>% log() ~ s %>% log(),data=df %>% filter(s_hat>0)) %>% summary()
ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.5) + scale_x_log10() + scale_y_log10() 


# How does it do for big firms (> median size).
q <- quantile(df$s,probs=0.5)
lm(s_hat %>% log() ~ s %>% log(),data=df %>% filter(df$s > q & s_hat>0)) %>% summary()
#ggplot(df %>% filter(df$s > q & s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.5) + scale_x_log10() + scale_y_log10() 

#ggplot(df %>% rownames_to_column() %>% gather(type,value,s_hat:s) %>% filter(value>0))+stat_density(aes(x=value,colour=type),position="dodge",geom="line") + scale_x_log10()

# How sparse are A and G combined?
nnz <- fit$df[length(fit$df)] / (R*N+N^2)
achk <- rowSums(A) %>% summary() # pretty far from correct, should all be 1.
gchk <- (rowSums(G) + args$beta) %>% summary() # pretty far from correct, should all be 1.
#sprintf("sparsity: %.4f, A: %.4f, G: %.4f",nnz,achk,gchk)
print(str_c("sparsity: ",nnz,"; "))
achk %>% print()
gchk %>% print()
