########################################################################
# main_x.R
# Clean, load and run relevant files.
# License: MIT

# Project: firm-network-lasso
# Jesse Tweedle
# Oct 18, 2016
########################################################################

source("R/helpers.R")
startup()

## Initialize things.
#set.seed(9) 
R <- 2  # Number of regions
N <- 1000 # Number of firms
Z <- 2  # Number of industries; each industry must have more than one firm, 
        # or glmnet fails (at least until I add more equations).
region_density <- 0.1
firm_density <- 0.1

# Create fake data.
args <- initialize_fake_s(R=R,Z=Z,N=N,region_density=region_density,firm_density=firm_density)

# Save returned arguments.
beta <- args$beta
iz <- args$iz
I <- args$I[,1]
s <- sinit <- args$s
Er <- args$Er
En <- args$En

# Er[Er>0] <- 0
# En[En>0] <- 0

# doesn't make any sense. with more options, it can't do as well?
# must have something to do with penalty scaling. says it scales up to nvars anyway,
# which means penalty prob really bad for non-sparse edges?

# penalty <- rbind(Er,En)
# rm(Er,En)
# dim(penalty) <- c((R+N)*N,1)
# 
# # A test.
# if (FALSE) { 
#   dim(penalty) <- c(R+N,N)
#   erx <- penalty[1:R,1:N]
#   erx + Er # all elements should be = 2
#   enx <- penalty[(R+1):(R+N),1:N]
#   enx + En # all elements should be = 2
# }

# need to switch Er, En around so that they're organized by industry.

# test.
rm(args)

# For input into glmnet_z
preprocess_args <- list(Er=Er,En=En,I=I,iz=iz,s=s,R=R,N=N)

# Apply glmnet_z to each industry, and return y,X_mc and X_ag that have been
# reduced to the possible of nonzero elements. Add a bit that doesn't penalize links I think
# should be there.
cl <- makeCluster(2)
res <- parLapply(cl=cl,X=1:Z,fun=glmnet_z,args=preprocess_args)
stopCluster(cl)
#res <- lapply(1:Z,FUN=glmnet_z,args=preprocess_args)
gc()

somey <- unlist(res)[seq(1,3*Z,3)]
# indexed by z. so 1st element is y_1, second is y_2, etc.
# now apply these to X_a and X_g, somehow.
yyy <- Reduce(rbind,somey) %>% summary() %>% tbl_df()

# Combine the returned X_mc by industry z
X_mc <- unlist(res)[seq(2,3*Z,3)] %>% bdiag()

# Combine the returned X_ag by industry z
X_ag <- Reduce(cbind, unlist(res)[seq(3,3*Z,3)]) #%>% str() #%>% bdiag()

# Create concordance between minimized coefficient vector (size < (R+N)*N) and original coefficient vector (size=(R+N)*N)
nonzero_vars <- (colSums(X_mc)>0) %>% 
  tbl_df() %>% 
  rownames_to_column(var="i_original") %>% 
  filter(value==TRUE) %>% 
  rownames_to_column(var="i_hat") %>% 
  select(-value) %>% 
  mutate(i_original=as.numeric(i_original),i_hat=as.numeric(i_hat))

# ^^ this could be faster ^^

# Reduce the matrices to match the number of parameters in y.
X_mc <- X_mc[,colSums(X_mc)>0]
X_ag <- X_ag[,colSums(X_ag)>0]

# Reorganize vector of firm sales to match industry ordering returned from glmnet_z lapply
cz_list <- 
  lapply(1:Z,FUN=function(z) {
    s[iz[z,]==1]
  })

# And apply them together to get vector of size N.
szs <- Reduce(c,cz_list)
c_mc <- szs # RHS for market clearing equations
c_a <- rep_len(1,R) # RHS for rowSums(A) = 1 equations
c_g <- rep_len(1,N) #1-beta # RHS for rowSums(G) = 1-beta equations.

# Apply it all together.
c <- c(c_mc,c_a,c_g)
X <- rbind(X_mc,X_ag)

rm(c_mc,X_mc,X_ag,res)
# deviance for some reason can go above 100%? I just want it to be as high as possible.
glmnet.control(devmax = 5) 

# fake an element of the last row, for some reason.
# X[dim(X)[1],sample.int(n=dim(X)[2],size=1)]<-1
# xgb.DMatrix(data=X, label = c)
# fit <- xgboost(data=as.matrix(X), label=c, max.depth=2, nthread=2, nround=20, objective = "reg:linear")

# have a lot of lambdas. just want to get as close as possible to the data,
# don't need to worry about overfitting. it should exactly fit the "training" data.
#,lower.limits=0,upper.limits=1) <- may or may not need this.
fit <- glmnet(X,c,alpha=1,nlambda=100,intercept=FALSE,lower.limits=0,upper.limits=1,lambda.min.ratio=0.0001)
y <- coef(fit) # could make this more efficient by calling coef(fit,s="correct lambda here")
print(fit)

pred <- predict(fit,newx=X,s=c(6.420e-05))
s_hat <- pred[1:N,1]
a_eq <- pred[(N+1):(R+N),1]
g_eq <- pred[(R+N+1):(R+2*N),1]

# compare pred to c.
lm(c ~ pred[,1]) %>% summary() # ok...
a_eq %>% summary()
g_eq %>% summary()
df <- tibble(s_hat=s_hat,s=szs)
ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()

xxx

# Convert coefficient x lambda matrix to a useful single vector from one lambda.
yx <- y %>% summary() %>% tbl_df() %>% 
  filter(j==max(j)) %>% # take the last iteration from glmnet
  mutate(i=i-1,j=1) # remove intercept row (i=i-1), and add a column so I can convert to matrix (j=1)

# Merge back on the original parameter indices; e.g., go from the minimized y to the original y.
yx <- yx %>% left_join(nonzero_vars,by=c("i"="i_hat")) %>% select(-i,i=i_original)

# Convert to sparse matrix.
yx <- yx %>% df_to_s(dims=c((R+N)*N,1))

# Calculate sparsity
nnz <- fit$df[length(fit$df)] / (R*N+N^2)
print(str_c("sparsity: ",nnz))
rm(fit,y)

# y---first R are a_{.1}, next N are g_{.1}, etc. so split by R+N.
# this works like magic. must be filled bycolumn (e.g., equivalent of byrow=FALSE)
dim(yx) <- c(R+N,N)

# Get region-firm expenditure matrix A
A <- yx[1:R,1:N]

# Get firm-firm expenditure share matrix G
G <- yx[(R+1):(R+N),1:N]

rm(yx)

# Solve for implied s; if firm i has no customers, s_hati will be zero; this may require some post-processing
# whoops, I forgot g probably includes beta?
#s_hat <- t(A) %*% (I %>% matrix(nrow=R,ncol=1)) + t(G) %*% (((1-beta)*s) %>% matrix(nrow=N,ncol=1)) 
df <- tibble(s_hat=s_hat[,1],s=szs)

# Compare the non-zero firm sizes to the original firm sizes on a log scale.
ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10() 

# How does it do for big firms (> median size).
q <- quantile(df$s,probs=0.5)
# lm(s_hat %>% log() ~ s %>% log(),data=df %>% filter(df$s > q & s_hat>0)) %>% summary()
#ggplot(df %>% filter(df$s > q & s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.5) + scale_x_log10() + scale_y_log10() 

#ggplot(df %>% rownames_to_column() %>% gather(type,value,s_hat:s) %>% filter(value>0))+stat_density(aes(x=value,colour=type),position="dodge",geom="line") + scale_x_log10()

achk <- rowSums(A) %>% summary() # pretty far from correct, should all be 1.
gchk <- (rowSums(G) + beta) %>% summary() # pretty far from correct, should all be 1.

# How sparse are A and G combined? Are the rowSums equations satisfied? s?
lm(s_hat %>% log() ~ s %>% log(),data=df %>% filter(s_hat>0)) %>% summary()
print(str_c("sparsity: ",nnz,"; "))
achk %>% print() # looks to be slightly off, about 60% too high., G too low.
gchk %>% print()
