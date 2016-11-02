########################################################################
# main.R
# Clean and load relevant files
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
########################################################################
rm(list=ls(all=TRUE))
gc()

source("R/helpers.R")
startup()

## Initialize things.
set.seed(9)

# data: create A, G, industry IO table
# get upper / lower bounds.

# need to make sure potential non-zero parameters get reshaped in the right way.

# create X_mc, X_ag, X_ind matrices
# and c_mc, c_ag, c_ind

# fit model, save parameters

# evaluate.


# Create fake data.
fakes <- TRUE
if (fakes) {
  print("Initialize fake data.")
  R <- 20  # Number of regions
  N <- 1000 # Number of firms
  K <- 5  # Number of industries; each industry must have more than one firm, 
  # or glmnet fails (at least until I add more equations).
  region_density <- 0.05
  firm_density <- 0.02
  scale <- N

#  args <- initialize_fake_s(R=R,K=K,N=N,region_density=region_density,firm_density=firm_density)
  args <- initialize_fake_economy(R=R,K=K,N=N,region_density=region_density,firm_density=firm_density,scale=scale)
} else {
  print("Load real data.")
  args <- load_data() 
}

# Save returned arguments.
beta <- args$beta
ik <- args$ik
I <- args$I[,1]
s <- args$s[,1]

## Fake A and G
A <- args$A
G <- args$G
kk <- args$kk # industry x industry expenditures

# Upper bound is A and G + random edges.
upper_bound <- upper_bound(A,G,region_density,firm_density) 

# Lower bound is a sub-sample of true data A and G.
lower_bound <- lower_bound(A,G,region_sample_frac=1,firm_sample_frac=1)

x <- upper_bound
dim(x) <- c((R+N)*N,1)
nonzero_vars <- x %>% summary() %>% 
  tbl_df() %>% 
  rownames_to_column(var="i_hat") %>% 
  rename(i_original=i) %>%
  select(-j,-x) %>% 
  mutate(i_original=as.numeric(i_original),i_hat=as.numeric(i_hat))
#rmgc(x)

dim(lower_bound) <- c((R+N)*N,1)
x <- cbind(x,lower_bound)
rmgc(lower_bound)
penalty <- 1 - x[rowSums(x)>0,2]
rmgc(x,lower_bound)

## instead: should create that in X_ind. can calculate hmm...
X_ind <- create_X_ind(s=s,upper_bound=upper_bound,ik=ik)
# choose from potential non-zero variables that go into glmnet
X_ind <- X_ind[,nonzero_vars[["i_original"]]]

X_mc <- create_X_mc(I,s,upper_bound)
c_mc <- s

X_ag <- create_X_ag(upper_bound)
# X_ag %>% dim()
# X_ind %>% dim()
# X_mc %>% dim()

# xxx
c_mc <- s # RHS for market clearing equations
c_a <- rep_len(1,R) # RHS for rowSums(A) = 1 equations
c_g <- rep_len(1,N) #1-beta # RHS for rowSums(G) = 1-beta equations.
c_ind <- kk[,1]

# Apply it all together.
c <- c(c_mc,c_a,c_g,c_ind)
X <- rbind(X_mc,X_ag,X_ind)
rmgc(X_mc,X_ag,X_ind)
# deviance for some reason can go above 100%? I just want it to be as high as possible.
glmnet.control(devmax = 5) 

nlambda <- 100

# fit; returns coefs, and prediction.
system.time(fit <- fit_glmnet(X,c,alpha=1,nlambda=nlambda,lambda.min.ratio=1e-2))

# do analysis.

# nnz <- fit$df[min(nlam,length(fit$df))] / (R*N+N^2)
# nnzp <- fit_p$df[min(nlam,length(fit_p$df))] / (R*N+N^2)
# Merge back on the original parameter indices; e.g., go from the minimized y to the original y.

yx <- fit$coefs %>% left_join(nonzero_vars,by=c("i"="i_hat")) %>% select(-i,i=i_original) %>% df_to_s(dims=c((R+N)*N,1))

# y---first R are a_{.1}, next N are g_{.1}, etc. so split by R+N.
# this works like magic. must be filled bycolumn (e.g., equivalent of byrow=FALSE)
dim(yx) <- c(R+N,N)

# Get region-firm expenditure matrix A
A_hat <- yx[1:R,1:N]

# Get firm-firm expenditure share matrix G
G_hat <- yx[(R+1):(R+N),1:N]

rm(yx)

# Solve for implied s; if firm i has no customers, s_hati will be zero; this may require some post-processing
s_hatx <- t(A_hat) %*% (I %>% matrix(nrow=R,ncol=1)) + t(G_hat) %*% (s %>% matrix(nrow=N,ncol=1)) 

df <- tibble(s=s,s_hatx=s_hatx[,1])
ggplot(df %>% filter(s_hatx>0),aes(x=s,y=s_hatx))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
lm(s %>% log() ~ s_hatx %>% log(),data=df %>% filter(s_hatx>0)) %>% summary() # ok...
lm(c ~ fit$pred[,1]) %>% summary() # ok...

# 1:N, (N+1):(R+2*N), (R+2*N+1):(R+2*N+K^2)
lm(c ~ fit$pred[,1]) %>% summary() # ok...
lm(c_mc ~ fit$pred[1:N,1]) %>% summary() # ok...
# lm(c_a ~ fit$pred[(N+1):(R+N),1]) %>% summary() # ok...
# lm(c_g ~ fit$pred[(R+N+1):(R+2*N),1]) %>% summary() # ok...
lm(c_ind ~ fit$pred[(R+2*N+1):(R+2*N+K^2),1]) %>% summary() # ok...

fit$pred[(N+1):(R+N),1] %>% summary()
fit$pred[(R+N+1):(R+2*N),1] %>% summary()
ggplot() + geom_point(aes(x=c_ind,y=fit$pred[(R+2*N+1):(R+2*N+K^2),1])) #+ 
#  scale_x_log10() + scale_y_log10()

xxx
  
  X <- G %>% summary() %>% tbl_df() %>% filter(x>0.26) %>% df_to_s(dims=dim(G))
  X@x <- rep_len(1,length(X@x))
  library(igraph)
## uhh, so use A, G to do stuff?
  x <- graph_from_adjacency_matrix(adjmatrix=X, mode = c("directed"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
  
  co <- layout_with_fr(x)
  co[1,]
  plot(x, layout=co, vertex.color="black", vertex.size=1.75, edge.color="grey50", edge.arrow.size=0.05, vertex.label=NA) 

# s_hat <- pred[1:N,1]
# # s_hat_p <- pred_p[1:N,1]
# 
# # s <- s/sum(s)
# # s_hat <- s_hat/sum(s_hat)
# # s_hat_p <- s_hat_p/sum(s_hat_p)
# 
# a_eq <- pred[(N+1):(R+N),1]
# g_eq <- pred[(R+N+1):(R+2*N),1]
# 
# # compare pred to c.
# lm(c ~ pred[,1]) %>% summary() # ok...
# 
# lm(s ~ s_hat) %>% summary() # ok...
# lm(s[s_hat>0] %>% log() ~ s_hat[s_hat>0] %>% log()) %>% summary() # ok...
# lm(s ~ pred_p[1:N,1]) %>% summary() # ok...
# 
# # scaling. not sure why. something about relative importance of equations.
# 
# # ayy, these don't have to add up to 1 anymore, just I and s.
# pred_p[(N+1):(R+N),1] %>% summary()
# pred_p[(N+1):(R+N),1] %>% summary()
# a_eq %>% summary()
# g_eq %>% summary()
# # ((a_eq-I)/I) %>% summary()
# # ((g_eq-s)/s) %>% summary()
# 
# df <- tibble(s_hat=s_hat,s=s)
# # ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.1)
# ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()

# xxx
# ggplot(df %>% filter(s_hat>0),aes(x=s_hat,y=s_hat_p))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
# ggplot(df %>% rownames_to_column() %>% gather(type,value,s_hat:s) %>% filter(value>0))+stat_density(aes(x=value,colour=type),position="dodge",geom="line") + scale_x_log10()
# How sparse are A and G combined?
# sprintf("sparsity: %.5f",nnz)

# if you give upper bound only, sparsity doesn't seem to matter much (only about /2 or /3 of original data).
# if you give lower bound, it likes to hit it (setting half of vector to 0 doubles the remaining penalties, because glmnet scales penalty
# to nvar).

# yx <- y %>% summary() %>% tbl_df() %>% 
#   filter(j==max(j)) %>% # take the last iteration from glmnet
#   mutate(i=i-1,j=1) # remove intercept row (i=i-1), and add a column so I can convert to matrix (j=1)

xxx


## test solution:
# ad <- (A-A_hat)/A
# gd <- (G-G_hat)/G
# dim(ad) <- c(R*N,1)
# dim(gd) <- c(N*N,1)
# ad <- ad[,1]
# gd <- gd[,1]
# ad %>% summary()
# gd %>% summary()

# things to do: compare nnz in A and A_hat, and how it depends on sparsity
# same for G and G_hat

## True sparsity vs. estimated sparsity.
a_sp <- length(A@x) / (R*N)
g_sp <- length(G@x) / (N*N)

ah_sp <- length(A_hat@x) / (R*N)
gh_sp <- length(G_hat@x) / (N*N)

print(a_sp,5)
ah_sp %>% print(5)
g_sp %>% print(5)
gh_sp %>% print(5)


# maybe better to compare expenditure shares instead.
aexp <- (I %>% to_sdiag()) %*% A
aexphat <- (I %>% to_sdiag()) %*% A_hat
dim(aexp) <- c(R*N,1)
dim(aexphat) <- c(R*N,1)
aexp <- aexp[,1]
aexphat <- aexphat[,1]
aexp <- tibble(a=aexp,ah=aexphat)

# Remember that G needs beta too, but G_hat doesn't, I don't think.
gexp <- (s %>% to_sdiag()) %*% G
gexphat <- (s %>% to_sdiag()) %*% G_hat
dim(gexp) <- c(N*N,1)
dim(gexphat) <- c(N*N,1)
gexp <- gexp[,1]
gexphat <- gexphat[,1]
gexp <- tibble(a=gexp,ah=gexphat)


# true zeros, true non-zeros, false zeros, false positives
tz <- gexp %>% filter(ah==0 & a==0)
tnz <- gexp %>% filter(a>0 & ah>0)
fz <- gexp %>% filter(ah==0 & a>0)
fnz <- gexp %>% filter(ah>0 & a==0)

sens <- dim(tnz)[1] / (dim(tnz)[1] + dim(fz)[1])
spec <- dim(tz)[1] / (dim(tz)[1] + dim(fnz)[1])
posp <- dim(tnz)[1] / (dim(tnz)[1] + dim(fnz)[1])

sprintf("sparsity: %.5f, lower_bound: %.5f, actual: %.5f",nnz,nnzp,(length(A@x)+length(G@x))/(R*N+N^2))
sprintf("sens: %.5f, spec: %.5f, posp: %.5f",sens,spec,posp)

# F-score / information retrieval measure.
2 * (posp * sens) / (posp + sens)

lm(ah %>% log()~ a %>% log(),data=tnz) %>% summary()
ggplot(tnz,aes(x=a,y=ah)) + geom_point(alpha=0.01) + scale_x_log10() + scale_y_log10()
xxx







# Then compare non-zero values.
