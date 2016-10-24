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
#set.seed(9)

# Create fake data.
fakes <- TRUE
if (fakes) {
  print("Initialize fake data.")
  R <- 50  # Number of regions
  N <- 3000 # Number of firms
  K <- 2  # Number of industries; each industry must have more than one firm, 
  # or glmnet fails (at least until I add more equations).
  region_density <- 0.01
  firm_density <- 0.005

#  args <- initialize_fake_s(R=R,K=K,N=N,region_density=region_density,firm_density=firm_density)
  args <- initialize_fake_economy(R=R,K=K,N=N,region_density=region_density,firm_density=firm_density)
} else {
  print("Load real data.")
  
  # need: 
  # 1. list of plants. possibly export/import corrected, nah, for now, just get it to work.
  # 2. beta.
  # 3. regional incomes, I.
  # 4. Er, En, the highest possible ones, I hope. Hopefully only 10m, but maybe more.
  args <- load_data() # write this func.
}

# then calculate 

# Save returned arguments.
beta <- args$beta
ik <- args$ik
I <- args$I[,1]
s <- args$s[,1]
Er <- args$Er
En <- args$En

## Fake A and G
A <- args$A
G <- args$G

# 100, seems like an okay scaling? Not sure why. Number of expenditure share equations is R+N.
# Maybe it depends on that. R+N->1000 equations, all equal to 1. But the N sales equations
# if normalized, equal to 1. May also depend on sparsity. Also depends on skewness in data.
# Bigger plants need possibility of more connections; could also scale the normalization
# equations so that the RHS terms are N total outputs, and then R+N total expenditures.
I <- I/sum(s)*N*10
s <- s/sum(s)*N*10

# This is to test the other scaling bit.
# Erx <- (I %>% to_sdiag()) %*% Er
# Enx <- (s %>% to_sdiag()) %*% En

# Want to add potential non-zeros that are wrong. See if it gets close to sparsity.
# Do that in init func.
Erx <- Er
Enx <- En

x <- ((s %>% to_sdiag()) %*% En)  %>% summary() %>% tbl_df() 

# merge on industries
k <- ik %>% summary() %>% tbl_df() %>% select(k=i,i=j,-x)
y <- x %>% left_join(k,by=c("i"="i")) %>% rename(ki=k) %>% left_join(k,by=c("j"="i")) %>% rename(kj=k)
X_ind <- with(y, sparseMatrix(i=((ki-1)*K+kj),j=(i-1)*N+i*R+j,x=x,dims=c(K^2,(R+N)*N)))
#with(y, sparseMatrix(i=((ki-1)*K+kj),j=(i-1)*N+j,x=1,dims=c(K^2,(N)*N)))

# still need to multiply by s. and get RHS things.
# elementwise multiply by some weird s thing. or merge it directly onto y or x.
# RHS: 


## so here try to add in info from Er, En.
## Get upper and lower bounds working. This is upper. Lower goes into the penalty.
## Which means: take this one, remove a bunch of links, then get penalty to work.

Erl <- A %>% summary() %>% tbl_df() %>% mutate(x=1) %>% sample_frac(size=1) %>% df_to_s(dims=c(R,N))
Enl <- G %>% summary() %>% tbl_df() %>% mutate(x=1) %>% sample_frac(size=1) %>% df_to_s(dims=c(N,N))

lower_bound <- rbind(Erl,Enl)
upper_bound <- rbind(Erx,Enx)

x <- upper_bound
dim(x) <- c((R+N)*N,1)

nonzero_vars <- x %>% summary() %>% 
  tbl_df() %>% 
  rownames_to_column(var="i_hat") %>% 
  rename(i_original=i) %>%
  select(-j,-x) %>% 
  mutate(i_original=as.numeric(i_original),i_hat=as.numeric(i_hat))
rm(x)

x <- upper_bound
dim(x) <- c((R+N)*N,1) # have to reduce this too.
dim(lower_bound) <- c((R+N)*N,1) # have to reduce this too.

x <- cbind(x,lower_bound)
penalty <- 1 - x[rowSums(x)>0,2]
rm(x,lower_bound)
rm(Er,En,Erx,Enx)

z <- c(I,s) 
erg <- (z %>% to_sdiag()) %*% upper_bound
erg <- erg %>% summary() %>% tbl_df()
X_mc <- lapply(erg %>% split(f=erg$j), function(l) l$x ) %>% bdiag() %>% t()
rm(erg)
gc()
c_mc <- s

erg2 <- upper_bound %>% summary() %>% tbl_df()
curr <- 0
erg3 <- lapply(erg2 %>% split(f=erg2$j), function(l) {
  # print(str_c("j: ",l[1,"j"]," dim(l): ",dim(l)[1]))
  l %>% select(-j) %>% rownames_to_column(var="j") %>% mutate(j=as.integer(j)) %>% df_to_s(dims=c(R+N,dim(l)[1]))
})


# Reduce takes way too long.
# another way: convert l to sparse diagonal, then elementwise multiply by sparseDiagonal, then remove need to elements errrg.
# X_ag <- Reduce(cbind, lapply(erg2 %>% split(f=erg2$j), function(l) {
#     print(l[1,"j"])
#     l %>% select(-j) %>% rownames_to_column(var="j") %>% mutate(j=as.integer(j)) %>% df_to_s(dims=c(R+N,dim(l)[1]))
#   }))

#X_ag <- erg3[[1]]
n_per <- 100 # cbind won't take the whole list at once. divide into 100 separate ~300 length lists.
NN <- floor(N/n_per)
# divide into 300 100 length lists?
erg4 <- lapply(1:NN, function(i) {
  # print(i)
  do.call(cbind,erg3[((i-1)*100+1):(i*100)])
})
erg4.5 <- do.call(cbind,erg3[(NN*100+1):N])
erg5 <- do.call(cbind,erg4) #,cbind(erg3[(NN*100+1):N]))
X_ag <- cbind(erg5,erg4.5)

rm(n_per,NN,erg4,erg5,erg4.5)
rm(erg2,erg3)
gc()
c_mc <- s # RHS for market clearing equations
# c_a <- I #rep_len(1,R) # RHS for rowSums(A) = 1 equations
# c_g <- s #rep_len(1,N) #1-beta # RHS for rowSums(G) = 1-beta equations.
c_a <- rep_len(1,R) # RHS for rowSums(A) = 1 equations
c_g <- rep_len(1,N) #1-beta # RHS for rowSums(G) = 1-beta equations.

# Apply it all together.
c <- c(c_mc,c_a,c_g)
X <- rbind(X_mc,X_ag)
rm(X_mc,X_ag)
# deviance for some reason can go above 100%? I just want it to be as high as possible.
glmnet.control(devmax = 2) 

nlam <- 100
# have a lot of lambdas. just want to get as close as possible to the data,
# don't need to worry about overfitting. it should exactly fit the "training" data.
#,lower.limits=0,upper.limits=1) <- may or may not need this.
fit <- glmnet(X,c,alpha=1,nlambda=nlam,intercept=FALSE,lower.limits=0,upper.limits=1,lambda.min.ratio=1e-6)
fit_p <- glmnet(X,c,alpha=1,nlambda=nlam,intercept=FALSE,lower.limits=0,upper.limits=1,lambda.min.ratio=1e-6,penalty.factor=penalty)
y <- coef(fit) # could make this more efficient by calling coef(fit,s="correct lambda here")
y_p <- coef(fit_p) # could make this more efficient by calling coef(fit,s="correct lambda here")
print(fit)
# xgboost could work if I predict a new value for every one of the million parameters etc. frig that.
pred <- predict(fit,newx=X,s=c(fit$lambda[min(length(fit$df),nlam)])) # some small lambda.
pred_p <- predict(fit_p,newx=X,s=c(fit$lambda[min(length(fit$df),nlam)])) # some small lambda.

nnz <- fit$df[min(nlam,length(fit$df))] / (R*N+N^2)
nnzp <- fit_p$df[min(nlam,length(fit_p$df))] / (R*N+N^2)
# rm(fit,fit_p)
gc()

s_hat <- pred[1:N,1]
# s_hat_p <- pred_p[1:N,1]

# s <- s/sum(s)
# s_hat <- s_hat/sum(s_hat)
# s_hat_p <- s_hat_p/sum(s_hat_p)

a_eq <- pred[(N+1):(R+N),1]
g_eq <- pred[(R+N+1):(R+2*N),1]

# compare pred to c.
lm(c ~ pred[,1]) %>% summary() # ok...

lm(s ~ s_hat) %>% summary() # ok...
lm(s[s_hat>0] %>% log() ~ s_hat[s_hat>0] %>% log()) %>% summary() # ok...
lm(s ~ pred_p[1:N,1]) %>% summary() # ok...

# scaling. not sure why. something about relative importance of equations.

# ayy, these don't have to add up to 1 anymore, just I and s.

pred_p[(N+1):(R+N),1] %>% summary()
a_eq %>% summary()
g_eq %>% summary()
# ((a_eq-I)/I) %>% summary()
# ((g_eq-s)/s) %>% summary()

df <- tibble(s_hat=s_hat,s=s)
# ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.1)
ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()

# ggplot(df %>% filter(s_hat>0),aes(x=s_hat,y=s_hat_p))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
# ggplot(df %>% rownames_to_column() %>% gather(type,value,s_hat:s) %>% filter(value>0))+stat_density(aes(x=value,colour=type),position="dodge",geom="line") + scale_x_log10()
# How sparse are A and G combined?
# sprintf("sparsity: %.5f",nnz)
sprintf("sparsity: %.5f, lower_bound: %.5f, actual: %.5f",nnz,nnzp,(length(A@x)+length(G@x))/(R*N+N^2))

# if you give upper bound only, sparsity doesn't seem to matter much (only about /2 or /3 of original data).
# if you give lower bound, it likes to hit it (setting half of vector to 0 doubles the remaining penalties, because glmnet scales penalty
# to nvar).

xxx
yx <- y %>% summary() %>% tbl_df() %>% 
  filter(j==max(j)) %>% # take the last iteration from glmnet
  mutate(i=i-1,j=1) # remove intercept row (i=i-1), and add a column so I can convert to matrix (j=1)

# Merge back on the original parameter indices; e.g., go from the minimized y to the original y.
yx <- yx %>% left_join(nonzero_vars,by=c("i"="i_hat")) %>% select(-i,i=i_original)

# Convert to sparse matrix.
yx <- yx %>% df_to_s(dims=c((R+N)*N,1))

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
#lm(s %>% log() ~ s_hatx[,1] %>% log()) %>% summary() # ok...
df <- tibble(s_hat=s_hat,s=s,s_hatx=s_hatx[,1])
ggplot(df %>% filter(s_hatx>0),aes(x=s,y=s_hatx))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()

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

a_sp
ah_sp
g_sp
gh_sp


xxx

# maybe better to compare expenditure shares instead.
aexp <- (I %>% to_sdiag()) %*% A
aexphat <- (I %>% to_sdiag()) %*% A_hat
dim(aexp) <- c(R*N,1)
dim(aexphat) <- c(R*N,1)
aexp <- aexp[,1]
aexphat <- aexphat[,1]
aexp <- tibble(a=aexp,ah=aexphat)
aexp %>% View()

xxx

lm(aexp ~ aexphat) %>% summary()







# Then compare non-zero values.
