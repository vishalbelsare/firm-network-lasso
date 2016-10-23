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
  N <- 1000 # Number of firms
  Z <- 2  # Number of industries; each industry must have more than one firm, 
  # or glmnet fails (at least until I add more equations).
  region_density <- 0.5
  firm_density <- 0.25

  args <- initialize_fake_s(R=R,Z=Z,N=N,region_density=region_density,firm_density=firm_density)
  
} else {
  print("Load real data.")
  
  # need: 
  # 1. list of plants. possibly export/import corrected, nah, for now, just get it to work.
  # 2. beta.
  # 3. regional incomes, I.
  # 4. Er, En, the highest possible ones, I hope. Hopefully only 10m, but maybe more.
  args <- load_data() # write this func.
  
  # and that's it? 
  
}


# Save returned arguments.
beta <- args$beta
iz <- args$iz
I <- args$I[,1]
s <- sinit <- args$s
Er <- args$Er
En <- args$En

# 100, seems like an okay scaling? Not sure why. Number of expenditure share equations is R+N.
# Maybe it depends on that. R+N->1000 equations, all equal to 1. But the N sales equations
# if normalized, equal to 1. May also depend on sparsity. Also depends on skewness in data.
# Bigger plants need possibility of more connections; could also scale the normalization
# equations so that the RHS terms are N total outputs, and then R+N total expenditures.
I <- I/sum(s)*N
s <- s/sum(s)*N

# This is to test the other scaling bit.
# Erx <- (I %>% to_sdiag()) %*% Er
# Enx <- (s %>% to_sdiag()) %*% En
Erx <- Er
Enx <- En

## so here try to add in info from Er, En.
## Get upper and lower bounds working. This is upper. Lower goes into the penalty.
## Which means: take this one, remove a bunch of links, then get penalty to work.

Erl <- Er %>% summary() %>% tbl_df() %>% sample_frac(size=1) %>% df_to_s(dims=c(R,N))
Enl <- En %>% summary() %>% tbl_df() %>% sample_frac(size=0.1) %>% df_to_s(dims=c(N,N))

lower_bound <- rbind(Erl,Enl)
upper_bound <- rbind(Erx,Enx)

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

# another way: convert l to sparse diagonal, then elementwise multiply by sparseDiagonal, then remove need to elements errrg.
X_ag <- Reduce(cbind, lapply(erg2 %>% split(f=erg2$j), function(l) l %>% select(-j) %>% rownames_to_column(var="j") %>% mutate(j=as.integer(j)) %>% df_to_s(dims=c(R+N,dim(l)[1]))))
rm(erg2)
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
#glmnet.control(devmax = 5) 

nlam <- 200
# have a lot of lambdas. just want to get as close as possible to the data,
# don't need to worry about overfitting. it should exactly fit the "training" data.
#,lower.limits=0,upper.limits=1) <- may or may not need this.
fit <- glmnet(X,c,alpha=0.1,nlambda=nlam,intercept=FALSE,lower.limits=0,upper.limits=1,lambda.min.ratio=1e-16)
# fit_p <- glmnet(X,c,alpha=1,nlambda=nlam,intercept=FALSE,lower.limits=0,upper.limits=1,lambda.min.ratio=0.0000001,penalty.factor=penalty)
#y <- coef(fit) # could make this more efficient by calling coef(fit,s="correct lambda here")
print(fit)
# xgboost could work if I predict a new value for every one of the million parameters etc. frig that.
pred <- predict(fit,newx=X,s=c(fit$lambda[nlam])) # some small lambda.
# pred_p <- predict(fit_p,newx=X,s=c(fit$lambda[nlam])) # some small lambda.

nnz <- fit$df[length(fit$df)] / (R*N+N^2)
# nnzp <- fit_p$df[length(fit_p$df)] / (R*N+N^2)
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

# scaling. not sure why. something about relative importance of equations.

# ayy, these don't have to add up to 1 anymore, just I and s.

# pred_p[(N+1):(R+N),1] %>% summary()
a_eq %>% summary()
g_eq %>% summary()
# ((a_eq-I)/I) %>% summary()
# ((g_eq-s)/s) %>% summary()

df <- tibble(s_hat=s_hat,s=s)
ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.1)
ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
xxx

# ggplot(df %>% filter(s_hat>0),aes(x=s_hat,y=s_hat_p))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
ggplot(df %>% rownames_to_column() %>% gather(type,value,s_hat:s) %>% filter(value>0))+stat_density(aes(x=value,colour=type),position="dodge",geom="line") + scale_x_log10()
# How sparse are A and G combined?
sprintf("sparsity: %.5f",nnz)
# sprintf("sparsity: %.5f, lower_bound: %.5f",nnz,nnzp)

