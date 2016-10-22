########################################################################
# main.R
# Clean and load relevant files
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
########################################################################

source("R/helpers.R")
startup()

## Initialize things.
#set.seed(9) 
R <- 2  # Number of regions
N <- 20 # Number of firms
Z <- 2  # Number of industries; each industry must have more than one firm, 
# or glmnet fails (at least until I add more equations).
region_density <- 0.5
firm_density <- 0.25

# Create fake data.
args <- initialize_fake_s(R=R,Z=Z,N=N,region_density=region_density,firm_density=firm_density)

# Save returned arguments.
beta <- args$beta
iz <- args$iz
I <- args$I[,1]
s <- sinit <- args$s
Er <- args$Er
En <- args$En

## so here try to add in info from Er, En.
## Get upper and lower bounds working. This is upper. Lower goes into the penalty.
## Which means: take this one, remove a bunch of links, then get penalty to work.

Erl <- Er %>% summary() %>% tbl_df() %>% sample_frac(size=0.25) %>% df_to_s(dims=c(R,N))
Enl <- En %>% summary() %>% tbl_df() %>% sample_frac(size=0.25) %>% df_to_s(dims=c(N,N))

lower_bound <- rbind(Erl,Enl)
upper_bound <- rbind(Er,En)
rm(Er,En)
dim(lower_bound) <- c((R+N)*N,1)

z <- c(I,s) 
erg <- (z %>% to_sdiag()) %*% penalty
erg <- erg %>% summary() %>% tbl_df()
X_mc <- lapply(erg %>% split(f=erg$j), function(l) l$x ) %>% bdiag() %>% t()
rm(erg)
c_mc <- s

erg2 <- penalty %>% summary() %>% tbl_df()

# another way: convert l to sparse diagonal, then elementwise multiply by sparseDiagonal, then remove need to elements errrg.
X_ag <- Reduce(cbind, lapply(erg2 %>% split(f=erg2$j), function(l) l %>% select(-j) %>% rownames_to_column(var="j") %>% mutate(j=as.integer(j)) %>% df_to_s(dims=c(R+N,dim(l)[1]))))
rm(erg2)
gc()

c_mc <- s # RHS for market clearing equations
c_a <- rep_len(1,R) # RHS for rowSums(A) = 1 equations
c_g <- rep_len(1,N) #1-beta # RHS for rowSums(G) = 1-beta equations.

# Apply it all together.
c <- c(c_mc,c_a,c_g)
X <- rbind(X_mc,X_ag)
rm(X_mc,X_ag)
# deviance for some reason can go above 100%? I just want it to be as high as possible.
glmnet.control(devmax = 5) 

# have a lot of lambdas. just want to get as close as possible to the data,
# don't need to worry about overfitting. it should exactly fit the "training" data.
#,lower.limits=0,upper.limits=1) <- may or may not need this.
fit <- glmnet(X,c,alpha=1,nlambda=50,intercept=FALSE,lower.limits=0,upper.limits=1,lambda.min.ratio=0.0001)
y <- coef(fit) # could make this more efficient by calling coef(fit,s="correct lambda here")
print(fit)
# xgboost could work if I predict a new value for every one of the million parameters etc. frig that.
pred <- predict(fit,newx=X,s=c(fit$lambda[50])) # some small lambda.

nnz <- fit$df[length(fit$df)] / (R*N+N^2)
rm(fit)
gc()

s_hat <- pred[1:N,1]
a_eq <- pred[(N+1):(R+N),1]
g_eq <- pred[(R+N+1):(R+2*N),1]

# compare pred to c.
lm(c ~ pred[,1]) %>% summary() # ok...
a_eq %>% summary()
g_eq %>% summary()
df <- tibble(s_hat=s_hat,s=s)
ggplot(df %>% filter(s_hat>0),aes(x=s,y=s_hat))+geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()
ggplot(df %>% rownames_to_column() %>% gather(type,value,s_hat:s) %>% filter(value>0))+stat_density(aes(x=value,colour=type),position="dodge",geom="line") + scale_x_log10()
# How sparse are A and G combined?
print(str_c("sparsity: ",nnz))

