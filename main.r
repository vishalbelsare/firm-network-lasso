########################################################################
# main.R
# Clean and load relevant files
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
# Notes: (1) it's important to pass around R and N to make sure the sparse matrix
# doesn't drop the last row or column if they're missing.
########################################################################
rm(list=ls(all=TRUE))
gc()

source("R/helpers.R")
startup()

## Initialize things.
# set.seed(8)

# Create fake data.
fakes <- TRUE
if (fakes) {
  print("Initialize fake data.")
  R <- 20  # Number of regions
  N <- 1000 # Number of firms
  K <- 10  # Number of industries; each industry must have more than one firm, 
  # or glmnet fails (at least until I add more equations).
  region_density <- 0.05
  firm_density <- 0.01
  scale <- N

#  args <- initialize_fake_s(R=R,K=K,N=N,region_density=region_density,firm_density=firm_density)
  args <- initialize_fake_economy(R=R,K=K,N=N,region_density=region_density,firm_density=firm_density,scale=scale)
} else {
  print("Load real data.")
  args <- load_data() 
}

# Save returned arguments.
I <- args$I[,1]
ik <- args$ik
kk <- args$kk
s <- args$s[,1]
beta <- args$beta

# Upper bound is A and G + random edges.
upper_bound <- upper_bound(args$A,args$G,region_density,firm_density) 

# Lower bound is a sub-sample of true data A and G.
lower_bound <- lower_bound(args$A,args$G,region_sample_frac=0.5,firm_sample_frac=0.1)

x <- upper_bound
dim(x) <- c((R+N)*N,1)
nonzero_vars <- x %>% summary() %>% 
  tbl_df() %>% 
  rownames_to_column(var="i_hat") %>% 
  rename(i_original=i) %>%
  select(-j,-x) %>% 
  mutate(i_original=as.numeric(i_original),i_hat=as.numeric(i_hat))

dim(lower_bound) <- c((R+N)*N,1)
x <- cbind(x,lower_bound)
rm(lower_bound)
penalty <- 1 - 1*x[rowSums(x)>0,2] # next, gonna try to change this to re-weight up firm expenditure relative to region expenditure (e.g., R/N or something)
rm(x)
gc()


# library(profvis)
# library(pryr)
## instead: should create that in X_ind. can calculate hmm...
# I think I Need to pass in nonzero_vars.
# profvis(X_ind <- create_X_ind(beta=beta,s=s,upper_bound=upper_bound,ik=ik,nonzero_vars=nonzero_vars))
X_ind <- create_X_ind(beta=beta,s=s,upper_bound=upper_bound,ik=ik,nonzero_vars=nonzero_vars)
# much better if I store it in transpose. but then I have to do that for all of them individually, then rbind instead of cbind, etc.

# choose from potential non-zero variables that go into glmnet.
# i_original are the column indices that are picked.
# X_ind <- X_ind[,nonzero_vars[["i_original"]]]

X_mc <- create_X_mc(I,beta,s,upper_bound)
X_ag <- create_X_ag(I,beta,s,upper_bound)

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

fit <- fit_glmnet(X,c,alpha=1,nlambda=nlambda,lambda.min.ratio=1e-12)
# fit <- fit_glmnet(X,c,alpha=1,nlambda=nlambda,lambda.min.ratio=1e-12,penalty=penalty)
pm <- fit %>% predicted_matrices(R=R,N=N,nonzero_vars=nonzero_vars)

prhs <- fit %>% predicted_rhs(R=R,N=N)

# these are so good that they aren't that useful anymore. but good for data.
# plot_rhs(prhs,var="mc",log=TRUE) # so output is all correct
# plot_rhs(prhs,var="g",log=TRUE) # but some firm-firm expenditure is zero,
# plot_rhs(prhs,var="a",log=FALSE)
# plot_rhs(prhs,var="ind",log=FALSE)

# Plot graph.
# plot_graph(G=pm$G_hat,edge_val_lower_bound=0.25)

fit_summary(prhs,c_mc,c_a,c_g,c_ind)
# fit_summary(prhs,c_mc,c_a,c_g,c_ind,stargazed=TRUE)

# Specificity, sensitivity.
# sensitivity_specificity(args$G,pm$G_hat) # input true and estimated matrices.

sparsity(fit)
# sparsity(fit,args$A,args$G)

# Q: region expenditures vs. firm expenditures? 

# plot rowSums(G) vs. total size?
rah <- rowSums(pm$A_hat)
rgh <- rowSums(pm$G_hat) #+beta
rah %>% summary()
rgh %>% summary()
ggplot() + geom_point(aes(x=s,y=rgh),alpha=0.25) + scale_x_log10(breaks=c(0.1,0.25,0.5,1,2,4,8)) + labs(x="Firm size",y="Sum of firm expenditure shares")
