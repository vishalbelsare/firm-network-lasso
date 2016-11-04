########################################################################
# fit_analysis.R
# License: MIT

# ""
# Jesse Tweedle
# Oct, 2016
########################################################################

fit_glmnet <- function(X,c,alpha=1,nlambda=100,lambda.min.ratio=1e-5,penalty=rep_len(1,dim(X)[2])) {
  fit <- glmnet(X,c,alpha=1,nlambda=nlambda,intercept=FALSE,lower.limits=0,upper.limits=1,lambda.min.ratio=lambda.min.ratio)
  
  # pick last iteration.
  t = min(length(fit$df),nlambda)
  
  y <- coef(fit) # could make this more efficient by calling coef(fit,s="correct lambda here")
  pred <- predict(fit,newx=X,s=c(fit$lambda[t])) # some small lambda.
  #df <- fit$df[t] #/ (R*N+N^2)
  
  coefs <- y %>% summary() %>% tbl_df() %>% 
    filter(j==t) %>% # take the last iteration from glmnet
    mutate(i=i-1,j=1) # remove intercept row (i=i-1), and add a column so I can convert to matrix (j=1)
  
  
  return(list(coefs=coefs,pred=pred))
}

predicted_matrices <- function(fit,R,N,nonzero_vars) {
  yx <- fit$coefs %>% left_join(nonzero_vars,by=c("i"="i_hat")) %>% select(-i,i=i_original) %>% df_to_s(dims=c((R+N)*N,1))
  
  # y---first R are a_{.1}, next N are g_{.1}, etc. so split by R+N.
  # this works like magic. must be filled bycolumn (e.g., equivalent of byrow=FALSE)
  dim(yx) <- c(R+N,N)
  
  # Get region-firm expenditure matrix A
  A_hat <- yx[1:R,1:N]
  
  # Get firm-firm expenditure share matrix G
  G_hat <- yx[(R+1):(R+N),1:N]
  
  return(list(A_hat=A_hat,G_hat=G_hat))
}

predicted_rhs <- function(fit,R,N) {
  # divided up like --- 1:N, (N+1):(R+2*N), (R+2*N+1):(R+2*N+K^2)

  c_mc_hat <- fit$pred[1:N,1] # predicted sizes
  c_a_hat <- fit$pred[(N+1):(R+N),1] # predicted region-final demand sums
  c_g_hat <- fit$pred[(R+N+1):(R+2*N),1] # predicted firm demand sums
  c_ind_hat <- fit$pred[(R+2*N+1):(R+2*N+K^2),1] # predicted industry-pair expenditures
  return(list(c_mc_hat=c_mc_hat,c_a_hat=c_a_hat,c_g_hat=c_g_hat,c_ind_hat=c_ind_hat))
}

fit_summary <- function(predicted_rhs,c_mc,c_a,c_g,c_ind) {
  # 1:N, (N+1):(R+2*N), (R+2*N+1):(R+2*N+K^2)

  print("Overall accuracy")
  lm(c(c_mc,c_a,c_g,c_ind) ~ Reduce(c,predicted_rhs)) %>% summary() %>% print()
  
  print("Firm size accuracy / market clearing equations")
  lm(c_mc ~ predicted_rhs$c_mc_hat) %>% summary() %>% print()

  print("Normalization equations")
  lm(c(c_a,c_g) ~ c(predicted_rhs$c_a_hat,predicted_rhs$c_g_hat)) %>% summary() %>% print()
  
  
  print("Industry expenditure equations")
  lm(c_ind ~ predicted_rhs$c_ind_hat) %>% summary() # ok...
}

sparsity <- function(fit,A_hat=NULL,G_hat=NULL) {
  # print nnz, etc etc.  
  # well, not really.

  nnz <- fit$coefs %>% length() / (R*N + N*N)
  
  if (!is.null(A_hat) & !is.null(G_hat)) {
    sprintf("sparsity: %.5f, actual: %.5f",nnz,(length(A@x)+length(G@x))/(R*N+N^2))
  } else {
    sprintf("sparsity: %.5f",nnz)
  }
}

sensitivity_specificity <- function(M, Mhat) { 
  # M, Mhat must be sparse matrices.
  # only useful for fake data. in real data, we don't know true zeros / non-zeros.

  size <- Reduce(`*`,dim(M))

  M <- M %>% summary() %>% tbl_df() %>% rename(m=x)
  Mhat <- Mhat %>% summary() %>% tbl_df() %>% rename(mhat=x)
  
  mm <- full_join(M,Mhat)
  
  # ahh, yeah, just use NA. and size
  # true zeros, true non-zeros, false zeros, false positives
  
  # A true zero isn't in the dataset. 
  true_zeros <- size - dim(mm)[1]
  
  # A true non-zero is when both are non-zero.
  true_nonzeros <- dim(mm %>% filter(m>0 & mhat>0 & !is.na(mhat) & !is.na(m)))[1]

  # A false zero is a zero/NA in mhat when it's non zero in the original.  
  false_zeros <- dim(mm %>% filter((mhat==0 | is.na(mhat)) & m>0))[1]
  
  # false non-zero is when mhat is positive but m is 0/NA
  false_nonzeros <- dim(mm %>% filter(mhat>0 & (m==0 | is.na(m))))[1]
  
  sensitivity <- true_nonzeros / (true_nonzeros + false_zeros) 
  specificity <- true_zeros / (true_zeros + false_nonzeros) 
  positive_predictive <- true_nonzeros / (true_nonzeros + false_nonzeros)
  sprintf("sensitivity: %.5f, specificity: %.5f, positive predictive value: %.5f",sensitivity,specificity,positive_predictive)
  
  return(list(true_zeros=true_zeros,true_nonzeros=true_nonzeros,false_zeros=false_zeros,false_nonzeros=false_nonzeros,sensitivity=sensitivity,specificity=specificity,positive_predictive=positive_predictive))
}

plot_rhs <- function(prhs,var,log=FALSE) {
  if (var=="mc") {
    lab="Firm output"
  } else if (var=="ind") {
    lab="Industry pair exp."
  } else if (var=="a") {
    lab="Region expenditure"
  } else if (var=="g") {
    lab="Firm expenditure"
  }
  
  x <- str_c("c_",var) %>% as.name() %>% eval()
  y <- prhs[[str_c("c_",var,"_hat")]]

  p <- ggplot() + geom_point(aes(
    x=x,
    y=y
  )) + labs(x=lab,y=str_c("Predicted ", lab))
  if (log) {
    lab <- str_c(lab," (log scale)")
    p <- p + scale_x_log10() + scale_y_log10() + labs(x=lab,y=str_c("Predicted ", lab))
  }
  plot(p)
}

# Figures:
plot_graph <- function(G,edge_val_lower_bound) { # G is a symmetric firm-firm trade matrix
  X <- G %>% summary() %>% tbl_df() %>% filter(x>edge_val_lower_bound) %>% df_to_s(dims=dim(G))
  X@x <- rep_len(1,length(X@x))
  library(igraph)
  x <- graph_from_adjacency_matrix(adjmatrix=X, mode = c("directed"), weighted = NULL, diag = TRUE,add.colnames = NULL, add.rownames = NA)
  
  co <- layout_with_fr(x)
  plot(x, layout=co, vertex.color="black", vertex.size=1.75, edge.color="grey50", edge.width=0.01, edge.arrow.size=0.05, vertex.label=NA) 
}