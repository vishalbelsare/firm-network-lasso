########################################################################
# prepare_matrices.R
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
########################################################################
r_binder <- function(l,r_bind=TRUE) {
  # r_bind -> rbind the matrices all together. 
  # !r_bind -> cbind the matrices all together. could use eval(parse) type code to eliminate if-statements?
  
  # new rbind function to avoid node stack overflow when rbinding N matrices together
  # N could be thousands or tens-of-thousands or more. same as c_binder, but transpose
  # for sparse matrices that are better stored in row form.
  N <- length(l)
  
  # recursively call lapply on half the list? no, too many.
  # divide into sqrt(N) lists of sqrt(N), plus leftovers.
  NN <- floor(N^(1/2))
  leftovers <- N-NN^2
  
  # do.call on a sublist of zz, determined by index i.
  # goes from NN*(i-1)+1 to NN*i
  if (r_bind) {
    g <- function(z) {
      library(Matrix)
      do.call(rbind,unlist(z))
    }
  } else if (!r_bind) {
    g <- function(z) {
      library(Matrix)
      do.call(cbind,unlist(z))
    }
  }

  # this is used to split the list into pieces; need a vector that says
  # 1,1,...,1 for the first NN, then 2,2,..,2, etc etc. 
  # Then split() uses that to split the given list l into NN+1 pieces, 
  # so that the first NN elements are in the first new list element, and so on.
  # then we pass each individual list into the g function, which rbinds it.
  x <- matrix(data=1:NN,nrow=NN,ncol=NN,byrow=TRUE)
  dim(x) <- c(NN^2,1)
  x <- c(x[,1],rep_len(NN+1,N-NN^2))
  z <- split(unlist(l),x)

  cl <- makeCluster(2)
  w <- parLapply(cl,z,g)
  stopCluster(cl)
  # we get back `w`, which is a list of rbinded matrices that need to be rbinded themselves.

  # rbind the sublists together.
  if (r_bind) {
    y <- do.call(rbind,unlist(w))
  } else if (!r_bind) {
    y <- do.call(cbind,unlist(w))
  }
  return(y)
}

c_binder <- function(l) {
  return(r_binder(l,r_bind=FALSE))
}

create_X_mc <- function(I,beta,s,upper_bound) {
## Create region-firm and firm-firm market clearing equations.
  w <- (c(I,(1-beta)*s) %>% to_sdiag()) %*% upper_bound # add income and size to equations.
  w <- w %>% summary() %>% tbl_df() # convert to dataframe
  
  # Idea: turn the matrix into a list of columns; use bdiag() to make a block diagonal matrix,
  # where each non-zero bit is c(I,s); need to transpose after, since it goes by column 
  # (since each part of the list is a column of upper_bound, I think)
  X_mc <- lapply(w %>% split(f=w$j), function(l) l$x ) %>% bdiag() %>% t()
  return(X_mc)
}

create_X_ag <- function(I,beta,s,upper_bound) {
## Create the region-firm and firm-firm expenditure matrix.
  
  library(parallel)
  
  # Save R and N. upper_bound is (R+N)xN
  dims <- dim(upper_bound)
  N <- dims[2]
  R <- dims[1] - N
  
  # Add region income and firm output to equations.
  x <- (c(I,(1-beta)*s) %>% to_sdiag()) %*% upper_bound
  
  # Convert that matrix to dataframe.
  x <- x %>% summary() %>% tbl_df()

  l <- x %>% split(x$j)
  rm(x) # attempt at memory management
  
  f <- function(l,R,N) {
    library(Matrix)
    library(dplyr)
    l %>% select(-j) %>% # remove the column id.
      tibble::rownames_to_column(var="j") %>% # add new column ids that represent the parameters
      mutate(j=as.integer(j)) %>%
      df_to_s(dims=c(R+N,dim(l)[1])) # make matrix (R+N)x(nnz parameters for firm/original column id `j`)
  }

  cl <- makeCluster(2)
  clusterExport(cl,c('df_to_s'),envir=environment())
  x_list <- parLapply(cl,l,fun=f,R=R,N=N)
  stopCluster(cl)
    
  X_ag <- c_binder(x_list) # manually split them up and cbind them together to avoid node stack overflow
  rm(x_list)
  
  return(X_ag)
}

create_X_ind <- function(beta, s, upper_bound, ik, nonzero_vars) {
## Create the industry-pair expenditure equation matrix.
  
  library(parallel)
  # may need ((1-beta) * s) here instead. nah.
 
  # Get (NxN) diagonal size matrix, left multiply by firm-firm section of upper_bound matrix
  x <- (((1-beta) * s %>% to_sdiag()) %*% # get sizes, multiply by...
          upper_bound[(R+1):(R+N),1:N]) %>%  # ...firm part of upper bound
    summary() %>% tbl_df()
  
  # Get industry information. tbl_df() format to merge.
  k <- ik %>% summary() %>% tbl_df() %>% filter(x>0) %>% select(k=i,i=j,-x)
  
  # merge industry information onto firm-firm trade
  y <- x %>% left_join(k,by=c("i"="i")) %>% rename(ki=k) %>% left_join(k,by=c("j"="i")) %>% rename(kj=k)
  
  # h <- y %>% split(y$j)
  # rm(y)

  f <- function(q,K,R,N) {
    library(Matrix)
    # Take the j-th column of upper bound, reshape it into (KKxR) matrix. KK is the number of industry pairs,
    # R are the number of regions; this makes an empty KKxR matrix for each firm (because I don't have
    # the region-industry expenditures yet.
    
    v <- sparseMatrix(i=NULL,j=NULL,x=0,dims=c(K^2,R))

    # Convert it to a KKxN sparse matrix; 
    z <- with(q, sparseMatrix(i=((ki-1)*K+kj),j=i,x=x,dims=c(K^2,N)))
    
    # Add v (empty industry^2xregion matrix) and z (industry^2xfirm matrix) together; 
    w <- cbind(v,z)
    return(t(w)) # transpose; sparse matrix stores things using column-orientation.
    # the w / X_ind matrices are large if stored that way---this stores them in row-orientation, then
    # removes the extra bits that aren't in nonzero_vars, then transposes back before it goes into the matrix.
    # This also means you need to rbind these together, instead of cbinding.
  }
  
  # Start with y; for each producer/firm/column (1:N), apply f
  cl <- makeCluster(2)
  zz <- parLapply(cl,y %>% split(y$j),fun=f,K=K,R=R,N=N)
  stopCluster(cl)

  # Special cbind program to avoid node stack overflow. Sticks all elements of list together.  
  # This should be a (KKx(R+N)N) matrix. I use nonzero_vars to remove 
  # parameters that aren't in upper_bound in order to pass it to glmnet with less parameters.
  X_ind <- r_binder(zz)
  X_ind <- X_ind[nonzero_vars[["i_original"]],] %>% t()
  return(X_ind)
}
