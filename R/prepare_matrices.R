########################################################################
# prepare_matrices.R
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
########################################################################

c_binder <- function(l) {
# new cbind function to avoid node stack overflow when cbinding N matrices together
# N could be thousands or tens-of-thousands or more. 
  
  N <- length(l)
  
  # recursively call lapply on half the list? no, too many.
  # divide into sqrt(N) lists of sqrt(N), plus leftovers.
  NN <- floor(N^(1/2))
  leftovers <- N-NN^2
  
  # do.call on a sublist of zz, determined by index i.
  # goes from NN*(i-1)+1 to NN*i
  g <- function(i,z) {
    do.call(cbind,z[(NN*(i-1)+1):(NN*i)] %>% unlist())
  }
  
  # cbind each of the elements of each sublist, to get a list
  # of cbinded elements. 
  w <- lapply(1:NN,g,z=l)
  
  # if leftovers, add them as well. 
  if (leftovers>0) {
    w <- list(w,l[(NN^2+1):N])
  }
  
  # cbind the sublists together.
  y <- do.call(cbind,w %>% unlist())
  return(y)
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
  
  # Save R and N. upper_bound is (R+N)xN
  dims <- dim(upper_bound)
  N <- dims[2]
  R <- dims[1] - N
  
  # Add region income and firm output to equations.
  x <- (c(I,(1-beta)*s) %>% to_sdiag()) %*% upper_bound
  
  # Convert that matrix to dataframe.
  x_df <- x %>% summary() %>% tbl_df()
  
  # Convert the dataframe to a list---otherwise it'll be too big to process.
  # Idea: split by columns, convert each column to a specific sparse matrix, then
  # cbind them all together afterward.
  x_list <- lapply(x_df %>% split(f=x_df$j), function(l) {
    l %>% select(-j) %>% # remove the column id.
      rownames_to_column(var="j") %>% # add new column ids that represent the parameters
      mutate(j=as.integer(j)) %>%
      df_to_s(dims=c(R+N,dim(l)[1])) # make matrix (R+N)x(nnz parameters for firm/original column id `j`)
  })
  rmgc(x_df) # attempt at memory management
  X_ag <- c_binder(x_list) # manually split them up and cbind them together to avoid node stack overflow
  rmgc(x_list)
  
  return(X_ag)
}

create_X_ind <- function(beta, s, upper_bound, ik) {
## Create the industry-pair expenditure equation matrix.
  
  # may need ((1-beta) * s) here instead. nah.
 
  # Get (NxN) diagonal size matrix, left multiply by firm-firm section of upper_bound matrix
  x <- (((1-beta) * s %>% to_sdiag()) %*% # get sizes, multiply by...
          upper_bound[(R+1):(R+N),1:N]) %>%  # ...firm part of upper bound
    summary() %>% tbl_df()

  # Get industry information. tbl_df() format to merge.
  k <- ik %>% summary() %>% tbl_df() %>% filter(x>0) %>% select(k=i,i=j,-x)
  
  # merge industry information onto firm-firm trade
  y <- x %>% left_join(k,by=c("i"="i")) %>% rename(ki=k) %>% left_join(k,by=c("j"="i")) %>% rename(kj=k)
  
  # region-firm upper bound data. use this as a placeholder. Later, need to add
  # industry-final demand IOT/SUT data. 
  rub <- upper_bound[1:R,1:N]

  f <- function(j,y) {
    temp <- j # producer/firm/column.
    
    # Take the j-th column of upper bound, reshape it into (KKxR) matrix. KK is the number of industry pairs,
    # R are the number of regions; this makes an empty KKxR matrix for each firm (because I don't have
    # the region-industry expenditures yet.
    v <- matrix(data=rub[,temp],nrow=K^2,ncol=R,byrow=TRUE) %>% as("sparseMatrix")
    v@x <- rep_len(0,length(v@x))
    
    # Pick the j-th "column" of y; producer information only.
    q <- y %>% filter(j==temp)
    
    # Convert it to a KKxN sparse matrix; 
    z <- with(q, sparseMatrix(i=((ki-1)*K+kj),j=i,x=x,dims=c(K^2,N)))
    
    # Add v (empty industry^2xregion matrix) and z (industry^2xfirm matrix) together; 
    w <- cbind(v,z)
    return(w)
  }
  
  # Start with y; for each producer/firm/column (1:N), apply f
  zz <- lapply(1:N,f,y=y)

  # Special cbind program to avoid node stack overflow. Sticks all elements of list together.  
  # This should be a (KKx(R+N)N) matrix. After I return it, I use nonzero_vars to remove 
  # parameters that aren't in upper_bound in order to pass it to glmnet with less parameters.
  X_ind <- c_binder(zz)

  return(X_ind)
}
