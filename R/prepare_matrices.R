########################################################################
# prepare_matrices.R
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
########################################################################

c_binder <- function(l) {
  N <- length(l)
  
  # recursively call lapply on half the list? no, too many.
  # divide into sqrt(N) lists of sqrt(N), plus leftovers.
  NN <- floor(N^(1/2))
  leftovers <- N-NN^2
  
  g <- function(i,zz) {
    do.call(cbind,zz[(NN*(i-1)+1):(NN*i)] %>% unlist())
  }
  
  zzz <- lapply(1:NN,g,zz=l)
  if (leftovers>0) {
    zzz <- list(zzz,l[(NN^2+1):N])
  }
  zzzz <- do.call(cbind,zzz %>% unlist())
  return(zzzz)
}

create_X_mc <- function(I,s,upper_bound) {
  z <- c(I,s)
  w <- (z %>% to_sdiag()) %*% upper_bound
  w <- w %>% summary() %>% tbl_df()
  X_mc <- lapply(w %>% split(f=w$j), function(l) l$x ) %>% bdiag() %>% t()
  return(X_mc)
}

create_X_ag <- function(I,s,upper_bound) {
  dims <- dim(upper_bound)
  N <- dims[2]
  R <- dims[1] - N
  
  x <- (c(I,s) %>% to_sdiag()) %*% upper_bound
  
  x_df <- x %>% summary() %>% tbl_df()
  x_list <- lapply(x_df %>% split(f=x_df$j), function(l) {
    l %>% select(-j) %>% # remove the column id.
      rownames_to_column(var="j") %>% # add new column ids that represent the parameters
      mutate(j=as.integer(j)) %>%
      df_to_s(dims=c(R+N,dim(l)[1])) # make matrix (R+N)x(nnz parameters for firm/original column id `j`)
  })
  rmgc(x_df)
  X_ag <- c_binder(x_list)
  rmgc(x_list)
  
  return(X_ag)
}

create_X_ind <- function(s, upper_bound, ik) {
  
  # may need ((1-beta) * s) here instead. but for now, the estimates are of g_hat = (1-beta)*g
  x <- ((s %>% to_sdiag()) %*% upper_bound[(R+1):(R+N),1:N]) %>% summary() %>% tbl_df()

  k <- ik %>% summary() %>% tbl_df() %>% filter(x>0) %>% select(k=i,i=j,-x)
  y <- x %>% left_join(k,by=c("i"="i")) %>% rename(ki=k) %>% left_join(k,by=c("j"="i")) %>% rename(kj=k)
  
  rub <- upper_bound[1:R,1:N]

  f <- function(j,y) {
    temp <- j
    v <- matrix(data=rub[,temp],nrow=K^2,ncol=R,byrow=TRUE) %>% as("sparseMatrix")
    v@x <- rep_len(0,length(v@x))
    q <- y %>% filter(j==temp)
    z <- with(q, sparseMatrix(i=((ki-1)*K+kj),j=i,x=x,dims=c(K^2,N)))
    w <- cbind(v,z)
    return(w)
  }
  zz <- lapply(1:N,f,y=y)
  
  X_ind <- c_binder(zz)

  return(X_ind)
}
