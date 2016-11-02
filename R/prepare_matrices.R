########################################################################
# prepare_matrices.R
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
########################################################################

create_X_mc <- function(I,s,upper_bound) {
  z <- c(I,s)
  erg <- (z %>% to_sdiag()) %*% upper_bound
  erg <- erg %>% summary() %>% tbl_df()
  X_mc <- lapply(erg %>% split(f=erg$j), function(l) l$x ) %>% bdiag() %>% t()
  return(X_mc)
}

create_X_ag <- function(upper_bound) {
  dims <- dim(upper_bound)
  N <- dims[2]
  R <- dims[1] - N
  
  erg2 <- upper_bound %>% summary() %>% tbl_df()
  erg3 <- lapply(erg2 %>% split(f=erg2$j), function(l) {
    l %>% select(-j) %>% rownames_to_column(var="j") %>% mutate(j=as.integer(j)) %>% df_to_s(dims=c(R+N,dim(l)[1]))
  })

  X_ag <- c_binder(erg3)

  return(X_ag)
}

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
  
  # X_ind <- zz[[1]]
  # for (i in seq_along(zz[-1])) {
  #   X_ind <- cbind(X_ind, zz[[i+1]])
  # }
  X_ind <- c_binder(zz)

  return(X_ind)
}
