########################################################################
# initialize_functions.R
# Functions to generate fake parameters and data
# License: MIT

# ""
# Jesse Tweedle
# Oct, 2016
########################################################################

initialize_fake_economy <- function(R,K,N,region_density,firm_density,scale=N) {
  
  # Plant value-added share
  beta <- runif(N,0.1,0.9)
  
  # Plant's region
  ir <- rbind(tibble(j=1:N,i=sample(1:R,N,replace=TRUE),x=1)) %>% df_to_s(dims=c(R,N))
  
  # Plant's industry, total of J.
  ik <- tibble(j=1:N,i=sample(1:K,N,replace=TRUE),x=1) %>% df_to_s(dims=c(K,N))
  # ik <- tibble(j=1:N,i=rep(1:K,ceiling(N/K))[1:N],x=1) %>% df_to_s(dims=c(K,N))
  
  Er <- rsparsematrix(R,N,density=region_density,rand.x=function(n) 1)
  # needs to be at least >=ir
  # Er <- Er + ir
  # Er@x <- rep_len(1,length(Er@x))

  while (min(colSums(Er))==0 | min(rowSums(Er))==0) {
    if (min(colSums(Er))==0) {
      q <- colSums(Er) == 0
      Er[,q] <- tibble(j=1:length(q[q]),
                       i=sample(1:R,length(q[q]),replace=TRUE),
                       x=1) %>%
        df_to_s(dims=c(R,length(q[q])))
    }
    if (min(rowSums(Er))==0) {
      q <- rowSums(Er) == 0
      Er[q,] <- tibble(i=1:length(q[q]),
                       j=sample(1:N,length(q[q]),replace=TRUE),
                       x=1) %>%
        df_to_s(dims=c(length(q[q]),N))
    }
    Er <- Er %>% summary() %>% tbl_df() %>% filter(x>0) %>% df_to_s(dims=c(R,N))
  }
  
  A <- Er
  A@x <- rlnorm(length(A@x),0,1) 
  # A@x <- runif(length(A@x)) 
  A <- A / rowSums(A)

  # Non-zero edges of plant-plant demand matrix
  # Possibly add services demand to each one that doesn't have a thing.
  En <- rsparsematrix(N,N,density=firm_density,rand.x=function(n) 1)

  while ( min(colSums(En))==0 | min(rowSums(En))==0 | sum(diag(En[1:N,1:N]))>0 ) {
    diag(En[1:N,1:N]) <- 0
    # for columns:
    if (min(colSums(En))==0) {
      q <- colSums(En) == 0
      
      En[,q] <- tibble(j=1:length(q[q]),
                       i=1+which(q),
                       x=1) %>%
        df_to_s(dims=c(N,length(q[q])))
    }
    if (min(rowSums(En))==0) {
      q <- rowSums(En) == 0
      En[q,] <- tibble(i=1:length(q[q]),
                       j=1+which(q),
                       x=1) %>%
        df_to_s(dims=c(length(q[q]),N))
    }
    En <- En %>% summary() %>% tbl_df() %>% filter(x>0) %>% df_to_s(dims=c(N,N))
  }
    
  G <- En
  G@x <- rlnorm(length(G@x),0,1) 
#  G@x <- runif(length(G@x))
  G <- G / rowSums(G)
  
  tol <- 1e-5
  obj <- tol + 1
  v1 <- v0 <- rep_len(1,N)
  while (obj > tol) {
    v0 <- v1
    I <- ir %*% (beta * v0)
    v1 <- t(A) %*% I + t(G) %*% (((1-beta) * v0) %>% matrix(nrow=N,ncol=1))
    v1 <- v1 / sum(v1)
    obj <- sqrt(sum((v1-v0)^2))
    # print(obj)
  }
  s <- v1
  I <- ir %*% (beta * s)
  
  # Scale
  I <- I/sum(s)*scale
  s <- s/sum(s)*scale
  
  # Industry IO table; should really include region final demand?
  x <- ((((1-beta)*s) %>% to_sdiag()) %*% G)  %>% summary() %>% tbl_df() 
  
  # merge on industries
  k <- ik %>% summary() %>% tbl_df() %>% select(k=i,i=j,-x)
  y <- x %>% left_join(k,by=c("i"="i")) %>% rename(ki=k) %>% left_join(k,by=c("j"="i")) %>% rename(kj=k)
  kk <- y %>% group_by(ki,kj) %>% summarize(x=sum(x)) %>% rename(i=ki,j=kj) %>% df_to_s(dims=c(K,K))
  dim(kk) <- c(K*K,1)
  
  return(list(A=A,beta=beta,G=G,I=I,ir=ir,ik=ik,kk=kk,s=s))
}

initialize_fake_s <- function(R,K,N,region_density,firm_density) {

  # Plant value-added share
  beta <- runif(N,0.1,0.9)
  
  # Plant output
  s <- rlnorm(N,0,3)

  # Plant's region
  ir <- rbind(tibble(j=1:N,i=sample(1:R,N,replace=TRUE),x=1)) %>% df_to_s(dims=c(R,N))
  
  # Plant's industry, total of J.
  # iz <- tibble(j=1:N,i=sample(1:Z,N,replace=TRUE),x=1) %>% df_to_s(dims=c(Z,N))
  ik <- tibble(j=1:N,i=rep(1:K,ceiling(N/K))[1:N],x=1) %>% df_to_s(dims=c(K,N))
  
    # Regional income
  I <- ir %*% (beta * s)
  
  Er <- rsparsematrix(R,N,density=region_density,rand.x=function(n) 1)

  # Non-zero edges of plant-plant demand matrix
  # Possibly add services demand to each one that doesn't have a thing.
  En <- rsparsematrix(N,N,density=firm_density,rand.x=function(n) 1)
  diag(En[1:N,1:N]) <- 0
  
  # Return fake data; don't need non-zero edge matrices anymore.
  return(list(beta=beta,Er=Er,En=En,I=I,ir=ir,ik=ik,s=s))
}
