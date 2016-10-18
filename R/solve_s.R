########################################################################
# solve_s.R
# License: MIT

# ""
# Jesse Tweedle
# Oct, 2016
########################################################################

solve_s <- function(R,N,args) {

  ir <- args$ir
  beta <- args$beta
  A <- args$A
  G <- args$G
  
# maybe this should have it's own initial level. yes, def.
  if (is.null(s)) {
    s0 <- s1 <- rep_len(1,N)
  } else {
    s0 <- s1 <- s 
  }

  tol <- 1e-5 / N
  obj <- tol + 1

  while (obj > tol) {
    s0 <- s1
    I <- ir %*% ((beta*s0) %>% matrix(nrow=N,ncol=1))
    s1 <- t(A) %*% I + t(G) %*% s0
    s1 <- s1 / sum(s1) # normalize weights.
    obj <- sum(abs(s1-s0))
  }

  return(s=s1)
}
