########################################################################
# initialize_functions.R
# Functions to generate fake parameters and data
# License: MIT

# ""
# Jesse Tweedle
# Oct, 2016
########################################################################

initialize_fake_s <- function(R,Z,N) {

  # Plant value-added share
  beta <- runif(N,0.1,0.9)
  
  # Plant output
  s <- rlnorm(N,0,1)
  
  # Plant's region
  ir <- rbind(tibble(j=1:N,i=sample(1:R,N,replace=TRUE),x=1)) %>% df_to_s(dims=c(R,N))
  
  # Plant's industry, total of J.
  # iz <- tibble(j=1:N,i=sample(1:Z,N,replace=TRUE),x=1) %>% df_to_s(dims=c(Z,N))
  iz <- tibble(j=1:N,i=rep(1:Z,ceiling(N/Z))[1:N],x=1) %>% df_to_s(dims=c(Z,N))
  
    # Regional income
  I <- ir %*% (beta * s)
  
  # Return fake data; don't need non-zero edge matrices anymore.
  return(list(beta=beta,I=I,ir=ir,iz=iz,s=s))
}
