########################################################################
# initialize_functions.R
# Functions to generate fake parameters and data
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

initialize_fake_s <- function(R,Z,N) {

  # Plant value-added share
  beta <- runif(N,0.1,0.9)
  
  # Plant output
  s <- rlnorm(N,0,1) * N
  
  # Plant's region
  ir <- rbind(tibble(j=1:N,i=sample(1:R,N,replace=TRUE),x=1)) %>% df_to_s(dims=c(R,N))
  
  # Plant's industry, total of J.
  iz <- rbind(tibble(j=1:N,i=sample(1:Z,N,replace=TRUE),x=1)) %>% df_to_s(dims=c(Z,N))

    # Regional income
  I <- ir %*% (beta * s)
  
  # Return fake data; don't need non-zero edge matrices anymore.
  return(list(beta=beta,I=I,ir=ir,iz=iz,s=s))
}
