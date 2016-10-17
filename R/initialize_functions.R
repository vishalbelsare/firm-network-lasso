########################################################################
# initialize_functions.R
# Functions to generate fake parameters and data
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

initialize_fake_links <- function(R,N) {

  # Plant value-added share
  beta <- runif(N,0.1,0.9)

  # Plant output
  s <- rlnorm(N,0,1) * N

  # Plant's region
  ir <- rbind(tibble(j=1:(N),i=sample(1:R,N,replace=TRUE),x=1)) %>% df_to_s(dims=c(R,N))

  # Regional income
  I <- ir %*% (beta * s)

  # Return generated parameters.
  return(list(beta=beta,I=I,ir=ir,s=s))
}

initialize_fakes <- function(R,N,args) {

  if (missing(args)) {
    args <- initialize_fake_links(R,N)
  }

  # Region income
  I <- args$I

  # Plant's region
  ir <- args$ir

  # Value-added share
  beta <- args$beta

  # Output
  s <- args$s

  # Return fake data; don't need non-zero edge matrices anymore.
  return(list(beta=beta,I=I,ir=ir,s=s))
}
