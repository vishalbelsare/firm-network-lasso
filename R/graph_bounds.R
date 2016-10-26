########################################################################
# graph_bounds.R
# Clean and load relevant files
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
########################################################################

######## SEPARATE CODE THAT DOES UPPER/LOWER BOUND


upper_bound <- function(A,G,region_density=0.5,firm_density=0.25) {
  
  Er <- rsparsematrix(R,N,density=min(region_density*2,1),rand.x=function(n) 1)
  Er <- Er + A
  Er@x <- rep_len(1,length(Er@x))

  En <- rsparsematrix(N,N,density=min(firm_density*2,1),rand.x=function(n) 1)
  En <- En + G
  En@x <- rep_len(1,length(En@x))

  return(rbind(Er,En))
}

lower_bound <- function(A,G,region_sample_frac=1,firm_sample_frac=1) {
  
  Er <- A %>% summary() %>% tbl_df() %>% mutate(x=1) %>% sample_frac(size=region_sample_frac) %>% df_to_s(dims=c(R,N))
  En <- G %>% summary() %>% tbl_df() %>% mutate(x=1) %>% sample_frac(size=firm_sample_frac) %>% df_to_s(dims=c(N,N))
  
  return(rbind(Er,En))
}