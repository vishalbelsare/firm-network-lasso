########################################################################
# main.R
# Clean and load relevant files
# License: MIT

# ""
# Jesse Tweedle
# Oct 18, 2016
########################################################################

library(dplyr)
library(tibble)
library(Matrix)
library(stringr)
library(tidyr)
library(ggplot2)
library(glmnet)

rm(list=ls(all=TRUE))
gc()

# Tip from http://jeromyanglim.tumblr.com/post/33418725712/how-to-source-all-r-files-in-a-directory
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir(paste0(getwd(),"/R"))

set.seed(9) # set.seed(10), R=2, N=300 screws up.
R <- 3
N <- 5

print("fake edges")
argsx <- initialize_fake_links(R,N)

print("fake shares")
args <- initialize_fakes(R,N,args=argsx)

Z <- c(args$I[,1],args$s) %>% matrix(nrow=1,ncol=R+N)
#Z %>% matrix(nrow=1,ncol=length(Z)) %>% bdiag()
# try tibble? rep_len?
#tibble(i=1:N,j=((i-1)*(R+N)+1):(i*(R+N)),x=1)
#expand.grid(i=1:N,j=((i-1)*(R+N)+1):(i*(R+N)))
# y = ?

#bdiag(list(Z,Z))
X <- Z %>% list() %>% rep(N) %>% bdiag()
# c = s
c <- args$s

# now have X, c.

fit <- glmnet(X,c,intercept=FALSE)

# s is lambda, higher s means more sparse.
y <- coef(fit,s=0.1)
#coef(fit)
# now turn this matrix back into A and G

# y <- a,g,a,g,etc.
y <- y[-1,1]
a1 <- y[1:R]
g1 <- y[(R+1):(R+N)]


# split into 1:(R*N+N^2)



#