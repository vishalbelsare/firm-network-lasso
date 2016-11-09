########################################################################
# helpers.R
# Functions to convert to and from sparse matrices
# License: MIT

# "Granularity, Network Asymmetry and Aggregate Volatility"
# Jesse Tweedle
# Sep 15, 2016
########################################################################

load_libraries <- function(inside_parallel=FALSE) {
  # Load relevant libraries:
  
  # For sparse matrices
  library(Matrix)
  
  # Tidy packages
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(tidyr)
  # library(tidyverse) # don't have tidyverse at work.
    
  # For Lasso.
  library(glmnet)
  
#  library(xgboost)
  
  if (!inside_parallel) {
    # Tidy packages
    library(tibble)
    library(stringr)
    library(tidyr)
    library(ggplot2)
    
    # for parallelization.
    library(parallel)
  }
}
startup <- function() {

  # Stop on warnings. Much easier for debugging.
#  options(warn = 2)
  
  sourceDir(paste0(getwd(),"/R"))
  
  load_libraries()
}

# Source all .R files to get functions. Make sure
# there are only functions in the R/ directory.
sourceDir <- function(path, trace = TRUE, ...) {
  # Tip from http://jeromyanglim.tumblr.com/post/33418725712/how-to-source-all-r-files-in-a-directory
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# Convert sparseMatrix to tbl_df.
s_to_df <- function(m) {
  m %>% as("sparseMatrix") %>% summary() %>% tbl_df()
}

# Convert data.frame to a sparseMatrix.
df_to_s <- function(tdf,dims) {
  with(tdf, sparseMatrix(i=i,j=j,x=x,dims=dims))
#  sparseMatrix(i=tdf$i,j=tdf$j,x=tdf$x,dims=dims) # make sure dims here; if last i or j is missing, dims might not match.
}

# Convert either a one-dimensional matrix or a vector to a sparse diagonal matrix.
# Useful alternative to perform elementwise multiplication on large sparse matrices.
to_sdiag <- function(x) {
  if (class(x)=="dgCMatrix" | class(x)=="dgeMatrix") {
    # If x is a matrix:
    return(.sparseDiagonal(n=length(x[,1]),x=x[,1]))
  } else if (is.vector(x)) {
    # If x is a vector:
    return(.sparseDiagonal(n=length(x),x=x))
  }
}
