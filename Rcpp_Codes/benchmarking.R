# load packages
library(microbenchmark)
library(Rcpp)

# Source scripts
sourceCpp("singlecore.cpp")
sourceCpp("multicore.cpp")
sourceCpp("halfvectorized.cpp")
sourceCpp("fullvectorized.cpp")

# set parameters
n <- 100
p <- .1
N <- 1e7
alpha <- .05
RcppParallel::setThreadOptions(numThreads = 12)

# test functions
CIrcpp1(n, p, N, alpha)
CIrcpp2(n, p, N, alpha)
CIrcpp3(n, p, N, alpha)
CIrcpp4(n, p, N, alpha)

# benchmark
microbenchmark(
  CIrcpp1(n, p, N, alpha),
  CIrcpp2(n, p, N, alpha),
  CIrcpp3(n, p, N, alpha),
  CIrcpp4(n, p, N, alpha),
  times=10)
