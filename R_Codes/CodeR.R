
  library(microbenchmark)
  library(parallel)
  library(rlist)
  
#-------------------------------- Loop version --------------------------------# 
  
  # Function for all CIs	
  ciallf = function(x, z, n, a1, a2){
    
    # Wald Intervall	
    c0 = x/n
    c1 = z*sqrt(c0*(1-c0)/n)
    
    # Wilson Intervall	
    b0 = (1/(1+ z^2/n))*(c0 + (z^2)/(2*n))
    b1 = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  )
    
    cplow = qbeta(a1, x, n - x + 1)
    cpup =  qbeta(a2, x + 1, n - x)
    
    return (c( (c0 - c1), (c0 + c1), 
               (b0 - b1), (b0 + b1), 
               cplow, cpup))
    
    
  }
  
  # Set parameters
    p          = 0.1
    N          = 1e7
    n          = 100
    a_level    = 0.05
    
    
    # Benchmarking
    microbenchmark(
    {
    z          = qnorm(1 - a_level/2)   
    ciall      = matrix(NA, N, 6)
    a1         = a_level/2 
    a2         = (1 - a_level/2)
   
    # Draw from binomial distribution
    X = rbinom(N, n, p)
    
    # Get confidence intervals  
    for(jj in 1:N) ciall[jj, ] = ciallf(X[jj], z, n, a1, a2)
  
    sum(ciall[, 1] <= p & p <=  ciall[, 2])/N 
    sum(ciall[, 3] <= p & p <=  ciall[, 4])/N 
    sum(ciall[, 5] <= p & p <=  ciall[, 6])/N 
    }, 
    
    {
    resconf <- mclapply(X, function(x) ciallf(x, z, n, a1, a2),
                        mc.cores = 12)
    resconf <- list.rbind(resconf)
    
    sum(resconf[, 1] <= p & p <=  resconf[, 2])/N 
    sum(resconf[, 3] <= p & p <=  resconf[, 4])/N 
    sum(resconf[, 5] <= p & p <=  resconf[, 6])/N }, 
    
    
    times = 1)
    
    
 #------------------------- Clopper-Pearson vectorized version ----------------#
  
  
# Function for all CIs	
  ciallf = function(x, z, n, cplow, cpup){
    
    # Wald Intervall	
    c0 = x/n
    c1 = z*sqrt(c0*(1-c0)/n)
    
    # Wilson Intervall	
    b0 = (1/(1+ z^2/n))*(c0 + (z^2)/(2*n))
    b1 = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  )

    return (c( (c0 - c1)    <= p & p <= (c0 + c1), 
               (b0 - b1)    <= p & p <= (b0 + b1), 
               cplow[x + 1] <= p & p <= cpup[x + 1]))
  
  }
  
  # Set parameters
    p          = 0.1
    N          = 1e7
    n          = 100
    a_level    = 0.05

# Benchmarking
  microbenchmark({
 
    z          = qnorm(1 - a_level/2)   
    ciall      = matrix(NA, N, 3)
    a1         = a_level/2 
    a2         = (1 - a_level/2)
    cplow      = sapply(0:n, function(x)  qbeta(a1, x, n - x + 1))
    cpup       = sapply(0:n, function(x)  qbeta(a2, x + 1, n - x) )

  # Draw from binomial distribution
    X = rbinom(N, n, p)
    
  # Get confidence intervals  
    for(jj in 1:N) ciall[jj, ] = ciallf(X[jj], z, n, cplow, cpup)
    
  # Compute coverage  
    sum(ciall[, 1])/N
    sum(ciall[, 2])/N
    sum(ciall[, 3])/N
                
        }, times = 1)
  
  
#------------------------- Fully Vectorized version -------------------------#  
  
  # Set parameters
  p          = 0.1
  N          = 1e7
  n          = 100
  a_level    = 0.05
  
  # Function for all CIs	
  ciallf = function(x, z, n, a1, a2){
    
    # Wald Intervall	
    c0 = x/n
    c1 = z*sqrt(c0*(1-c0)/n)
    
    # Wilson Intervall	
    b0 = (1/(1+ z^2/n))*(c0 + (z^2)/(2*n))
    b1 = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  )
    
    cplow = qbeta(a1, x, n - x + 1)
    cpup  = qbeta(a2, x + 1, n - x)
  
    return (c( (c0 - c1) <= p & p <=  (c0 + c1), 
               (b0 - b1) <= p & p <=  (b0 + b1),
               cplow     <= p & p <= cpup) )
    
  }
  
  # Benchmarking
  microbenchmark(
    
    {
      
    z          = qnorm(1 - a_level/2)   
    a1         = a_level/2 
    a2         = (1 - a_level/2)
    covall     = sapply(0:n, function(x) ciallf(x, z, n, a1, a2))
    cowi       = covall[1, ]
    cowa       = covall[2, ]
    cocp       = covall[3, ]
  
    # Draw from binomial distribution
      X = rbinom(N, n, p)

      cowi  = cowi[X + 1]
      cowa  = cowa[X + 1]
      cocp  = cocp[X + 1] 
      
    # Compute coverage 
      sum(cowi)/N
      sum(cowa)/N
      sum(cocp)/N
      
     }, 
    
 times = 10)
  

# Use binom package
  library(binom)
  n <- 100
  p <- 0.1
  N <- 1e7
  microbenchmark(
      {X <- rbinom(N, n, p)   
      cp = binom.confint(X, n, methods = "exact")[, 5:6] 
      wi = binom.confint(X, n, methods = "wilson")[, 5:6]
      wa = binom.confint(X, n, methods = "asymptotic")[, 5:6]
      sum(cp[, 1] <= p & p <= cp[, 2])/N
      sum(wi[, 1] <= p & p <= wi[, 2])/N
      sum(wa[, 1] <= p & p <= wa[, 2])/N}, times = 1)

# Compute asymptotic coverage rates
  covall     = sapply(0:n, function(x) ciallf(x, z, n, a1, a2))
  cowi       = covall[1, ]
  cowa       = covall[2, ]
  cocp       = covall[3, ]
  sum(dbinom(which(cowi == T) - 1, 100, 0.1))
  sum(dbinom(which(cowa == T) - 1, 100, 0.1))
  sum(dbinom(which(cocp == T) - 1, 100, 0.1))  
  
  
  
  