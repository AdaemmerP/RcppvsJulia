"""
* Experiment: Draw 10^7 times from Bin(n = 100, p = 0.1) and compute coverage rates 
  for Wald, Wilson and Clopper-Pearson confidence intervals (95 %)		
* Codes 1 - 4 also use multiple threads but are not vectorized (slow versions)
* Code 5 computes the confidence intervals for Clopper-Pearson before the main loop. 
  In addition, it uses StaticArrays. (faster version)
* Code 5 also contains the function rbinom (line 14) from R which is faster than drawing the random numbers
  with the 'Distributions.jl' package
* Code 6 is fully vectorized (fastest version)
* The script 'PackagesF.jl' loads packages and helping functions
"""

# Load packages and functions
  include("PackagesF.jl")

# Include scripts to benchmark
  include("Code1.jl")
  include("Code2.jl")
  include("Code3.jl")
  include("Code4.jl")
  include("Code5.jl")
  include("Code6.jl")

# Set paramters
  p  = 0.1;
  n  = 100;     
  N  = 10^7; 
  α	 = 0.05;
	

# Run functions once for compilation	
  f1MC(n, p, 100, α) ;	f1SC(n, p, 100, α) 
  f2MC(n, p, 100, α) ;	f2SC(n, p, 100, α) 
  f3MC(n, p, 100, α) ;  f3SC(n, p, 100, α)   
  f4MC(n, p, 100, α) ;  f4SC(n, p, 100, α) 
  f5(n, p, 100, α, false) 	
  f6(n, p, 100, α, false) 


# Benchmarks
  # Part I
  @benchmark f1SC($n, $p, $N, $α)  seconds = 400 samples = 10
  @benchmark f1MC($n, $p, $N, $α)  seconds = 80  samples = 10
  @benchmark f2SC($n, $p, $N, $α)  seconds = 400 samples = 10
  @benchmark f2MC($n, $p, $N, $α)  seconds = 70  samples = 10

  # Part II
  @benchmark f3SC($n, $p, $N, $α)  seconds = 400 samples = 10
  @benchmark f3MC($n, $p, $N, $α)  seconds = 70  samples = 10
  @benchmark f4SC($n, $p, $N, $α)  seconds = 400 samples = 10
  @benchmark f4MC($n, $p, $N, $α)  seconds = 70  samples = 10

  # Part III
  rb = true; # Use rbinom() function from R? (Code5 & Code6)
  @benchmark f5($n, $p, $N, $α, $rb) seconds = 10 samples = 10
  @benchmark f6($n, $p, $N, $α, $rb) seconds = 10 samples = 10

#-------------------------------------------------------------------#
#--- Compare computation of Wald and Clopper-Pearson interval ------#	
#-------------------------------------------------------------------#

# Wald interval	
  function c1(x::Int64, z::Float64, n::Int64)

    # Wilson CI
      c0    = x/n;
      b0    = (1/(1 +  z^2/n))*(c0 + (z^2)/(2*n)) 
      b1    = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  )
    
    # Return CIs	
      [ (b0 - b1) (b0 + b1) ]

	end

# Clopper-Pearson interval
  function c2(x::Int64, n::Int64, a1::Float64, a2::Float64)

		
  # Clopper-Pearson CI
    cplow = (x == 0 ? 0.0 : quantile(Beta(x,     n - x + 1), a1))    
    cpup  = (x == n ? 1.0 : quantile(Beta(x + 1, n - x),     a2))  

  # Return CIs	
    [cplow cpup]

  end

# Benchmark and compare
  b1 = @benchmark c1(10, 1.65, 100)         samples = 100
  b2 = @benchmark c2(10, 100, α/2, 1 - α/2) samples = 100
  median(b2).time/median(b1).time 


  