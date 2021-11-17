library(JuliaCall)
library(microbenchmark)
julia <- julia_setup()
julia_library("Distributions")
julia_library("StaticArrays")


#-----------------------  Vectorize CP intervals only  ----------------------#


# Function for to compute and extract CIs
  julia_command("
       function	ciallf(x::Int64, z::Float64, n::Int64, 
  						cplow::Vector{Float64}, cpup::Vector{Float64}, p::Float64)
  
  					# Wald
  						c0 = x/n;
  						c1 = z*sqrt(c0*(1 - c0)/n);
  
  					# Wilson
  						b0  = (1/(1+ z^2/n))*(c0 + (z^2)/(2*n))
  						b1  = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  )
  
  					#Return CIs
  					SVector{3, Bool}((c0 - c1)   <= p <= (c0 + c1),  
  					 								(b0 - b1)  	 <= p <= (b0 + b1), 
  													cplow[x + 1] <= p <= cpup[x + 1])   
  
  		end")

# Function which uses external vector of random numbers
  julia_command("
      function fR(n::Int64, p::Float64, N::Int64, α::Float64, randnrs::Vector{Int64})  
     
      # Set parameters 
        p  = 0.1;  n  = 100;	   
        N  = 10^7; α  = 0.05;      
        z  = quantile(Normal(), 1 - α/2); a1 = α/2; 
        a2 = (1 - α/2)
    
      # Compute all possible lower and upper bound for Clopper-Pearson interval
        cplow =  map(x -> x == 0 ? 0.0 : quantile(Beta(x, n - x + 1), a1), 0:1:n)
        cpup  =  map(x -> x == n ? 1.0 : quantile(Beta(x + 1, n - x), a2), 0:1:n)
      
      # Get intervals
        ciall =  map(x -> ciallf(x, z, n, cplow, cpup, p),  randnrs)  
      
      # Compute coverage rates 
      [	mapreduce(x -> x[1], + , ciall)/N
        mapreduce(x -> x[2], + , ciall)/N
        mapreduce(x -> x[3], + , ciall)/N  ] 
  
  end") 

# Function which uses 'Distributions' to draw from Bin(n, p)
  julia_command("
      function fJ(n::Int64, p::Float64, N::Int64, α::Float64)  
     
      # Set parameters 
        z  = quantile(Normal(), 1 - α/2) 
        a1 = α/2; 
        a2 = (1 - α/2)
      
     # Use Julia to draw from Bin(n, p)
       dist    = Binomial(n, p)
       randnrs = rand(dist, N)
      
     # Compute all possible lower and upper bound for Clopper-Pearson interval
       cplow =  map(x -> x == 0 ? 0.0 : quantile(Beta(x, n - x + 1), a1), 0:1:n)
       cpup  =  map(x -> x == n ? 1.0 : quantile(Beta(x + 1, n - x), a2), 0:1:n)
      
     # Get intervals
       ciall =  map(x -> ciallf(x, z, n, cplow, cpup, p),  randnrs)  
      
     # Compute coverage rates 
      [	mapreduce(x -> x[1], + , ciall)/N
        mapreduce(x -> x[2], + , ciall)/N
        mapreduce(x -> x[3], + , ciall)/N  ] 
  
    end") 
  
  # Set paramters
    n     <- 100L
    p     <- 0.1
    N     <- 1e7L
    alpha <- 0.05

# Run functions once for compilations
  julia_eval("fJ")(n, p, N, alpha)
  x <- rbinom(10^7, 100, 0.1)
  julia_eval("fR")(n, p, N, alpha, x)

# Benchmark
  microbenchmark(
      julia_eval("fJ")(n, p, N, alpha), 
      {x <- rbinom(10^7, 100, 0.1)
      julia_eval("fR")(n, p, N, alpha, x)},
    times = 10)
  
  
#----------------------------- Fully vectorized -------------------------------#
           
# Compute all quantiles beforehand
  julia_command("
    function ciallfv(x::Int64, z::Float64, n::Int64, 
											 p::Float64, a1::Float64, a2::Float64)::SMatrix{3, 1, Bool} 

					# Wald
						c0 = x/n;
						c1 = z*sqrt(c0*(1 - c0)/n);

					# Wilson
						b0  = (1/(1+ z^2/n))*(c0 + (z^2)/(2*n))
						b1  = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  )

					# Clopper Pearson
						cplow =  (x == 0 ? 0.0 : quantile(Beta(x, n - x + 1), a1)) 
						cpup  =  (x == n ? 1.0 : quantile(Beta(x + 1, n - x), a2)) 

					# Return CIs
						SMatrix{3, 1, Bool}((c0 - c1)  <= p <= (c0 + c1),  
																(b0 - b1)  <= p <= (b0 + b1), 
																cplow 		 <= p <= cpup)   

  end") 
  
 # Use Distributions package  
  julia_command("
    function vecfR(n::Int64, p::Float64, N::Int64, α::Float64, randnrs::Vector{Int64}) 
	
      	# Set parameters 
      		z  = quantile(Normal(), 1 - α/2) 
      		a1 = α/2
      		a2 = (1 - α/2)
      
      	# Get all possible CIs	
      		cico = map(x ->  ciallfv(x, z, n,	p, a1, a2), 0:1:n)
      		cico = reduce(hcat, cico)'[randnrs .+ 1, :]
      
      	# Compute coverage rates 		
      		[	sum(cico[:, 1])/N   
      	 		sum(cico[:, 2])/N   
      	 	  sum(cico[:, 3])/N ] 
		
    end") 
  
 # Use external vector with binomial draws    
  julia_command("
    function vecf(n::Int64, p::Float64, N::Int64, α::Float64) 
	
      	# Set parameters 
      		z  = quantile(Normal(), 1 - α/2) 
      		a1 = α/2
      		a2 = (1 - α/2)
      	
      	# Use Julia to draw from Bin(n, p)
      		 dist    = Binomial(n, p)
      		 randnrs = rand(dist, N)
      
      	# Get all possible CIs	
      		cico = map(x ->  ciallfv(x, z, n,	p, a1, a2), 0:1:n)
      		cico = reduce(hcat, cico)'[randnrs .+ 1, :]
      
      	# Compute coverage rates 		
      		[	sum(cico[:, 1])/N   
      	 		sum(cico[:, 2])/N   
      	 	  sum(cico[:, 3])/N ] 
		
    end") 
  
  # Set paramters
  n     <- 100L
  p     <- 0.1
  N     <- 1e7L
  alpha <- 0.05
  
 microbenchmark(
   julia_eval("vecf")(n, p, N, alpha), 
   {x    <- rbinom(10^7, 100, 0.1)
   julia_eval("vecfR")(n, p, N, alpha, x)},
    times = 10)
 
 
