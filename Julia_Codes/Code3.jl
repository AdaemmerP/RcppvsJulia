
# Parallelized
  function f3MC(n::Int64, p::Float64, N::Int64, α::Float64) 
	
  # Compute quantile
    z  			= quantile(Normal(), 1 - α/2)
	
  # Draw # of successes from Bin(n, p)
    dist    = Binomial(n, p)
    randnrs = rand(dist, N)
	
  # Compute CIs 
    ciall::Vector{Matrix{Float64}} = ThreadsX.map(x -> ciallf(x, z, n,  α/2, (1 - α/2)), randnrs)
			 
  # Compute coverage rates 
    [ sum(getindex.(ciall, 1) .<= p .<= getindex.(ciall, 2))/N
      sum(getindex.(ciall, 3) .<= p .<= getindex.(ciall, 4))/N
      sum(getindex.(ciall, 5) .<= p .<= getindex.(ciall, 6))/N  ]

  end

# Single core
  function f3SC(n::Int64, p::Float64, N::Int64, α::Float64) 
	
  # Compute quantile
    z  			= quantile(Normal(), 1 - α/2)
	
  # Draw # of successes from Bin(n, p)
    dist    = Binomial(n, p)
    randnrs = rand(dist, N)
	
  # Compute CIs 
    ciall::Vector{Matrix{Float64}} = map(x -> ciallf(x, z, n,  α/2, (1 - α/2)), randnrs)
			 
	# Compute coverage rates 
    [ sum(getindex.(ciall, 1) .<= p .<= getindex.(ciall, 2))/N
      sum(getindex.(ciall, 3) .<= p .<= getindex.(ciall, 4))/N
      sum(getindex.(ciall, 5) .<= p .<= getindex.(ciall, 6))/N  ]

  end	
  

