
# Full vectorized version
  function f6(n::Int64, p::Float64, N::Int64, α::Float64, rb::Bool) 
	
  # Set parameters 
	  z  = quantile(Normal(), 1 - α/2) 
    a1 = α/2
    a2 = (1 - α/2)
	
  # Used Distributions package or rbinom?	
    if rb == 0	

    # Use Julia to draw from Bin(n, p)
      dist    = Binomial(n, p)
      randnrs = rand(dist, N)

    else	

    # Use R to draw from Bin(n, p)	
      randnrs::Vector{Int64} = (rcopy(Vector{Int64}, R"rbinom($N, $n, $p)"))

    end

  # Get all possible CIs	
    cico::Vector{SMatrix{3, 1, Bool, 3}} = map(x ->  ciallfv(x, z, n,	p, a1, a2), 0:1:n)

  # Compute coverage rates 		
  [ sum(getindex.(cico, 1)[randnrs .+ 1])/N
    sum(getindex.(cico, 2)[randnrs .+ 1])/N
    sum(getindex.(cico, 3)[randnrs .+ 1])/N ] 
		
 end 