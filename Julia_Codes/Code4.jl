# Parallelized
  function f4MC(n::Int64, p::Float64, N::Int64, α::Float64) 
		
    # Draw # of successes 
      dist    = Binomial(n, p)
      randnrs = rand(dist, N)
				
    # Compute CIs 
      waldci::Vector{Tuple{Float64, Float64}}    = ThreadsX.map(x -> confint(BinomialTest(x, n, p), level = (1 - α),  
                                                                             method = :wald),   randnrs)	
      wilsonci::Vector{Tuple{Float64, Float64}}  = ThreadsX.map(x -> confint(BinomialTest(x, n, p), level = (1 - α),  
                                                                             method = :wilson), randnrs)
      cpci::Vector{Tuple{Float64, Float64}}      = ThreadsX.map(x -> confint(BinomialTest(x, n, p), level = (1 - α),  
                                                                             method = :clopper_pearson), randnrs)

    # Compute coverage rates
      [ sum(first.(waldci)   .<= p .<=  last.(waldci))/N
        sum(first.(wilsonci) .<= p .<=  last.(wilsonci))/N
        sum(first.(cpci)     .<= p .<=  last.(cpci))/N  ]
					
   end

# Single core
  function f4SC(n::Int64, p::Float64, N::Int64, α::Float64) 
			
    # Draw # of successes 
      dist    = Binomial(n, p)
      randnrs = rand(dist, N)
				
    # Compute CIs 
      waldci::Vector{Tuple{Float64, Float64}}    = map(x -> confint(BinomialTest(x, n, p), level = (1 - α),  
                                                                    method = :wald),            randnrs)	
      wilsonci::Vector{Tuple{Float64, Float64}}  = map(x -> confint(BinomialTest(x, n, p), level = (1 - α),  
                                                                    method = :wilson),          randnrs)
      cpci::Vector{Tuple{Float64, Float64}}      = map(x -> confint(BinomialTest(x, n, p), level = (1 - α),  
                                                                    method = :clopper_pearson), randnrs)

    # Compute coverage rates
      [ sum(first.(waldci)   .<= p .<=  last.(waldci))/N
        sum(first.(wilsonci) .<= p .<=  last.(wilsonci))/N
        sum(first.(cpci)     .<= p .<=  last.(cpci))/N     ]
					
  end



