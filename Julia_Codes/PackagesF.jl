
# Load packages
  using Distributions
  using HypothesisTests
  using ThreadsX
  using BenchmarkTools
  using StaticArrays
  using RCall
		
# Manual Functions 

# Function for 'Code3.jl'
  function ciallf(x::Int64, z::Float64, n::Int64, a1::Float64, a2::Float64)

      # Wald CI
        c0  = x/n;
        c1  = z*sqrt(c0*(1 - c0)/n);

      # Wilson CI
        b0    = (1/(1 +  z^2/n))*(c0 + (z^2)/(2*n))
        b1    = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  )

      # Clopper-Pearson CI
        cplow = (x == 0 ? 0.0 : quantile(Beta(x,     n - x + 1), a1))    
        cpup  = (x == n ? 1.0 : quantile(Beta(x + 1, n - x),     a2))  

      # Return CIs	
        [(c0 - c1) (c0 + c1) (b0 - b1) (b0 + b1)  cplow cpup]

  end

# Function for 'Code5.jl'
  function ciallf(x::Int64, z::Float64, n::Int64, 
                  cplow::Vector{Float64}, cpup::Vector{Float64}, p::Float64)

      # Wald
        c0 = x/n;
        c1 = z*sqrt(c0*(1 - c0)/n);

      # Wilso
        b0  = (1/(1+ z^2/n))*(c0 + (z^2)/(2*n))
        b1  = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  )

      # Return CIs
        SVector{3, Bool}((c0 - c1)     <= p <= (c0 + c1),
                         (b0 - b1)     <= p <= (b0 + b1), 
                          cplow[x + 1] <= p <= cpup[x + 1])   

  end

# Function for 'Code6.jl' (fully vectorized version)
  function ciallfv(x::Int64, z::Float64, n::Int64, 
                   p::Float64, a1::Float64, a2::Float64)
      
          # Wald
            c0  = x/n;
            c1  = z*sqrt(c0*(1 - c0)/n);

          # Wilson
            b0  = (1/(1+ z^2/n))*(c0 + (z^2)/(2*n))
            b1  = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  )

          # Clopper Pearson
            cplow =  (x == 0 ? 0.0 : quantile(Beta(x, n - x + 1), a1)) 
            cpup  =  (x == n ? 1.0 : quantile(Beta(x + 1, n - x), a2)) 

          # Return CIs
            SMatrix{3, 1, Bool}((c0 - c1)  <= p <= (c0 + c1),  
                                (b0 - b1)  <= p <= (b0 + b1), 
                                cplow      <= p <= cpup)   

  end



