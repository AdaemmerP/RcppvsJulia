# Function to benchmark
  function f5(n::Int64, p::Float64, N::Int64, α::Float64, rb::Bool)

  # Set parameters
    z  = quantile(Normal(), 1 - α/2)
    a1 = α/2
    a2 = (1 - α/2)

  # Use rbinom from R?  
    if rb == 0

    # Use Julia to draw from Bin(n, p)
      dist    = Binomial(n, p)
      randnrs = rand(dist, N)

    else

    # Use R to draw from Bin(n, p)
      randnrs::Vector{Int64} = rcopy(Vector{Int64}, R"rbinom($N, $n, $p)") #

    end

  # Compute all possible lower and upper bound for Clopper-Pearson interval 
    cplow =  map(x -> x == 0 ? 0.0 : quantile(Beta(x, n - x + 1), a1), 0:1:n)
    cpup  =  map(x -> x == n ? 1.0 : quantile(Beta(x + 1, n - x), a2), 0:1:n)

  # Get intervals
    ciall =  map(x -> ciallf(x, z, n, cplow, cpup, p),  randnrs)

  # Compute coverage rates
   [ mapreduce(x -> x[1], + , ciall)/N
     mapreduce(x -> x[2], + , ciall)/N
     mapreduce(x -> x[3], + , ciall)/N   ]

  end



