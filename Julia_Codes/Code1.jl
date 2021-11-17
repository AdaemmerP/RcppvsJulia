# Parallelized
  function f1MC(n, p, N, α) 

  # Pre-allocate		
    waldci    = zeros(N, 2)
    wilsonci 	= zeros(N, 2)
    cpci      = zeros(N, 2)

    # Loop 
    Threads.@threads for jj = 1:N

      # Compute boolean vector
        randnrs = rand(n) .<= p 

      # Compute CIs
        waldtemp         = confint(BinomialTest(randnrs, p), level = (1 - α), method = :wald)
        wilsontemp       = confint(BinomialTest(randnrs, p), level = (1 - α), method = :wilson)
        cptemp           = confint(BinomialTest(randnrs, p), level = (1 - α), method = :clopper_pearson) 

      # Fill matrices with CIs				  
        waldci[jj, 1]    = waldtemp[1]
        waldci[jj, 2]    = waldtemp[2]

        wilsonci[jj, 1]	 = wilsontemp[1]
        wilsonci[jj, 2]  = wilsontemp[2]

        cpci[jj, 1]      = cptemp[1]
        cpci[jj, 2]      = cptemp[2]

    end
	
# Compute coverage rates 
  [ sum(waldci[:, 1]   .<= p .<= waldci[:, 2])/N 
    sum(wilsonci[:, 1] .<= p .<= wilsonci[:, 2])/N 
    sum(cpci[:, 1]     .<= p .<= cpci[:, 2])/N 			]

end

#-----------------------------------------------------------------------------#

# Single core 
  function f1SC(n, p, N, α) 

    # Pre-allocate		
      waldci    = zeros(N, 2)
      wilsonci 	= zeros(N, 2)
      cpci      = zeros(N, 2)

      # Loop 
      for jj = 1:N

       # Compute boolean vector
         randnrs = rand(n) .<= p 

       # Compute CIs
         waldtemp         = confint(BinomialTest(randnrs, p), level = (1 - α), method = :wald)
         wilsontemp       = confint(BinomialTest(randnrs, p), level = (1 - α), method = :wilson)
         cptemp           = confint(BinomialTest(randnrs, p), level = (1 - α), method = :clopper_pearson) 

       # Fill matrices with CIs				  
         waldci[jj, 1]  = waldtemp[1]
         waldci[jj, 2]  = waldtemp[2]

         wilsonci[jj, 1]  = wilsontemp[1]
         wilsonci[jj, 2]  = wilsontemp[2]

         cpci[jj, 1]      = cptemp[1]
         cpci[jj, 2]      = cptemp[2]

     end
		
  # Compute coverage rates 
    [ sum(waldci[:, 1]   .<= p .<= waldci[:, 2])/N 
      sum(wilsonci[:, 1] .<= p .<= wilsonci[:, 2])/N 
      sum(cpci[:, 1]     .<= p .<= cpci[:, 2])/N      ]

  end
