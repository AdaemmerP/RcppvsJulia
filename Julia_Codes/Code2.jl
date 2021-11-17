# Parallelized
  function f2MC(n, p, N, α) 

  # Compute boolean Matrix
    randnrs   = rand(n, N) .<= p 

  # Compute CIs
    waldci    = ThreadsX.map(x -> confint(BinomialTest(x, p), level = (1 - α),  method = :wald),    eachcol(randnrs))
    wilsonci  = ThreadsX.map(x -> confint(BinomialTest(x, p), level = (1 - α),  method = :wilson),  eachcol(randnrs))
    cpci      = ThreadsX.map(x -> confint(BinomialTest(x, p), level = (1 - α),  method = :clopper_pearson), eachcol(randnrs))

  # Compute coverage rates 
    [ sum(first.(waldci)   .<=  p .<=  last.(waldci))/N
      sum(first.(wilsonci) .<=  p .<=  last.(wilsonci))/N
      sum(first.(cpci)     .<=  p .<=  last.(cpci))/N      ]

  end

#-----------------------------------------------------------------------------#

# Single core 	
  function f2SC(n, p, N, α) 

    # Compute boolean Matrix
      randnrs   = rand(n, N) .<= p 
	
    # Compute CIs
      waldci    = map(x -> confint(BinomialTest(x, p), level = (1 - α),  method = :wald),   eachcol(randnrs))
      wilsonci  = map(x -> confint(BinomialTest(x, p), level = (1 - α),  method = :wilson), eachcol(randnrs))
      cpci      = map(x -> confint(BinomialTest(x, p), level = (1 - α),  method = :clopper_pearson), eachcol(randnrs))
	
    # Compute coverage rates 
      [ sum(first.(waldci)   .<=  p .<=  last.(waldci))/N
        sum(first.(wilsonci) .<=  p .<=  last.(wilsonci))/N
        sum(first.(cpci)     .<=  p .<=  last.(cpci))/N     ]
	
  end

