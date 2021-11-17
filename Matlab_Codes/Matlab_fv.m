
time_vec = zeros(10, 1)
% Rinit('stats')

for tt = 1:10
    
tic()	
% Set parameters and pre-allocate matrices
  p         = 0.1;	    
  n         = 100;	   
  N         = 10^7;
  alpha     = 0.05;      
  z         = norminv(1 - alpha/2);
  a1        = alpha/2; 
  a2        = (1 - alpha/2); 
  walow     = zeros(n, 1); 
  waup      = zeros(n, 1); 
  wilow     = zeros(n, 1); 
  wiup      = zeros(n, 1); 
  cplow     = zeros(n, 1); 
  cpup      = zeros(n, 1); 
  % Rpush('p', p, 'n', n, 'N', N)
    
% Main loop to compute CIs    
  for ii = 1:(n + 1) 
       
      x          = ii - 1;
      cplow(ii)  = betainv(a1, x,     n - x + 1); 
      cpup(ii)   = betainv(a2, x + 1, n - x); 
     
     % Wald
      c0            = x/n;
      c1            = z*sqrt(c0*(1 - c0)/n);
      walow(ii)     = c0 - c1 ; 
      waup(ii)      = c0 + c1 ; 

	 % Wilson
      b0            = (1/(1+ z^2/n))*(c0 + (z^2)/(2*n));
      b1            = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  );
      wilow(ii)     = b0 - b1;
      wiup(ii)      = b0 + b1;
          
     % CP
      cplow(ii)  = betainv(a1, x,     n - x + 1); 
      cpup(ii)   = betainv(a2, x + 1, n - x); 
        
  end
  
  % Replace NaN 
    cplow(1)  = 0;
    cpup(end) = 1;
   
  % Draw N times from Bin(n, p) with RCall
  %  Rrun('randnrs <- rbinom(N, n, p)')
  %  randnrs = Rpull('randnrs');
    
  % Draw N times from Bin(n, p) with Matlab function
    randnrs = binornd(n, 0.1, 1, N);
     
  % Replace NaN 
    sum(walow(randnrs + 1) <= p & p <= waup(randnrs + 1))/N;
    sum(wilow(randnrs + 1) <= p & p <= wiup(randnrs + 1))/N;
    sum(cplow(randnrs + 1) <= p & p <= cpup(randnrs + 1))/N;
   
   
    % Save time 
     time_vec(tt) = toc()
end

% Average time 
mean(time_vec)

      
      