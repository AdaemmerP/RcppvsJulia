
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
  walow     = zeros(N, 1); 
  waup      = zeros(N, 1); 
  wilow     = zeros(N, 1); 
  wiup      = zeros(N, 1); 
  cplow     = zeros(N, 1); 
  cpup      = zeros(N, 1); 
% Rpush('p', p, 'n', n, 'N', N)

  cpveclow  = zeros(n + 1, 1); 
  cpvecup   = zeros(n + 1, 1); 
 
% Compute lower and upper bounds for Clopper-Pearson interval
  for  ii = 1:(n + 1) 

      x             = ii - 1;
      cpveclow(ii)  = betainv(a1, x,     n - x + 1); 
      cpvecup(ii)   = betainv(a2, x + 1, n - x); 

  end

% Replace NaN 
  cpveclow(1)  = 0;
  cpvecup(end) = 1;


% Draw N times from Bin(n, p) with RCall
%  Rrun('randnrs <- rbinom(N, n, p)')
%  randnrs = Rpull('randnrs');
%  randnrs = double(randnrs);
  
% Draw from Bin(n, p) with Matlab function
  randnrs = binornd(n, 0.1, 1, N);

% Main loop to compute CIs    
  for jj = 1:N
       
      x             = randnrs(jj);   
     
     % Wald
      c0            = x/n;
      c1            = z*sqrt(c0*(1 - c0)/n);
      walow(jj)     = c0 - c1 ; 
      waup(jj)      = c0 + c1 ; 

	 % Wilson
      b0        = (1/(1+ z^2/n))*(c0 + (z^2)/(2*n));
      b1        = (z/(1 + (z^2/n)))*sqrt((c0*(1 - c0)/n) + z^2/(4*n^2)  );
      wilow(jj) = b0 - b1;
      wiup(jj)  = b0 + b1;
          
     % CP
      cplow(jj) = cpveclow(x + 1);
      cpup(jj)  = cpvecup(x + 1);
        
  end
     
    % Compute coverage rates
     sum(walow <= p & p <=  waup)/N;
     sum(wilow <= p & p <=  wiup)/N;
     sum(cplow <= p & p <=  cpup)/N;
     
    % Save time 
     time_vec(tt) = toc()
end

% Average time 
mean(time_vec)
      
      