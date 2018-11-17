  function [tau,f,termcode] = line_srch(tau0,phi0,f0,g0,p,H,b);
  
%  Perform backtracking projected linesearch.

  mu = 0.25;
  tau = tau0;
  k = 0;
  k_max = 40;
  termcode = 0;
  phi0_prime = g0'*p;
  if phi0_prime >= 0
    f = f0;
    termcode = 2;
    disp('  *** p is not descent direction in line search.')
    return
  end
  
  %  Compute beta_1.
  
  indx = (p<0) & (f0>0);
  if sum(indx)==0
    beta1 = 0;
  else
    tmp =  -f0(indx)./p(indx);
    beta1 = min(tmp);
    if beta1 >= tau0
      f = f0 + tau0*p;
      return
    end
  end

  %  Line search iteration.
  
  while termcode == 0
    k = k + 1;
    f = proj(f0 + tau*p);
    phi = 0.5*f'*H*f - f'*b;
    psi = mu*g0'*(f - f0);
    if phi < phi0 + psi
      return
    elseif k > k_max
      f = f0;
      termcode = 1;
      disp('  *** Too many iterations in line search.')
      return
    else    %  Quadratic backtrack.
      tau1 = -.5*phi0_prime*tau^2 / (phi - phi0 - phi0_prime*tau);
      tau = max(median([.1*tau,tau1,.5*tau]),beta1);
    end
  end