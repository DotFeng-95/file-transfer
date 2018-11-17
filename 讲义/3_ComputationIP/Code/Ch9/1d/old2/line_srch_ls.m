  function [tau,f] = line_srch(tau0,phi0,f0,g0,p,H,b);
  
%  Perform backtracking projected linesearch.

  mu = 0.25;
  tau = tau0;
  k = 0;
  k_max = 40;
  termcode = 0;
  
  while termcode == 0
    k = k + 1;
    f = f0 + tau*p;
    Pf = proj(f);
    phi = 0.5*Pf'*H*Pf - Pf'*b;
    psi = mu*g0'*(Pf - f0);
    if phi < phi0 + psi
      f = Pf;
      return
    elseif k > k_max
      disp('  *** Too many iterations in line search.')
      return
    else
      tau = 0.5 * tau;
    end
  end