  function [tau,f] = line_srch_lhd(tau0,phi0,f0,g0,p, K,d,L,alpha,sigsq);
  
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
    Kf = K*Pf;
    phi = sum(Kf+sigsq - (d+sigsq).*log(Kf+sigsq)) + alpha/2*Pf'*L*Pf;
    psi = mu*g0'*(Pf - f0);
    if phi < phi0 - psi
      f = Pf;
      return
    elseif k > k_max
      disp('  *** Too many iterations in line search.')
      return
    else
      tau = 0.5 * tau;
    end
  end