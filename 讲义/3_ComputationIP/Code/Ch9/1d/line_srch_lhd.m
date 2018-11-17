  function [f,phi,tau,ls_termcode] = line_srch_lhd(f0,phi0,tau0,g0,p, ...
    K,d,L,alpha,sigsq)
  
%  Perform backtracking projected line search. Compute an "acceptable"
%  minimizer tau for 
%      phi(tau) = J( P(f0 + tau*p) ).
%  Here P denotes the projection onto the feasible set f>=0, p is a
%  descent direction, and
%      J(v) = sum(K*v+sigsq - (d+sigsq).*log(K*v+sigsq)) + alpha/2 * v'*L*v

  mu = 0.25;
  tau = tau0;
  ls_iter = 0;
  ls_iter_max = 40;
  ls_termcode = 0;
  phi0_prime = g0'*p;
  if phi0_prime >= 0
    f = f0;
    tau = 0;
    phi = phi0;
    ls_termcode = 2;
    disp('  *** p is not descent direction in line search.')
    return
  end
  
  %  Line search iteration.
  
  while ls_termcode == 0
    ls_iter = ls_iter + 1;
    f = proj_nonneg(f0 + tau*p);
    Kf = K*f;
    phi = sum(Kf+sigsq - (d+sigsq).*log(Kf+sigsq)) + alpha/2*f'*L*f;
    psi = mu*g0'*(f - f0);  %  More'-Toraldo sufficient decrease
%%%    psi = -mu / tau * norm(f-f0)^2;  %  Bertsekas sufficient decrease
    if phi < phi0 + psi
      ls_termcode = 1;    %  Normal line search termination
      return
    elseif ls_iter > ls_iter_max
      f = f0;
      phi = phi0;
      tau = 0;
      ls_termcode = 4;
      disp('  *** Too many iterations in line search.')
      return
    else    %  Quadratic backtrack.
      tau1 = -.5*phi0_prime*tau^2 / (phi - phi0 - phi0_prime*tau);
      tau = median([.1*tau,tau1,.5*tau]);
    end
  end
  