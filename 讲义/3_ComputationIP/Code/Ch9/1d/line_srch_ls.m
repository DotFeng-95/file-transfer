  function [f,phi,tau,ls_termcode] = line_srch_ls(f0,phi0,tau0,g0,p,H,b);
  
%  Perform backtracking projected line search. Compute an "acceptable"
%  minimizer tau for 
%      phi(tau) = J( P(f0 + tau*p) ).
%  Here P denotes the projection onto the feasible set f>=0, p is a
%  descent direction, and
%      J(v) = 0.5*v'*H*v - b'*v.

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
  
  %  Compute beta_1. This is the shortest step-length parameter which
  %  increases the active set.
  
  indx = (p<0) & (f0>0);
  if sum(indx(:)) == 0
    beta1 = 0;
  else
    tmp =  -f0(indx) ./ p(indx);
    beta1 = min(tmp);
    if tau0 <= beta1   %  Unconstrained step-len param is smaller than beta1
      tau = tau0;
      f = f0 + tau*p;
      phi = 0.5*f'*H*f - f'*b;
      ls_termcode = 3;
      return
    end
  end

  %  Line search iteration.
  
  while ls_termcode == 0
    ls_iter = ls_iter + 1;
    f = proj_nonneg(f0 + tau*p);
    phi = 0.5*f'*H*f - f'*b;
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
      if tau < beta1
	tau = beta1;
	f = f0 + tau*p;
	phi = 0.5*f'*H*f - f'*b;
	ls_termcode = 5;
	return
      end
    end
  end
  