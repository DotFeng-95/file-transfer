  function [f,k,termcode] = grad_proj_lhd(f,K,d,L,alpha,sigsq,max_iter)
  
%  Apply gradient projection to compute minimizer of quadratic
%    J(f) = 0.5 * f'*H*f - f'*b
%  subject to nonnegativity constraints, f >= 0.

  sigma = 0.25;
  Kf = K*f;
  J0 = sum(Kf+sigsq - (d+sigsq).*log(Kf+sigsq)) + alpha/2*f'*L*f;
  max_delJ = 0;
  Active0 = (f==0);  %  Current active set of indices.
  k = 0;
  termcode = 0;
  
  while termcode == 0
    k = k + 1;    %  Iteration count.
    g = K'*diag(1./(Kf+sigsq))*(Kf-d) + alpha*L*f;
    H = K'*diag((d+sigsq)./(Kf+sigsq).^2)*K + alpha*L;
    p = -g;
    tau0 = norm(p)^2 / (p'*H*p);
    [tau,f] = line_srch_lhd(tau0,J0,f,g,p, K,d,L,alpha,sigsq);
    Kf = K*f;
    J = sum(Kf+sigsq - (d+sigsq).*log(Kf+sigsq)) + alpha/2*f'*L*f;
    delta_J = J0 - J;
    max_delJ = max(delta_J,max_delJ);
    Active = (f==0);
    
    %  Check stopping criteria
    
    if delta_J < sigma * max_delJ    %  Insufficient decrease in J
      termcode = 1;
      return
    elseif Active == Active0         %  Active set unchanged.
      termcode = 2;
      return
    elseif k >= max_iter             %  Max iteration count exceeded.
      termcode = 3;
      return
    end

    J0 = J;
    Active0 = Active;
    
  end