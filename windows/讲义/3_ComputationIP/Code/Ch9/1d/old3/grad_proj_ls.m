  function [f,k,termcode] = grad_proj_ls(f,H,b,max_iter)
  
%  Apply gradient projection to compute minimizer of quadratic
%    J(f) = 0.5 * f'*H*f - f'*b
%  subject to nonnegativity constraints, f >= 0.

  sigma = 0.25;
  J0 = 0.5 * f'*H*f - f'*b;
  max_delJ = 0;
  Active0 = (f==0);  %  Current active set of indices.
  k = 0;
  termcode = 0;
  
  while termcode == 0
    k = k + 1;    %  Iteration count.
    g = H*f - b;
    p = -g.*((1-Active0)+Active0.*(g<0));
    
    tau0 = norm(p)^2 / (p'*H*p);
    [tau,f,ls_termcode] = line_srch_ls(tau0,J0,f,g,p,H,b);
    J = 0.5 * f'*H*f - f'*b;
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
    elseif ls_termcode ~= 0          %  Line search failed.
      termcode = -ls_termcode;
      return
    end

    J0 = J;
    Active0 = Active;
    
  end