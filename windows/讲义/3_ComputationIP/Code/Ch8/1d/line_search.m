  function [x,Jnew,gnew,termcode,Jcalls,lambda] = ...
                         line_search(x,s,J,g,cost_params,ls_params)
		     
%  [x,Jnew,gnew,termcode,Jcalls,lambda] = ...
%                         line_search(x,s,J,g,cost_params,ls_params)
%
%  Minimize
%       y(lambda) = J(x + lambda*s)
%  using a quadratic backtracking line search. J(x) is the cost 
%  functional, g = grad J(x), and s is the search direction.

  maxiter = ls_params.max_ls_iter;  %  Max. no. line search iterations.
  GAconst = ls_params.GA_const;     %  Armejo constant.
  minstep = ls_params.min_step;     %  Minimum step size.

  lambda = 1;
  Jcalls = 0;   
  J0 = J;
  Jprime0 = g'*s;
  termcode = -1;
  if Jprime0 >= 0,
    disp(' Direction s is not a descent direction.');
    termcode = 1;
  end  
  
  while termcode < 0,
    Jcalls = Jcalls + 1;
    xnew = x + lambda*s;
    [Jnew,gnew] = cost_functional(xnew,cost_params);  
    
    %%% Check Armejio condition for sufficient decrease in J.
    
    if Jnew - J0 <= GAconst * Jprime0 * lambda,
      x = xnew;
      termcode = 0;
    else
      lambda = - Jprime0 * lambda^2 / (2*(Jnew - Jprime0*lambda - J0));
    end
    if Jcalls >= maxiter,
      termcode = 2;
      disp(' Max. no. iterations exceeded in line srch.')
      gnew = g;
      Jnew = J;
    end
  end
  
    
