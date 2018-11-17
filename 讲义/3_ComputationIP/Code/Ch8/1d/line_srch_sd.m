  function [x,Jnew,gnew,termcode] = ...
                         line_srch_sd(x,lambda,p,J,g,cost_params,ls_params)
%
%  Minimize
%       y(lambda) = J(x + lambda*p)
%  using a quadratic line search.
%  g = grad J(x) and p is the search direction.
%

  minstep = ls_params.min_step;  %  Minimum step size.
  GA1 = ls_params.GA_const(1);
  GA2 = ls_params.GA_const(2);
  maxiter = ls_params.max_ls_iter;

%  Initialization phase.

  Jcalls = 0;   
  J0 = J;
  Jprime0 = g'*p;
  if Jprime0 >= 0,
    disp(' Direction p is not a descent direction.');
    termcode = 1;
    initflag = -1;
  else
    initflag = 0;
  end  

  while initflag == 0,
    xnew = x + lambda * p;
    [Jnew,gnew] = cost_functional(xnew,cost_params);
    if Jnew >= J0 + GA1 * lambda * Jprime0;
      initflag = 1;
      termcode = -1;
    elseif gnew'*p >= GA2 * Jprime0,
      initflag = 1;
      termcode = 0;
      x = xnew;
    else
      lambda = lambda * 2;
    end
  end

%  Backtracking phase.  

  lambda = - Jprime0 * lambda^2 / (2*(Jnew - Jprime0*lambda - J0));
  while termcode < 0,
    xnew = x + lambda*p;
    [Jnew,gnew] = cost_functional(xnew,cost_params);  
    Jcalls = Jcalls + 1;
    
    %%% Check Goldstein Armejio condition.
    
    if Jnew - J0 <= GA1 * Jprime0 * lambda,
      x = xnew;
      termcode = 0;
    else
      lambda = - Jprime0 * lambda^2 / (2*(Jnew - Jprime0*lambda - J0));
    end
    if Jcalls >= maxiter,
      termcode = 2;
      disp(' Max. no. iterations exceeded in line srch.')
    end
  end
    
