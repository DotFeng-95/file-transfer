%  Newton.m
%
%  For a fixed value of the regularization parameter alpha, compute 
%  the minimizer u_alpha of the penalized least squares functional
%       ||K*u - d||^2 / 2 + alpha * J(u)
%  using Primal Newton's method with a line search.
%  For beta > 0,
%       J(u) = sum_i psi'(|[D*u]_i|^2,beta) * Delta_x
%  is a smooth version of the discrete total variation functional
%       TV(u) = sum_i |[D*u]_i|^2 * Delta_x.

  alpha = input(' Regularization parameter alpha = ');
  beta = input(' Smoothing parameter beta = ');
  Newt_max = input(' Max. no. of Newton iterations = ');
  Newt_grad_tol = input(' Gradient stopping tolerance = ');

  GAconst = 1e-2;
  Delta_x = 1/n;
  LSiter = 10;
  
  %  Set up discretization of first derivative operator.
  
  D = spdiags([-ones(n-1,1) ones(n-1,1)], [0 1], n-1,n) / h;

  %  Store parameters and data in structure arrays.
  
  cost_params.size = n;
  cost_params.Delta_x = Delta_x;
  cost_params.smoothing_operator = K;
  cost_params.derivative_operator = D;
  cost_params.data_vec = d;
  cost_params.reg_param = alpha;
  cost_params.smoothing_param = beta;
  
  ls_params.max_ls_iter = LSiter;
  ls_params.GA_const = GAconst;
  ls_params.min_step = eps^(1/3);
  
  u_newt = zeros(n,1);
  newt_snorm = [];
  newt_gradnorm = [];
  Newton_conv_history = [];        %  Newton convergence history
  if exist('u_alpha','var')
    newt_enorm = norm(u_newt - u_alpha) / norm(u_alpha);
  end
  
%%% Newton iteration.
  
  for Newt_iter = 1:Newt_max
    
    %  Evaluate cost functional and its gradient.
    
    Du_sq = (D*u_newt).^2;
    [T,gradT] = cost_functional(u_newt,cost_params);
    normgrad = norm(gradT);

    %  Evaluate Hessian.
    
    diff_coef = psi_prime(Du_sq,beta) + 2*psi_doubleprime(Du_sq,beta).*Du_sq;
    H_reg = D' * diag(diff_coef) * D * Delta_x;

    H = K'*K + alpha*H_reg;
    p = -H \ gradT;

    %%% Perform backtracking line search.

    [unew,Tnew,gradT,LSflag,costfn_calls,lambda] = ...
	line_search(u_newt,p,T,gradT,cost_params,ls_params); 
    s = unew - u_newt;
    u_newt = unew;
    Newton_conv_history = [Newton_conv_history u_newt];
    T = Tnew;
    snorm = norm(s);
    gradnorm = norm(gradT);
    newt_snorm = [newt_snorm; snorm];
    newt_gradnorm = [newt_gradnorm; gradnorm];
    if exist('u_alpha','var')
      newt_enorm = [newt_enorm; norm(u_newt-u_alpha)/norm(u_alpha)];
    end
    
    %  Display information
    
    fprintf(' Newton iter = %4.0f, ||grad|| = %6.4e,  ||step|| = %6.4e.\n', ...
       Newt_iter, gradnorm, snorm);
   
    figure(1)
      plot(x,u_newt,'-',x,f_true,'--'), xlabel('x axis'), ylabel('u(x)'),
      title('Exact Solution (--) and TV Reconstruction')
    figure(2)
      subplot(221)
        semilogy(newt_gradnorm,'o')
        xlabel('Primal Newton Iteration')
        ylabel('||g(u^(k))||_1')
        title('Norm of Newton Gradient')
      subplot(222)
        semilogy(newt_snorm,'o')
        xlabel('Primal Newton Iteration')
        ylabel('||u^(k+1)-u^(k)||_1')
        title('Norm of Newton Step')
    if exist('u_alpha','var')
      subplot(223)
        semilogy(newt_enorm,'o')
        xlabel('Primal Newton Iteration')
	title('Relative Solution Error')
    end
    drawnow
      
    if gradnorm <= Newt_grad_tol, break, end

  end 					%%% End of Newton iteration loop.
    
  clear Newt_max;