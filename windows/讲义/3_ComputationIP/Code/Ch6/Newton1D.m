%  Newton1D.m
%
%  Use Newton or Gauss-Newton iterations to minimize the regularized
%  least squares functional,
%      J(q) = ||F(q) - d||^2 / 2 + alpha/2 * q'*L*q,
%  where L is the negative Laplacian with Neumann BC.

  fprintf('\n Run code once with lots of iterations to generate');
  fprintf('\n the "true" minimizer. Then run code again to obtain');
  fprintf('\n the iterative solution error.\n\n');
  
  alpha = 1e-2; %%%input(' Regularization parameter alpha = ');
  max_iter = input(' Max. no. of quasi-Newton iterations = ');
  qnflag = input(' Enter 0 for Newton iter; 1 for Gauss-Newton iter: ');
  nq = n_nodes + 1;
  L = laplacian(nq);
  q = ones(nq,1);    %  Initial guess.
  
  gnormvec = [];
  errorvec = [];
  if exist('q_alpha','var')==0
    q_alpha = q_true;  %  Set q_alpha to true q if q_alpha has not been set.
  end
  
  for k = 1:max_iter
    [g_ls,J_ls,A,u] = eval_gradls(q,b,d);
    if qnflag == 1
      H_ls = eval_HGN(q,b,d);  %  Gauss-Newton Hessian approximation.
    else
      H_ls = eval_Hess(g_ls,q,b,d);  %  True Hessian.
    end
    g = g_ls + alpha * L * q;
    H = H_ls + alpha * L;
    s = -H\g;
    q = q + s;
    
    figure(2)
      plot(x_mid,q,'-', x_mid,q_true,'--')
      title('True (--) and Estimated (-) Parameters')
    disp('  Hit any key to continue.')
    pause
    gnormvec = [gnormvec; norm(g)];
    errorvec = [errorvec; norm(q-q_alpha)/norm(q_alpha)];
  end

  figure(3)
    semilogy(gnormvec,'o-')
    xlabel('Quasi-Newton Iteration')
    title('Norm of Gradient')
    
  figure(4)
    semilogy(errorvec,'o-')
    xlabel('Quasi-Newton Iteration')
    title('Norm of Iterative Solution Error')
    
  q_alpha = q;   %  Reset q_alpha to current estimate for q_alpha.


