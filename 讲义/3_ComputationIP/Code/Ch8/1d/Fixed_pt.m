%  Fixed_pt.m
%
%  Use "lagged diffusivity" fixed point iteration to minimize
%      T(u) = ||K*u - d||^2/2 + alpha*J(u),
%  where K is a discretized integral operator, d is discrete data, 
%  ||.|| denotes the l^2 norm, alpha is a positive regularization 
%  parameter, and J is a smooth approximation to the 
%  Total Variation functional.
%      J(u) = sum_i psi(|[D*u]_i|^2,beta) * Delta_x,
%  where D is a discretization of the first derivative operator and
%  beta is a positive smoothing parameter.
%
%  At each iteration, replace u by u + Delta_u, where Delta_u solves
%    (K'*K + alpha*L(u)) * Delta_u = K'(K*u-d) + alpha*L(u)*u,
%  where 
%    L(u)*v = D'* diag(psi'(|[D*u]_i|^2,beta) * D * Delta_x. 

  alpha = input(' Regularization parameter alpha = ');
  beta = input(' TV smoothing parameter beta = ');
  fp_iter = input(' No. of fixed point iterations = ');

  %  Set up discretization of first derivative operator.
  
  D = spdiags([-ones(n-1,1) ones(n-1,1)], [0 1], n-1,n) / h;

  %  Initialization.

  Delta_x = 1 / n;
  fp_gradnorm = [];
  snorm_vec = [];
  FPconv_history = [];
  u_fp = zeros(n,1);
  if exist('u_alpha','var')
      fp_enorm = norm(u_fp-u_alpha)/norm(u_alpha);
  end
  
  for k = 1:fp_iter
      
    %  Set up regularization operator L.

    Du_sq = (D*u_fp).^2;
    L = D' * diag(psi_prime(Du_sq,beta)) * D * Delta_x;
    H = K'*K + alpha*L;
    g = H*u_fp - K'*d;
    s = -H \ g;
    u_fp = u_fp + s;
    FPconv_history = [FPconv_history u_fp];
    
    snorm = norm(s);
    gradnorm = norm(g);
    snorm_vec = [snorm_vec; snorm];
    fp_gradnorm = [fp_gradnorm; gradnorm];
    if exist('u_alpha','var')
      fp_enorm = [fp_enorm; norm(u_fp-u_alpha)/norm(u_alpha)];
    end
    
    %  Display solution and gradient norm
 
    fprintf(' FP iter = %3.0f, ||gradient|| = %6.4e,  ||step|| = %6.4e.\n', ...
       k, gradnorm, snorm);
   
    figure(1)
      plot(x,f_true,'--', x,u_fp,'-')
      xlabel('x axis')
      title('True Solution (--) and TV Regularized Solution (-)')
    figure(2)
      subplot(221)
        semilogy([1:k],fp_gradnorm,'o', [1:k],fp_gradnorm,'-')
        xlabel('Fixed Point Iteration')
        title('Norm of FP Gradient')
      subplot(222)
        semilogy(snorm_vec,'o-')
        xlabel('Fixed Point Iteration')
        title('Norm of FP Step')
    if exist('u_alpha','var')
      subplot(223)
        semilogy(fp_enorm,'o-')
        xlabel('Fixed Point Iteration')
        title('Relative Solution Error')
    end

    pause(.1);
    drawnow
      
  end  %  for fp_iter
      
  clear fp_iter;