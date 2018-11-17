%  Steepest_d.m
%
%  Use the Steepest Descent method to minimize the functional
%      T(u) = ||K*u - d||^2/2 + alpha*J(u),
%  where K is a discretized integral operator, d is discrete data, 
%  ||.|| denotes the l^2 norm, alpha is a positive regularization 
%  parameter, and J is a smooth approximation to the 
%  Total Variation functional.
%      J(u) = sum_i 2*psi(|[D*u]_i|^2,beta) * Delta_x,
%  where D is a discretization of the first derivative operator and
%  beta is a positive smoothing parameter.
   
%  Initialization

  alpha = input(' Regularization parameter alpha = ');
  beta = input(' Smoothing parameter beta = ');
  SD_max = input(' Max. no. of steepest descent iterations = ');

  GA1 = .1;
  GA2 = .5;
  nx = n;
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
  ls_params.GA_const = [GA1; GA2];
  ls_params.min_step = eps^(1/3);
  
  u_sd = zeros(n,1);
  sd_snorm = [];
  if exist('u_alpha','var')
    sd_enorm = norm(u_sd - u_alpha) / norm(u_alpha);
  end
  sd_gradnorm = [];
  SD_conv_history = [];
  [T,g] = cost_functional(u_sd,cost_params);
  lambda = 1;
  
%  Steepest descent iteration.
  
  for SD_iter = 1:SD_max

    p = -g/norm(g);
    [unew,Tnew,g,termcode] = ...
	line_srch_sd(u_sd,lambda,p,T,g,cost_params,ls_params);
    SD_conv_history = [SD_conv_history u_sd];
    du = unew - u_sd;
    tau = norm(du);
    gradnorm = norm(g);
    snorm = norm(du);
    sd_gradnorm = [sd_gradnorm; gradnorm];
    sd_snorm = [sd_snorm; snorm];
    u_sd = unew;
    T = Tnew;
    if exist('u_alpha','var')
      sd_enorm = [sd_enorm; norm(u_sd-u_alpha)/norm(u_alpha)];
    end
      
    %  Display info.
    
    fprintf(' SD iter = %4.0f, ||grad|| = %6.4e,  ||step|| = %6.4e.\n', ...
       SD_iter, gradnorm, snorm);
    
    figure(1)
      plot(x,u_sd,x,f_true,'--'), 
      title('Exact Soln (--) and TV Reconstruction (-)')
      xlabel('x axis'), ylabel('u(x)')
    figure(2)
      subplot(221)
        semilogy(sd_gradnorm,'o')
        xlabel('Steepest Descent Iteration')
        title('Norm of Steepest Descent Gradient')
      subplot(222)
        semilogy(sd_snorm,'o')
        xlabel('Steepest Descent Iteration')
        title('Norm of Steepest Descent Step')
    if exist('u_alpha','var')
      subplot(223)
        semilogy(sd_enorm,'o')
        xlabel('Steepest Descent Iteration')
        title('Norm of Relative Solution Error')
    end

  end          %%% End of Steepest Descent loop
    
  clear SD_max;