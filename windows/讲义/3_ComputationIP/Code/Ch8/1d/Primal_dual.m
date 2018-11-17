%  Primal_dual.m
%
%  Use primal-dual Newton's method to minimize the functional
%      T(u) = ||K*u - d||^2/2 + alpha*J(u),
%  where K is a discretized integral operator, d is discrete data, 
%  ||.|| denotes the l^2 norm, alpha is a positive regularization 
%  parameter, and J is a smooth approximation to the 
%  Total Variation functional.
%      J(u) = sum_i 2*psi(|[D*u]_i|^2,beta) * Delta_x,
%  where D is a discretization of the first derivative operator and
%  beta is a positive smoothing parameter.
%  The primal-dual system
%      g1 = K'*(Ku-d) + alpha*D'*v         = 0
%      g2 = D*u - psi'(|[D*u]_i|^2,beta)*v = 0
%      max_i |v_i| <= 1
%  is solved using Newton's method. At each iteration, the size of 
%  the v-step is controlled with a line search to ensure that the 
%  constraint max_i |v_i| <= 1 is maintained.

  alpha = input(' Regularization parameter alpha = ');
  beta = input(' TV smoothing parameter beta = ');
  PDnewt_iter = input(' Max. no. of primal-dual Newton iterations = ');

  %  Set up discretization of first derivative operator.
  
  Delta_x = 1/n;
  D = spdiags([-ones(n-1,1) ones(n-1,1)], [0 1], n-1,n) / Delta_x;

  %  Initialize primal variable u_PD and dual variable v.
  
  u_PD = zeros(n,1);
  v = zeros(n-1,1);
  pd_gradnorm = [];
  snorm_vec = [];
  if exist('u_alpha','var')
      pd_enorm = norm(u_PD-u_alpha)/norm(u_alpha);
  end
  PDconv_history = [];        %  Primal-dual convergence history
  
  %  Newton iteration.
  
  for k = 1:PDnewt_iter
    Du = D*u_PD;
    Du_sq = Du.^2;
    psi_1 = psi_prime(Du_sq,beta);
    psi_2 = psi_doubleprime(Du_sq,beta);
    B = diag(psi_1);
      
    %  Compute components G1, G2 of nonlinear system.
      
    G1 = K'*(K*u_PD - d) + alpha*D'*v * Delta_x;
    G2 = Du - B\v;
    G = [G1; G2];
    gradnorm = norm(G);
    pd_gradnorm = [pd_gradnorm; gradnorm];

    %  Compute Jacobian and compute Newton step.
      
%    J11 = K'*K;
%    J12 = alpha*D' * Delta_x;
%    J21 = M*D;
%    J22 = -inv(B);
%    J = [J11 J12; J21 J22];
%    Delta = -J\G;
%    Delta_u = Delta(1:n);
%    Delta_v = Delta(n+1:2*n-1);
    
    M = diag(1 + 2*psi_2 ./ psi_1.^2 .* Du .* v);
    L = D' * B * M * D * Delta_x;
    rhs = -K'*(K*u_PD - d) - alpha * Delta_x * D' * B * Du;
    Delta_u = (K'*K + alpha*L) \ rhs;
    Delta_v = -v + B * (Du + M*D*Delta_u);
    
    %  Update primal variable
    
    u_PD = u_PD + Delta_u;
    PDconv_history = [PDconv_history u_PD];
    
    %  Perform line search in the dual variable v.  Calculate
    %  minimum rho in (0,1] for which |v(i) + rho*Delta_v(i)| = 1.
      
    rho = (1 - sign(Delta_v).*v) ./ (abs(Delta_v) + eps);
    rho_min = min(min(rho(:)),1);
    if rho_min < 1
      rho_min = .9*rho_min;
    end
    v = v + rho_min*Delta_v;
      
    %  Display results.
    
    snorm = sqrt(norm(Delta_u)^2 + norm(Delta_v)^2);
    snorm_vec = [snorm_vec; snorm];
    if exist('u_alpha','var')
      pd_enorm = [pd_enorm; norm(u_PD-u_alpha)/norm(u_alpha)];
    end
    fprintf(' PD Newt iter=%4.0f, ||grad||=%6.4e,  ||step||=%6.4e.\n', ...
       k, gradnorm, snorm);
    
    figure(1)
      plot(x,f_true,'--', x,u_PD,'-')
      xlabel('x axis')
      title('True Solution (--) and TV Regularized Solution (-)')
      
    figure(2)
      indx = [1:max(size(pd_gradnorm))]';
      subplot(221)
        semilogy(indx,pd_gradnorm,'o', indx,pd_gradnorm,'-')
        xlabel('Primal-Dual Newton Iteration')
        title('Norm of Primal-Dual Gradient')
      subplot(222)
        semilogy(indx,snorm_vec,'o', indx,snorm_vec,'-')
        xlabel('Primal-Dual Newton Iteration')
        title('Norm of Primal-Dual Step')
    if exist('u_alpha','var')
      subplot(223)
        semilogy(pd_enorm,'o-')
        xlabel('Primal/Dual Newton Iteration')
        title('Relative Solution Error')
    end
    drawnow;
    pause(.1);

    if gradnorm < grad_tol, break, end
    
  end  %  for k
  
  clear PDnewt_iter;