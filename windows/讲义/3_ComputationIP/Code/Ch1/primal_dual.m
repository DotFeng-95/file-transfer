  function [u_PD,pd_gradnorm] = primal_dual(alpha,beta,PDnewt_iter,K,d,ioflag)

%  Use primal-dual Newton's method to minimize the TV-penalized least
%  squares functional
%      T(u) = ||K*u - d||^2/2 + alpha*J(u).
%  Here K is a discretized integral operator, d is discrete data, 
%  ||.|| denotes the l^2 norm, alpha is a positive regularization 
%  parameter, and J is a smooth approximation to the 
%  Total Variation functional.
%      J(u) = sum_i 2*psi(|[D*u]_i|^2,beta),
%  where D is a discretization of the first derivative operator and
%  beta is a positive smoothing parameter.
%  The primal-dual system
%      g1 = K'*(Ku-d) + alpha*D'*v         = 0
%      g2 = D*u - psi'(|[D*u]_i|^2,beta)*v = 0
%      max_i |v_i| <= 1
%  is solved using Newton's method. At each iteration, the size of 
%  the v-step is controlled with a line search to ensure that the 
%  constraint max_i |v_i| <= 1 is maintained.

  %  Set up discretization of first derivative operator.
  
  n = max(size(d));
  x = [0:1/(n-1):1]';
  Delta_x = 1/n;
  D = spdiags([-ones(n-1,1) ones(n-1,1)], [0 1], n-1,n) / Delta_x;

  %  Initialize primal variable u_PD and dual variable v.
  
  u_PD = zeros(n,1);
  v = zeros(n-1,1);
  pd_gradnorm = [];
  snorm_vec = [];
  
  %  Newton iteration.
  
  for k = 1:PDnewt_iter
    Du = D*u_PD;
    Du_sq = Du.^2;
    psi_1 = 1 ./ sqrt(Du_sq + beta^2);
    psi_2 = -.5 * (Du_sq + beta^2).^(-1.5);
    B = diag(psi_1);
      
    %  Compute components G1, G2 of nonlinear system.
      
    G1 = K'*(K*u_PD - d) + alpha*D'*v;
    G2 = Du - B\v;
    G = [G1; G2];
    gradnorm = norm(G);
    pd_gradnorm = [pd_gradnorm; gradnorm];

    %  Compute primal-dual Newton step.
      
    M = diag(1 + 2*psi_2 ./ psi_1.^2 .* Du .* v);
    L = D' * B * M * D;
    rhs = -K'*(K*u_PD - d) - alpha * D' * B * Du;
    Delta_u = (K'*K + alpha*L) \ rhs;
    Delta_v = -v + B * (Du + M*D*Delta_u);
    
    %  Update primal variable
    
    u_PD = u_PD + Delta_u;
    
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
    
    if ioflag
      fprintf(' PD Newt iter=%4.0f, ||grad||=%6.4e,  ||step||=%6.4e.\n', ...
        k, gradnorm, snorm);
      figure(1)
        plot(x,u_PD,'-')
        xlabel('x axis')
        title('TV Regularized Solution (-)')
      
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
      drawnow;
      pause(.1);
    end
    
    if gradnorm < 1e-12, break, end
    
  end  %  for k
