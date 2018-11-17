  function [T,g] = cost_functional(u,cost_params)
%
%  Evaluate penalized least squares function 
%      T(u) = ||K*u - d||^2/2 + alpha * J_reg
%  and its gradient
%      g(u) = K'*(K*u-d) + alpha * grad J_reg
%  The penalty functional takes the form
%      J_reg = 2 * sum_i psi([Du]_i^2,beta) * Delta_x.
%  The gradient of the penalty functional takes the form
%      grad J_reg = D' * diag(psi'(Du^2,beta)) * D * Delta_x.

%  Get parameters and data.

  alpha = cost_params.reg_param;
  beta = cost_params.smoothing_param;
  Delta_x = cost_params.Delta_x;
  K = cost_params.smoothing_operator;
  D = cost_params.derivative_operator;
  d = cost_params.data_vec;

  Du = D*u;
  Du_sq = Du.^2;
  resid = K*u - d;
  J_reg = .5 * sum(psi(Du_sq,beta)) * Delta_x;
  grad_J_reg = D' * diag(psi_prime(Du_sq,beta)) * Du * Delta_x;
  T = norm(resid)^2 / 2 + alpha * J_reg;
  g = K'*resid + alpha * grad_J_reg;
