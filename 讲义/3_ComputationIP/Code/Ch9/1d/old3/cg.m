  function f = cg(f0,H,b,max_iter_CG);
  
%  Use CG iteration to solve H*f = b, or equivalently, to minimize
%      J(f) = 0.5*f'*H*f - f'*b.

  g = H*f - b;
  d = -g;
  delta = norm(g)^2;
  cg_iter = 0;
  
  while term_code == 0
    cg_iter = cg_iter + 1;
    Hd = H*d;
    dHd = d(:)'*Hd(:);              %  Compute d'*A*d.
    tau = delta / dHd;              %  Line search parameter.
    f = f + tau*d;                  %  Update approximate solution.
    g = g + tau*Ad;                 %  Update gradient g = H*f-b.

    delta_new = norm(g)^2;
    my_beta = delta_new / delta;    %  Note: beta is a MATLAB function.
    d = -g + my_beta*d;             %  Update CG search direction.
    delta = delta_new;
  