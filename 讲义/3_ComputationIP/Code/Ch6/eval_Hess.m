  function Hess = eval_Hess(g0,q0,b,d)
  
%  Evaluate Hessian of least squares functional.  
  
  nq = max(size(q0));
  Hess = zeros(nq,nq);
  g0 = eval_gradls(q0,b,d);
  for i = 1:nq
    tau = eps^(1/3) * max(abs(q0(i)),1);
    qi = q0;
    qi(i) = q0(i) + tau;
    gi = eval_gradls(qi,b,d);
    Hess(:,i) = (gi - g0) / tau;
  end
  