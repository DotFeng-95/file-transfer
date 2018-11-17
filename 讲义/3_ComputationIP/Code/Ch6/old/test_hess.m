%  test_hess.m
%

%  Evaluate Hessian of least squares cost functional by finite
%  differencing the gradient.

  q0 = q_true;
  nq = max(size(q0));
  Hess = zeros(nq,nq);
  [g0,J0,A0,u0] = eval_gradls(q0,b,d);
  for i = 1:nq
    tau = eps^(1/3) * max(abs(q0(i)),1);
    qi = q0;
    qi(i) = q0(i) + tau;
    gi = eval_gradls(qi,b,d);
    Hess(:,i) = (gi - g0) / tau;
  end
  
%  Gauss-Newton approximation to the Hessian.

  H_GN = zeros(nq,nq);
  [J0,A0,u0] = eval_Jls(q0,b,d);
  Du = diff(u0)/h;
  for i = 1:nq
    Aiu = zeros(n,1);
    if i == 1
      Aiu(1) = u0(1)/h;
    elseif i == nq
      Aiu(n) = u0(n)/h;
    else
      Aiu(i-1) = -Du(i-1);
      Aiu(i)   = Du(i-1);
    end
    si = A0 \ (A0 \ Aiu);

    H_GN(1,i) = u0(1) * si(1) / h;
    for j = 2:n
      H_GN(j,i) = -Du(j-1)*si(j-1) + Du(j-1)*si(j);
    end
    H_GN(n+1,i) = u0(n) * si(n) / h;
  end
      
