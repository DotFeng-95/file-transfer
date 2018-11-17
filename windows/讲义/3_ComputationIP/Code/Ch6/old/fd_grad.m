%  fd_grad.m
%
%  Evaluate gradient of least squares cost functional using finite
%  differences.

  J0 = eval_Jls(q0,b,d);
  g_fd = zeros(nq,1);
  for i = 1:nq
    tau = eps^(1/3) * max(abs(q0(i)),1);
    qi = q0;
    qi(i) = q0(i) + tau;
    Ji = eval_Jls(qi,b,d);
    g_fd(i) = (Ji - J0) / tau;
  end

  figure(1)
    plot(x_mid,g_fd)
    title('Finite Difference Gradient')
