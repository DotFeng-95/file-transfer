%  costate.m

%  Evaluate gradient of least squares cost functional using costate,
%  or adjoint, method.

  [J0,A0,u_q0] = evalJ(q0,b,d);
  resid = u_q0 - d;
  z = -A0 \ resid;
  g_costate = zeros(n+1,1);
  
  Du = diff(u_q0)/h;
  g_costate(1) = u_q0(1) * z(1) / h;
  for i = 2:n
    g_costate(i) = -Du(i-1)*z(i-1) + Du(i-1)*z(i);
  end
  g_costate(n+1) = u_q0(n) * z(n) / h;
  
  figure(2)
    plot(x_mid,g_costate)
    title('Costate Gradient')
