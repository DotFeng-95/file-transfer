  function [g,J,A,u] = eval_gradls(q,b,d)

%  Evaluate gradient of least squares cost functional using costate,
%  or adjoint, method.

  nq = max(size(q));
  [n,m] = size(b);
  g = zeros(nq,m);
  h = 1/nq;
  
  [J,A,u] = eval_Jls(q,b,d);
  resid = u - d;
  z = -A \ resid;
  Du = diff(u)/h;
  g(1,:) = u(1,:).*z(1,:) / h;
  for i = 2:n
    g(i,:) = -Du(i-1,:).*z(i-1,:) + Du(i-1,:).*z(i,:);
  end
  g(n+1,:) = u(n,:).*z(n,:) / h;
  
  g = sum(g')';   %  Sum along rows.