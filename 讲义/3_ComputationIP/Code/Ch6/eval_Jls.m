  function [J_ls,A,u] = eval_Jls(q,b,d)

%  Evaluate parameter-to-solution map u(q) = inv(A(q)) * b. Then 
%  evaluate the least squares cost functional,
%      J_ls(q) = ||u(q) - d||^2 / 2.

%  Assemble stiffness matrix A(q).

  [n,m] = size(b);
  h = 1 / (n+1);
  Adiag = (q(1:n) + q(2:n+1)) / h;
  Asub = -q(2:n) / h;
  Asuper = Asub;
  A = spdiags([[Asub;0] Adiag [0;Asuper]], [-1 0 1], n,n);

%  Solve systems A*u(:,l) = b(:,l), l = 1,...,m.

  u = A\b;

%  Evaluate least squares cost functional.

  J_ls = norm(u-d,'fro')^2/2;
