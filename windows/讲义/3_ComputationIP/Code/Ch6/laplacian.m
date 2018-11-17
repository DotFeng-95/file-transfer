  function A = laplacian(n)

%  Set up 1-D discrete Laplacian with Neumann B.C.

  Asub = -ones(n,1);
  Adiag = 2*ones(n,1);
  Adiag(1) = 1;
  Adiag(n) = 1;
  A = spdiags([Asub Adiag Asub], [-1 0 1], n,n);