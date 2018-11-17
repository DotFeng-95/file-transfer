  function [c,C] = optimal_circ(A)
  
%  Compute optimal circulant approximation C to n X n matrix A.
%    C = argmin {||B - A||_fro : B is n X n circulant}
%      = circulant(c)

  [m,n] = size(A);
  if m ~= n
    fprintf('\n *** Input A must be a square matrix.\n');
    return
  end
  
  c = zeros(n,1);
  c(1) = sum(diag(A));
  for j=1:n-1
    c(j+1) = sum(diag(A,-j)) + sum(diag(A,n-j));
  end
  c = c / n;

  if nargout == 2
    C = circulant(c);
  end
  
  