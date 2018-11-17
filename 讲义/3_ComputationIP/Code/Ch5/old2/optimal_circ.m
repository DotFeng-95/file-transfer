  function [c,C] = optimal_circ(A)
  
%  c = optimal_circ(A)
%  [c,C] = optimal_circ(A)
%  
%  Compute circulant approximation C to square matrix A
%  which minimizes the Frobenius norm.
%    C = argmin {||B - A||_fro : B is n X n & circulant}
%  c is the representer for C, i.e., C = circulant(c).

%  A should be a square matrix.

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
  
  