  function [c,C] = optimal_circ_toep(t)

%  c = optimal_circ_toep(t)
%  [c,C] = optimal_circ_toep(t)
%
%  Compute the circulant approximation C to a square Toeplitz 
%  matrix T which minimizes the Frobenius norm. 
%       C = argmin {||B - T||_fro : B is n X n & circulant}
%  t is the representer for T, i.e., T = my_toeplitz(t).
%  c is the representer for C, i.e., C = circulant(c).

%  t should be a vector of odd dimension.

  [m,n] = size(t);
  if min(m,n) ~= 1
    fprintf('\n *** Input t must be a 1-d array.\n');
    return
  end
  nn = max(m,n);
  if mod(nn,2) == 0
    fprintf('\n *** Length of input t must be odd.\n');
    return
  end

  n = ceil(nn/2);
  c = zeros(n,1);
  c(1) = n*t(n);
  for j=1:n-1
    c(j+1) = (n-j)*t(n+j) + j*t(j);
  end
  c = c / n;

  if nargout == 2
    C = circulant(c);
  end
  
  