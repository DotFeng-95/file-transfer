  function [v] = prolong(u)
%
%  [v] = prolong(u)
%
%  2-D prolongation operator corresponding to piecewise constant 
%  basis functions.

  nsq = max(size(u));
  n = sqrt(nsq);
  umat = reshape(u,n,n);
  vmat = zeros(2*n,2*n);
  for i = 1:n
    ui = umat(i,:);
    vi = [ui; ui];
    vi = vi(:)';
    vmat(2*i-1:2*i,:) = [vi; vi];
  end
  v = vmat(:);
  