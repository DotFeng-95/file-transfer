  function [v] = prolongs(u,p)   %  Symmetric (unscaled) prolongation.
%
%  [v] = prolongs(u)
%  [v] = prolongs(u,p)
%
%  Apply 2-D multilevel prolongation based on piecewise constant 
%  heirarchical basis functions.
%
%  u is the grid function (2-D array) being prolongated.
%  p is the number of levels of prolongation.

  if nargin == 1
    p = 1;
  end
  [m,n] = size(u);
  if m ~= n
    disp(' *** Error in restrict.m.  u must be square array ***');
    return 
  elseif 2^(round(log2(n))) ~= n
    disp(' *** Error in restrict.m.  size(u) must be power of 2 ***');
    return
  elseif p < 1
    disp(' *** Error in restrict.m.  p must be positive ***');
    return
  end
  
  for level = 1:p
    n2 = 2 * n;
    v = zeros(n2,n2);
    for i = 1:n
      ui = u(i,:);
      vi = [ui; ui];
      vi = vi(:)';
      v(2*i-1:2*i,:) = [vi; vi];
    end
    u = v;
    n = n2;
  end