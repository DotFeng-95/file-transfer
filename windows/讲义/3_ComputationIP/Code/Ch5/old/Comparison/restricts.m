  function [v] = restricts(u,p)   %  Symmetric (unscaled) restriction
%
%  [v] = restricts(u)
%  [v] = restricts(u,p)
%
%  Apply 2-D multilevel restriction based on piecewise constant 
%  heirarchical basis functions.
%
%  u is the grid function (2-D array) being restricted.
%  p is the number of levels of restriction.

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
  elseif 2^p > n
    disp(' *** Error in restrict.m.  2^p > size(u) ***');
    return
  end

  for level = 1:p
    nd2 = n/2;
    v = zeros(nd2,nd2);
    for i = 1:nd2
      vi = sum([u(2*i-1,:); u(2*i,:)]);
      v(i,:) = vi(1:2:n-1) + vi(2:2:n);
    end
    u = v;
    n = nd2;
  end
  