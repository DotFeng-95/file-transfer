  function y = proj(x)
  
%  Compute projection of x onto C = {x | x >= 0}.
  
  y = max(x,0);