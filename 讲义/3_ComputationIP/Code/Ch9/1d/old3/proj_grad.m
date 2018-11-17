  function [pg] = proj_grad(g,f)
  
%  Compute projected gradient.

  indx0 = (f==0);
  if max(size(indx0)) > 0
    pg(indx0) = min(g(indx0),0);
  end
  
  indx1 = (f>0);
  if max(size(indx1)) > 0
    pg(indx1) = g(indx1);
  end
  