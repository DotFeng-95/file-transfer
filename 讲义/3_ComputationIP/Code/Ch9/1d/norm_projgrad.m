  function [norm_Pg] = norm_projgrad(g,x)
  
%  Compute the Euclidean norm of the projected gradient at x,
%    Pg_i = g_i,        if x_i > 0
%    Pg_i = min(g_i,0), if x_i = 0
%  Here g denotes the gradient and the feasible set is {x | x_i >= 0}.

  indx0 = (x==0);
  if max(size(indx0)) > 0
    n1sq = norm( min(g(indx0), 0) )^2;
  else
    n1sq = 0;
  end
  
  indx1 = (x>0);
  if max(size(indx1)) > 0
    n2sq = norm( g(indx1) )^2;
  else
    n2sq = 0;
  end
  
  norm_Pg = sqrt(n1sq + n2sq);
  