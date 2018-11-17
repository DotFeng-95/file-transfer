  function y = mat_inv(x,params)

%  y = mat_inv(x,params)
%
%  Compute matrix-vector product M^{-1}*x, where the matrix M is 
%  stored as M = params.preconditioner.
  
  M = params.preconditioner;
  y = M\x;