  function Ax = mat_mult(x,params)

%  Ax = mat_mult(x,params)
%
%  Compute matrix-vector product A*x, where the matrix A is 
%  stored as A = params.matrix.
  
  A = params.matrix;
  Ax = A*x;