  function y = m_inv(b,params)
  
  A = params.matrix;
  y = (A + .001*eye(size(A)))\b;
  