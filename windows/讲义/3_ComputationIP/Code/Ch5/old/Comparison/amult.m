  function Av = amult(v,k_hat,alpha,L)
%
%  Compute A*v, where A = K'*K + alpha*L.

  [nx,ny] = size(v);
  Kv = convolve_2d(k_hat,v);
  Kv = Kv(1:nx,1:ny);
  KstarKv = convolve_2d(conj(k_hat),Kv);
  Av = KstarKv(1:nx,1:ny) + alpha*reshape(L*v(:),nx,ny);