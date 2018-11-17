  function y_array = a_mult(x_array,params);

%
%
%  Compute array(y), where y = (T'*T + alpha*L)*x. T is assumed to be
%  block Toeplitz with Toeplitz blocks (BTTB). Block circulant
%  extensions are used to compute matrix-vector products involving T
%  and T'. The products have O(n log n) complexity through the use of
%  2-D FFT's. 

  t_ext_hat = params.t_extension_hat;
  alpha = params.reg_param;
  L = params.penalty_matrix;
  nx = params.array_size(1);
  ny = params.array_size(2);
  nx2 = 2*nx;
  ny2 = 2*ny;
  
  %  Compute T'*(T*x) via circulant extension.
  
  tmp = zeros(nx2,ny2);
  tmp(1:nx,1:ny) = x_array;
  Tx_array = real(ifft2( t_ext_hat .* fft2(tmp) ));
  tmp = zeros(nx2,ny2);
  tmp(1:nx,1:ny) = Tx_array(1:nx,1:ny);
  y_array = real(ifft2( conj(t_ext_hat) .* fft2(tmp) ));
  y_array = y_array(1:nx,1:ny);
  
  %  Add on alpha*L*x.
  
  Lx_array = reshape(L*x_array(:),nx,ny);
  y_array = y_array + alpha*Lx_array;
  