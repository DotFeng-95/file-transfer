  function y_array = circ_ext_inverse(x_array,params);

%  y_array = circ_ext_inverse(x_array,params);
%
%  Compute array(y), where y = inverse(M)*x. M is the matrix
%  representer for the block circulant extension preconditioner.
%  The matrix inverse is computed with O(n log n) complexity through
%  the use of 2-D FFT's.  

  m_inv_hat = params.precond_representer;
  nx = params.array_size(1);
  ny = params.array_size(2);
  nx2 = 2*nx;
  ny2 = 2*ny;
  
  tmp = zeros(nx2,ny2);
  tmp(1:nx,1:ny) = x_array;
  y_array = real(ifft2( m_inv_hat .* fft2(tmp) ));
  y_array = y_array(1:nx,1:ny);
