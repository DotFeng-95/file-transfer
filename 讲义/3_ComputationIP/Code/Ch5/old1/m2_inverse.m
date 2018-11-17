  function y_array = m2_inverse(x_array,params)

%  y_array = m2_inverse(x_array,params);
%
%  Compute array(y), where y = inverse(M2)*x. M2 = bccb(m2) is the 
%  block circulant-circulant block (BCCB) level 2 approximation to 
%  the matrix A = T'*T + alpa*L. The matrix inverse is computed with
%  O(n log n) complexity through the use of 2-D FFT's. 

  m2_hat = params.level2_fourier_representer;
  y_array = real(ifft2( fft2(x_array) ./ m2_hat ));
