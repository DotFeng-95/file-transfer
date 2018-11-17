%  Demo_fourier.m
%
%  Demonstrate connection between the Fourier matrix and 
%  the discrete Fourier transform.

  nx = 5;
  hx = 1/nx;
  x = [0:hx:1-hx]';
  f = sin(3*pi*x);
  
  %  Construct sqrt(nx)*F, where F is the nx X nx Fourier matrix.

  Fx = sqrt(nx) * fourier(nx);
  
  %  Compare Fx*f with FFT(f).
  
  Fxf = Fx*f;
  FFTf = fft(f);
  Fxf_minus_FFTf = norm(Fxf-FFTf) / norm(FFTf)
  
  %  Construct ny X ny Fourier matrix, which is the Kronacker
  %  product of 1-D Fourier matrices.
  
  ny = 4;
  Fy = sqrt(ny) * fourier(ny);
  F = kron(Fy,Fx);
  
  %  Construct 2-D "image" tau.
  
  hy = 1/ny;
  y = [0:hy:1-hy]';
  [xx,yy] = meshgrid(y,x);   %  size(xx) = size(yy) = [nx,ny]
  
  tau = 2*xx + 3*yy - xx.*yy + xx.^2 -4*yy.^2;
  
  Ftau = F*tau(:);
  
  %  Compute DFT using the FFT.
  
  dft_tau = fft2(tau);        %  This is a 2-D array.
  vec_dft_tau = dft_tau(:);  %  This is a column vector.

  %  Compare.
  
  Ftau_minus_dft_tau = norm(Ftau - vec_dft_tau) / norm(dft_tau)
  