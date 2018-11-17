%  Deblur_fourier.m
%
%  Implement Tikhonov regularization using Fourier representations. 
%  Solve
%    (K'*K + alpha*L)f = K'*d,
%  where d denotes image data, K is the blurring operator, L is the
%  regularization operator, and alpha is a positive regularization 
%  parameter. It is assumed that BCCB(PSF) is the matrix representing 
%  K, with lexicographical column ordering of the unknowns. L can be
%  selected to be either the identity or the discrete Laplacian with
%  periodic boundary conditions.

  k_hat = fft2(PSF);
  d_hat = fft2(dat);
  alpha = input(' Regularization parameter alpha = ');
  Lflag = input(' Regularization type (0=I, 1=discrete Laplacian): ');
  
%  Construct Fourier representer for L.

  if Lflag == 0
    lhat = ones(nx,ny);
  elseif Lflag == 1
    l = zeros(nx,ny);
    l(1, 1) =  4;
    l(2 ,1) = -1;
    l(nx,1) = -1;
    l(1 ,2) = -1;
    l(1,ny) = -1;
    lhat = fft2(l);
  end
  
%  Compute regularized solution.
  
  f_hat = conj(k_hat) .* d_hat ./ (abs(k_hat).^2 + alpha*lhat);
  f_alpha = real(ifft2(f_hat));  %  Imaginary part is tiny, but is removed.
  f_alpha = f_alpha(1:nfx,1:nfy);
  
%  Display solution.

  mesh(f_alpha)
%  imagesc(f_alpha), colorbar
  title('Reconstructed Object')
  