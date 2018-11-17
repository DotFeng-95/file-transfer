%  Rand_trace.m
%
%  Compute true and estimated values for the trace of the influence
%  matrix for Tikhonov regularization,
%      A(alpha) = K * (K^*K+alpha*I)^{-1} * K'.
%  Here K is a blurring operator, and alpha is the regularization
%  parameter. K has a matrix representation BCCB(PSF), i.e., it is block 
%  circulant with circulant blocks. Then the exact trace is
%      trace A(alpha) = sum(A_hat(alpha)),
%  where
%      A_hat(alpha) = abs(k_hat)^2 / (abs(k_hat)^2 + alpha),
%  and k_hat is the discrete Fourier transform of the PSF.
%  Use Hutchinson's trace estimator,
%      t(alpha) = u' * A(alpha) * u,
%  where u is realization of a random vector whose entries are
%  independent, have zero mean and unit variance, and take on values
%  +1 or -1.

  fprintf('\n  Be sure to run Gen_data first!\n');
  k_hat = dft(PSF);
  alphavec = logspace(-6,1,60);  %  Range of alpha is 10^{-6} to 10^1.
  n = nx*ny;
  n_rand = input(' No. of random vectors used for each trace estimate = ');
  
  true_trace = zeros(size(alphavec));
  rand_trace = zeros(size(alphavec));

  for ii = 1:max(size(alphavec))
    alpha = alphavec(ii);
    A_hat = abs(k_hat).^2 ./ (abs(k_hat).^2 + alpha);
    true_trace(ii) = sum(A_hat(:));
    
  %  Generate random vectors u with entries +1 or -1.
  
    rand('state',1);     %  Initialize random number generator.
    tbar = 0;
    for k = 1:n_rand
      v = rand(nx,ny);   %  v ~ Uniform(0,1)
      u = (v >= .5) - (v < .5);
      Au = real(ifft2(A_hat .* fft2(u)));
      tbar = tbar + u(:)' * Au(:);
    end
    rand_trace(ii) = tbar / n_rand;
  end

%  Display results.

  figure(1)
    loglog(alphavec,true_trace,'-', alphavec,rand_trace,'--')
    xlabel('Regularization Parameter \alpha')
    title('True trace is solid line; randomized estimate is dashed line.') 
  figure(2)
    relative_error = abs( rand_trace - true_trace )./ true_trace;
    loglog(alphavec,relative_error)
    xlabel('Regularization Parameter \alpha')
    title('Relative error in randomized trace estimate.')

