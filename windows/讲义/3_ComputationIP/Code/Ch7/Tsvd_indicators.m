%  Tsvd_indicators.m
%
%  Implement the UPRE and GCV methods to selecting the regularization  
%  parameter in the TSVD solution of K*f=d,
%      f_alpha = DFT^{-1}( sum_{|k_hat|_i^2 > alpha} d_hat/k_hat).
%  Here d denotes image data, K is the blurring operator, and alpha is 
%  a positive regularization parameter. It is assumed that BCCB(PSF) 
%  is the matrix representing K, with lexicographical column ordering 
%  of the unknowns. 
%    A Fourier representation is used, yielding for the UPRE
%      U(alpha) = sum [(A_hat(alpha)-1)*d_hat]^2 + 2*sigma^2/n 
%                    * sum(A_hat(alpha)),
%  where
%      A_hat(alpha) = 1, if |k_hat|_i^2 > alpha,
%                     0, otherwise.

  k_hat = dft(PSF);
  d_hat = dft(dat);
  f_extend = zeros(nx,ny);
  f_extend(1:nfx,1:nfy) = f_true;
  f_true_hat = dft(f_extend);
  alphavec = logspace(-6,1,60);
  n = nx*ny;
  
%  Construct Fourier representer for I.

  lhat = ones(nx,ny);
  
  U = zeros(size(alphavec));
  GCV = zeros(size(alphavec));
  enorm = zeros(size(alphavec));
  pnorm = zeros(size(alphavec));

  for ii = 1:max(size(alphavec))
    alpha = alphavec(ii);
    A_hat = (abs(k_hat).^2 > alpha);
    resid = abs(((A_hat - 1) .* d_hat));
    rms_resid = sum(resid(:).^2) / n;
    trace = sum(A_hat(:));
    Vdenom = (1 - trace/n)^2;
    
    %  Compute U(alpha) and GCV(alpha).

    U(ii) = rms_resid + 2*stdev^2/n * trace - stdev^2;
    GCV(ii) = rms_resid / Vdenom;

    %  Compute regularized solution, estimation (solution) error, and
    %  prediction error.
    
    f_hat = A_hat .* d_hat ./ (k_hat + eps);
    enorm(ii) = norm(f_hat(:) - f_true_hat(:))^2 / n^2;
    p_error = k_hat .* (f_hat - f_true_hat);
    pnorm(ii) = norm(p_error(:))^2 / n;
  end

  [U_min,U_indx] = min(U);
  [GCV_min,GCV_indx] = min(GCV);
  [p_min,p_indx] = min(pnorm);
  [e_min,e_indx] = min(enorm);
  
  %  Display results.
  
  fprintf(' U(alpha) is minimized at alpha = %6.4e\n', alphavec(U_indx));
  fprintf(' GCV(alpha) is minimized at alpha = %6.4e\n', alphavec(GCV_indx));
  fprintf(' P(alpha)  is minimized at alpha = %6.4e\n', alphavec(p_indx));
  fprintf(' E(alpha)  is minimized at alpha = %6.4e\n', alphavec(e_indx));
  
  loglog(alphavec,U,'-', alphavec,pnorm,'--', alphavec,GCV,'-.', alphavec,100*enorm,'o-')
  xlabel('Regularization Parameter \alpha')
  title('U(\alpha) is solid line; P(\alpha) is dashed line; GCV(\alpha) is dot-dashed line.') 
