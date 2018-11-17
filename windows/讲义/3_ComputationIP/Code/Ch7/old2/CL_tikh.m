%  CL_tikh.m
%
%  Implement the CL method for select the regularization parameter
%  in the Tikhonov solution of K*f=d,
%    f_alpha = (K^*K+alpha*I)^{-1} * K'*d.
%  Here d denotes image data, K is the blurring operator, and alpha is 
%  a positive regularization parameter. It is assumed that BCCB(PSF) 
%  is the matrix representing K, with lexicographical column ordering 
%  of the unknowns. 
%  A Fourier representation is used, yielding 
%      CL(alpha) = sum [(A_hat(alpha)-1)*d_hat]^2 + 2*sigma^2/n 
%                    * sum(A_hat(alpha)),
%  where
%      A_hat(alpha) = abs(k_hat)^2 / (abs(k_hat)^2 + alpha).

  k_hat = dft(PSF);
  d_hat = dft(dat);
  f_extend = zeros(nx,ny);
  f_extend(1:nfx,1:nfy) = f_true;
  f_true_hat = dft(f_extend);
  alphavec = logspace(-10,-2,30);
  n = nx*ny;
  
%  Construct Fourier representer for I.

  lhat = ones(nx,ny);
  
  CL = zeros(size(alphavec));
  enorm = zeros(size(alphavec));
  pnorm = zeros(size(alphavec));

  for ii = 1:max(size(alphavec))
    alpha = alphavec(ii);
    A_hat = abs(k_hat).^2 ./ (abs(k_hat).^2 + alpha);
    S1 = abs(((A_hat - 1) .* d_hat)).^2 / n;
    S2 = 2*stdev^2/n * A_hat;
    
    %  Compute CL(alpha).

    CL(ii) = sum(S1(:) + S2(:)) - stdev^2;

    %  Compute regularized solution, estimation (solution) error, and
    %  prediction error.
    
    f_hat = conj(k_hat) .* d_hat ./ (abs(k_hat).^2 + alpha);
    enorm(ii) = norm(f_hat(:) - f_true_hat(:))^2 / n^2;
    p_error = k_hat .* (f_hat - f_true_hat);
    pnorm(ii) = norm(p_error(:))^2 / n;
  end

  [CL_min,CL_indx] = min(CL);
  [p_min,p_indx] = min(pnorm);
  [e_min,e_indx] = min(enorm);
  
  %  Display results.
  
  fprintf(' CL(alpha) is minimized at alpha = %6.4e\n', alphavec(CL_indx));
  fprintf(' P(alpha)  is minimized at alpha = %6.4e\n', alphavec(p_indx));
  fprintf(' E(alpha)  is minimized at alpha = %6.4e\n', alphavec(e_indx));
  
  loglog(alphavec,CL,'-.', alphavec,pnorm,'--', alphavec,100*enorm,'-')
  xlabel('Regularization Parameter \alpha')
  title('CL(\alpha) is dot-dashed line; P(\alpha) is dashed line; E(\alpha) is solid line.')
  

