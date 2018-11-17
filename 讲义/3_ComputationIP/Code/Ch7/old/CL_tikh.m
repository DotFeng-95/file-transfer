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
  alphavec = logspace(-6,1,60);
  n = nx*ny;
  
%  Construct Fourier representer for I.

  lhat = ones(nx,ny);
  
  CL = zeros(size(alphavec));
  GCV = zeros(size(alphavec));
  enorm = zeros(size(alphavec));
  pnorm = zeros(size(alphavec));

  for ii = 1:max(size(alphavec))
    alpha = alphavec(ii);
    A_hat = abs(k_hat).^2 ./ (abs(k_hat).^2 + alpha);
    resid = abs(((A_hat - 1) .* d_hat));
    rms_resid = sum(resid(:).^2) / n;
    trace = sum(A_hat(:));
    Vdenom = (1 - trace/n)^2;
   
    %  Compute CL(alpha) and GCV(alpha).

    CL(ii) = rms_resid + 2*stdev^2/n * trace - stdev^2;
    GCV(ii) = rms_resid / Vdenom;

    %  Compute regularized solution, estimation (solution) error, and
    %  prediction error.
    
    f_hat = conj(k_hat) .* d_hat ./ (abs(k_hat).^2 + alpha);
    enorm(ii) = norm(f_hat(:) - f_true_hat(:))^2 / n^2;
    p_error = k_hat .* (f_hat - f_true_hat);
    pnorm(ii) = norm(p_error(:))^2 / n;

    %  Display regularized reconstruction.

    figure(2)
      f_alpha = real(idft(f_hat));
      imagesc(f_alpha(1:65,1:65)), colorbar
      title('Reconstructed Image')
  end

  [CL_min,CL_indx] = min(CL);
  [GCV_min,GCV_indx] = min(GCV);
  [p_min,p_indx] = min(pnorm);
  [e_min,e_indx] = min(enorm);
  
  %  Display results.
  
  fprintf(' CL(alpha)  is minimized at alpha = %6.4e\n', alphavec(CL_indx));
  fprintf(' GCV(alpha) is minimized at alpha = %6.4e\n', alphavec(GCV_indx));
  fprintf(' P(alpha)   is minimized at alpha = %6.4e\n', alphavec(p_indx));
  fprintf(' E(alpha)   is minimized at alpha = %6.4e\n', alphavec(e_indx));
  
  loglog(alphavec,CL,'-', alphavec,pnorm,'--', alphavec,GCV,'-.', alphavec,100*enorm,'o-')
  xlabel('Regularization Parameter \alpha')
  title('CL(\alpha) is solid line; P(\alpha) is dashed line; GCV(\alpha) is dot-dashed line.') 
  

