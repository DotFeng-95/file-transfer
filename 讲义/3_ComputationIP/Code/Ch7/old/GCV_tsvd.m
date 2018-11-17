%  GCV_tsvd.m
%
%  Implement the method method of generalized cross validation (GCV) for 
%  selecting the regularization parameter in the Tikhonov solution 
%  of K*f=d,
%    f_alpha = DFT^{-1}( sum_{|k_hat|_i^2 > alpha} d_hat/k_hat).
%  Here d denotes image data, K is the blurring operator, and alpha is 
%  a positive regularization parameter. It is assumed that BCCB(PSF) 
%  is the matrix representing K, with lexicographical column ordering 
%  of the unknowns. "hat" denotes Fourier transform.
%
%    A Fourier representation is used, yielding 
%      V(alpha) = sum|(A_hat(alpha)-1)*d_hat|^2/n ...
%                   / [1 - sum(A_hat(alpha))/n]^2,
%  where
%      A_hat(alpha) = abs(k_hat)^2 / (abs(k_hat)^2 + alpha).

  k_hat = fft2(PSF);
  d_hat = fft2(dat);
  f_extend = zeros(nx,ny);
  f_extend(1:nfx,1:nfy) = f_true;
  f_true_hat = fft2(f_extend);
  alphavec = logspace(-8,0,30);
  n = nx*ny;
  
  CL = zeros(size(alphavec));
  V = zeros(size(alphavec));
  enorm = zeros(size(alphavec));
  pnorm = zeros(size(alphavec));

  for ii = 1:max(size(alphavec))
    alpha = alphavec(ii);
    A_hat = (abs(k_hat).^2 > alpha);
    Vnum = ( abs(((A_hat - 1) .* d_hat)).^2 ) / n;
    Vnum = sum(Vnum(:));
    Vdenom = (1 - sum(A_hat(:))/n)^2;
    CL_S2 = 2*stdev^2 * sum(A_hat(:));
    
    %  Compute V(alpha).

    V(ii) = Vnum / Vdenom;
    CL(ii) = Vnum + CL_S2;

    %  Compute regularized solution, estimation (solution) error, and
    %  prediction error.
    
    f_hat = A_hat .* d_hat ./ (k_hat + eps);
    enorm(ii) = norm(f_hat(:) - f_true_hat(:))^2 / n^2;
    p_error = k_hat .* (f_hat - f_true_hat);
    pnorm(ii) = norm(p_error(:))^2 / n;
  end

  [V_min,V_indx] = min(V);
  [CL_min,CL_indx] = min(CL);
  [p_min,p_indx] = min(pnorm);
  [e_min,e_indx] = min(enorm);
  
  %  Display results.
  
  fprintf(' V(alpha) is minimized at alpha  = %6.4e\n', alphavec(V_indx));
  fprintf(' CL(alpha) is minimized at alpha = %6.4e\n', alphavec(CL_indx));
  fprintf(' P(alpha)  is minimized at alpha = %6.4e\n', alphavec(p_indx));
  fprintf(' E(alpha)  is minimized at alpha = %6.4e\n', alphavec(e_indx));
  
  loglog(alphavec,V,'-', alphavec,CL,'--', alphavec,pnorm,'-.')
  xlabel('Regularization Parameter \alpha')
  title('V(\alpha) is solid line; CL(\alpha) is dashed line; P(\alpha) is dot-dashed line.')
  

