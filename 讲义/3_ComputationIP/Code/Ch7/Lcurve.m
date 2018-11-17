%  Lcurve.m
%
%  Implement the L-curve methods for selecting the regularization 
%  parameter in the Tikhonov solution of K*f=d,
%      f_alpha = (K^*K+alpha*I)^{-1} * K'*d.
%  Here d denotes image data, K is the blurring operator, and alpha is 
%  a positive regularization parameter. It is assumed that BCCB(PSF) 
%  is the matrix representing K, with lexicographical column ordering 
%  of the unknowns. 
%
%  The L-curve is the parameterized curve (X(alpha),Y(alpha)), where 
%      Y(alpha) = log ||f_alpha||,
%      X(alpha) = log ||r_alpha||,
%      r_alpha = K*f_alpha - d is the residual.
%  The L-curve method is to maximize the curvature,
%      kappa(alpha) = (X''Y' - X'Y'') / ((X')^2 + (Y')^2)^(3/2).


  k_hat = dft(PSF);
  khatsq = abs(k_hat).^2;
  d_hat = dft(dat);
  dhatsq = abs(d_hat).^2;
  f_extend = zeros(nx,ny);
  f_extend(1:nfx,1:nfy) = f_true;
  f_true_hat = dft(f_extend);
  alphavec = logspace(-6,1,60);  %  Range of alpha is 10^{-6} to 10^1.
  n = nx*ny;
  
  R = zeros(size(alphavec));
  S = zeros(size(alphavec));
  kappa = zeros(size(alphavec));
  
  for i = 1:max(size(alphavec))
    alpha = alphavec(i);
    
    %  Compute regularized solution, estimation (solution) error, and
    %  components of the L-curve.
    
    f_hat = conj(k_hat) .* d_hat ./ (abs(k_hat).^2 + alpha);
    enorm(i) = norm(f_hat(:) - f_true_hat(:))^2 / n^2;
    A_hat = abs(k_hat).^2 ./ (abs(k_hat).^2 + alpha);
    Ri = norm( (1 - A_hat) .* d_hat,'fro')^2;
    Si = norm(f_hat(:))^2;
    Sprime = khatsq .* dhatsq ./ (khatsq + alpha).^3;
    Sprime = -2 * sum(Sprime(:));
    kappa(i) = -((Ri*Si)*(alpha*Ri + alpha^2*Si) + Ri^2*Si^2/Sprime) / ...
      (Ri^2 + alpha^2 * Si^2)^(3/2);
    R(i) = Ri;
    S(i) = Si;
  end
  
  [e_min,e_indx] = min(enorm);
  [k_max,k_indx] = max(kappa);
  
%  Display results.

  fprintf(' Estimation  error is minimized at alpha = %6.4e\n', ...
    alphavec(e_indx(1)));
  fprintf(' L-curve curvature is minimized at alpha = %6.4e\n', ...
    alphavec(k_indx(1)));

  figure(3)
    loglog(R,S, R(k_indx(1)),S(k_indx(1)),'o', R(e_indx(1)),S(e_indx(1)),'*')
    xlabel('X(\alpha) = ||K f_{\alpha} - d||^2')
    ylabel('Y(\alpha) = ||f_{\alpha}||^2')
    title('L-curve')
    
  figure(4)
    semilogx(alphavec,kappa, alphavec(k_indx(1)),kappa(k_indx(1)),'o',...
	alphavec(e_indx(1)),kappa(e_indx(1)),'*')
    xlabel('Regularization Parameter \alpha')
    ylabel('\kappa(\alpha)')
    title('Curvature of the L-curve')
    