%  Filter_fn.m
%
%  Plot filter functions for TSVD and Tikhonov regularization.

  sigma_sq = logspace(-5,1,100);
  alpha = 1e-2; %%%input(' Regularization parameter alpha = ');
  w_TSVD = (sigma_sq > alpha);
  w_Tikh = sigma_sq ./ (sigma_sq + alpha);
  
  semilogx(sigma_sq,w_TSVD,'--', sigma_sq,w_Tikh,'-')
  title('Filter Functions for TSVD and Tikhonov Regularization')
  xlabel('s^2')
  ylabel('w_{\alpha}(s^2)')
  ylim([0 1.2])
  