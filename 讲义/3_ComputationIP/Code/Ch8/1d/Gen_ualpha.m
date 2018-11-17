%  Gen_ualpha.m
%
%  Generate "true" TV reconstruction corresponding to given 
%  regularization parameter values using a primal-dual Newton's method.
%
%  This "true" TV reconstruction is stored in the file u_alpha.mat 
%  and is named u_alpha.

%  Generate data for 1-D test problem

  Setup;
  close all;

%  Set parameters.

  alpha = 1e-4;
  beta = .1;
  grad_tol = eps;  
  PDnewt_iter = 20;
  
%  Generate TV reconstruction using primal-dual method.

  Primal_dual;

  [dum,n_iter] = size(PDconv_history);
  u_alpha = PDconv_history(:,n_iter);
  save u_alpha u_alpha



