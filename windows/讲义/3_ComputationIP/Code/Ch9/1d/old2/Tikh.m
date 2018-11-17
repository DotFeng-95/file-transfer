%  Tikh.m
%
%  Compute unconstrained and nonnegatively constrained minimizers of
%  the Tikhonov functional
%    T(f) = ||K*f-d||^2 + alpha*||f||^2
%         = ||A*f - b||^2,
%  where
%    A = [     K       ],  b = [d]
%        [sqrt(alpha)*I]       [0]

  alpha = input(' Regularization parameter alpha = ');
  
  %  Compute unconstrained minimizer.
  
  f_alpha = (K'*K + alpha*eye(n)) \ (K'*d);
  
  %  Compute nonnegatively constrained minimizer.
  
  A = [K; sqrt(alpha)*eye(n)];
  b = [d; zeros(n,1)];
  f_alpha_nonneg = nnls(A,b);
  
  %  Display results.
  
  figure(1)
    plot(x,f_alpha, x,f_true,'--')
    xlabel('x axis')
    title('Reconstruction without Nonnegativity Constraints')
    
  figure(2)
    plot(x,f_alpha_nonneg, x,f_true,'--')
    xlabel('x axis')
    title('Reconstruction with Nonnegativity Constraints')
    