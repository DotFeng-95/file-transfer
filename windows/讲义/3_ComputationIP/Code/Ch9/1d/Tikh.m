%  Tikh.m
%
%  Compute unconstrained minimizers of the Tikhonov functional
%    T(f) = ||K*f-d||^2 + alpha*||f||^2

  alpha = input(' Regularization parameter alpha = ');
  f_alpha = (K'*K + alpha*eye(n)) \ (K'*d);
  
  figure(1)
    plot(x,f_alpha, x,f_true,'--')
    xlabel('x axis')
    title('Reconstruction without Nonnegativity Constraints')
    
  rel_error_tikh = norm(f_alpha-f_true) / norm(f_true)