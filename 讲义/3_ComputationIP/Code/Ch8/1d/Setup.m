%  setup.m
%
%  MATLAB code written by Curt Vogel, Dept of Mathematical Sciences,
%  Montana State University, for Chapter 1 of the SIAM Textbook,
%  "Computational Methods for Inverse Problems".
%
%  Set up a discretization of a convolution integral operator K with a 
%  Gaussian kernel. Generate a true solution and convolve it with the 
%  kernel. Then add random error to the resulting data. 
%  Also, compute an eigendecomposition of the discretized operator.

%  Set up parameters.

  n = 80; %%%input(' No. of grid points = ');
  sig = .05; %%%input(' Kernel width sigma = ');
  err_lev = 2; %%%input(' Percent error in data = ');
  
%  Set up grid.

  h = 1/n;
  x = [h/2:h:1-h/2]';
  
%  Compute matrix K corresponding to convolution with Gaussian kernel.

  kernel = (1/sqrt(pi)/sig) * exp(-(x-h/2).^2/sig^2);
  K = toeplitz(kernel)*h;

%  Set up true solution f_true and data d = K*f_true + error.
 
  f_true = .75*(.1<x&x<.25) + .25*(.3<x&x<.32) + (.5<x&x<1).*sin(2*pi*x).^4;
  Kf = K*f_true;
  randn('state',0);
  eta = err_lev/100 * norm(Kf) * randn(n,1)/sqrt(n);
  d = Kf + eta;
  
%  Display the data.

  figure(1)
    plot(x,f_true,'-', x,d,'o',x,Kf,'--')
    xlabel('x axis')
    title('True Solution and Discrete Noisy Data')
    
%  Compute an eigendecomposition of K. K is symmetric, so this is
%  equivalent to an SVD.

  [V,svals] = eig(K);
  [svals,indx] = sort(-diag(svals));  %  Sort -eigenvalues in decreasing order.
  svals = -svals;                     %  +eigs are in increasing order.
  V = V(:,indx);                %  Corresponding eigenvectors.
  
%  Plot the eigenvalues (which are the singular values) and a few of the
%  corresponding eigenvectors.
  
  figure(2)
    semilogy(svals,'o')
    xlabel('index i')
    ylabel('\sigma_i')
    title('Singular Values of K')
  figure(3)
    subplot(221), plot(x,V(:,1)), xlabel('x axis')
    title('Singular Vector v_1')
    subplot(222), plot(x,V(:,4)), xlabel('x axis')
    title('v_4')
    subplot(223), plot(x,V(:,10)), xlabel('x axis')
    title('v_{10}')
    subplot(224), plot(x,V(:,20)), xlabel('x axis')
    title('v_{20}')
    