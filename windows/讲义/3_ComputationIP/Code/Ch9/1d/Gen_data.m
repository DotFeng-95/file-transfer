%  Gen_data.m
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
  stdev = 1; %%%input(' Std. deviation of Gaussian error in data = ');
  f_scale = 1000; %%%input(' Scaling parameter for object = ');
  
%  Set up grid.

  h = 1/n;
  x = [h/2:h:1-h/2]';
  
%  Compute matrix K corresponding to convolution with Gaussian kernel.

  kernel = (1/sqrt(pi)/sig) * exp(-(x-h/2).^2/sig^2);
  K = toeplitz(kernel)*h;

%  Set up true solution f_true and data d = K*f_true + error.
 
  a = .75;  b = .85;
  f_true = .75*(.1<x&x<.25) + .25*(.3<x&x<.32) + 500*(a<x&x<b).*(x-a).*(b-x);
  f_true = f_scale * f_true;
  Kf = K*f_true;
  randn('state',0);
  d = randpoisson(Kf) + stdev * randn(n,1);
  
%  Display the data.

  figure(1)
    plot(x,f_true,'-', x,d,'o',x,Kf,'--')
    xlabel('x axis')
    title('True Solution and Discrete Noisy Data')
    
