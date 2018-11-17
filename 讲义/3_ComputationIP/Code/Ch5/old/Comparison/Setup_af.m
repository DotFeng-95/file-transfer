%  setup.m  
%
%  Load Air Force satellite exact solution u_exact, kernel K, 
%  and data z on cell-centered mesh, where we assume the model
%
%      z(x,y) = \int_Omega k(x-x',y-y') u_exact(x',y') dx' dy'
%                         + error.
%
%  Omega is the unit square in R^2.
%
%  Also, extend the kernel by zero (to reduce edge effects), compute its 
%  Fourier transform khat, compute the Fourier transform abs(khat).^2 
%  of and K^*K, and evaluate K^* z.

  load af 
  n_reduce = input(' No. of levels of data reduction = ');

  if n_reduce > 0
    for kk = 1:n_reduce

    %  Sum over 4 neighbors.
  
      n = max(size(z));
      kx = kernel(1:2:n-1,:) + kernel(2:2:n,:);
      kernel = kx(:,1:2:n-1) + kx(:,2:2:n);
      kernel = kx(:,1:2:n-1) + kx(:,2:2:n);
      ux = u_exact(1:2:n-1,:) + u_exact(2:2:n,:);
      u_exact = ux(:,1:2:n-1) + ux(:,2:2:n);
      zx = z(1:2:n-1,:) + z(2:2:n,:);
      z = zx(:,1:2:n-1) + zx(:,2:2:n);
      clear kx ux zx
    end
    n = n/2
  else
    n  = max(size(z));
  end
  
  nx = n;
  ny = n;
  hx = 1 / nx;
  hy = 1 / ny;
  
%  Plot kernel, exact solution, and blurred, noisy data

  figure(1)
  colormap(jet)
  subplot(221), imagesc(kernel), title('Point Spread Function'), colorbar
  subplot(223), imagesc(u_exact), title('True Image'), colorbar
  subplot(224), imagesc(z), title('Noisy, Blurred Data'), colorbar

%  Extend kernel and compute its 2-d Fourier transform. Then use this to 
%  compute K'*z and kstark_hat, the 2-d Fourier transform for K'*K. 

  m2 = 2*n;
  nd2 = n / 2;
  kernele = zeros(m2, m2) ;
  kernele(nd2+1:n+nd2,nd2+1:n+nd2) = kernel ;

  k_hat = fft2(fftshift(kernele)) ;
  Kstar_dat = convolve_2d(conj(k_hat),z); 
  Kstar_dat = Kstar_dat(1:nx,1:ny);  %%%integral_op(z,conj(k_hat));

  dat = z;
  f_true = u_exact;

  clear z u_exact