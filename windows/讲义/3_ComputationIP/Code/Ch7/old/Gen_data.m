%  Gen_data.m
%
%  Generate atmospheric image data.

%  Generate object f_true.

  nfx = 65; %%%input(' Object size nx = ');
  nfy = nfx;
  h = 2/(nfx-1);
  x = [-1:h:1]';
  [X,Y] = meshgrid(x);
  R = sqrt(X.^2 + Y.^2);
  f1 = (R<.25);
  f2 = (-.2<X-Y) & (X-Y<.2) & (-1.4<X+Y) & (X+Y<1.4);
  f3 = (-.2<X+Y) & (X+Y<.2) & (-1.4<X-Y) & (X-Y<1.4);
  panel = f2.*(1-f1) + f3.*(1-f1);
  f4 = (-.25<X) & (X<.25) & (0<Y) & (Y<.5);
  f5 = (sqrt(X.^2+(Y-.5).^2) < .25);
  body = .5*f1 + (f4.*(1-f1) + f5.*(1-f4));
  body_support = (body>0);
  f_true = .75*panel + body.*(1-panel);
  f_true = f_true';
  c_f = 1e4;
  f_true = c_f * f_true;
  
%  Generate pupil function.

  nx = 128;
  ny = nx;
  icen = nx/2 + 1;  % Pupil is centered at icen. 
  ix = [1:nx]' - icen;
  iy = [1:ny]' - icen;
  [iX,iY] = meshgrid(ix,iy);
  R = nx/2;
  r = sqrt(iX.^2 + iY.^2) / R;
  pupil = (.1 < r & r < .5);

%  Generate phase phi. This is a realization of stochastic process 
%  with Von Karman statistics. See p. 60 of Roggemann and Welsh's
%  "Imaging Through Turbulence", CRC Press, 1996.

  L0m2 = 1e-4 / R;
  c_phi = 3e2;
  indx = [0:nx/2, ceil(nx/2)-1:-1:1]';
  [indx1,indx2] = meshgrid(indx,indx');
  randn('state',0);   %  Reset random number generator to initial state.
  phi = real(ifft2( (sqrt(indx1.^2+indx2.^2).^2 + L0m2).^(-11/12) .* ...
      fft2(randn(nx,ny)) ));
  phi = phi - mean(phi(:));
  phi = c_phi * phi .* pupil;
  
%  Generate point spread function, PSF = |F^{-1}(p*exp(i*phi))|^2,
%  where p is the pupil function, phi is the phase, and F denotes 2-D
%  Fourier transform.

  fprintf(' ... generating image data ...\n');
  H = fftshift( pupil.*exp(sqrt(-1)*phi));
  PSF = abs(ifft2(H)).^2;

  %  Rescale PSF so the convolution operator is mean preserving.
  %  OTF stands for "Optical Transfer Function", the Fourier transform
  %  of the PSF. OTF(1,1) is the constant component of the transform.

  OTF = dft(PSF);
  PSF = PSF / OTF(1,1);
  
%  Generate simulated image data according to the model
%      data = convolution(PSF,f_true) + noise. 

  dat = convolve_2d(dft(PSF), f_true);  %  Convolve PSF with f_true.
  error_percnt = input(' Percent data error = ');
  stdev = .01 * error_percnt * norm(dat(:)) / sqrt(nx*ny);
  dat = dat + stdev*randn(nx,ny);         %  Add Gaussian noise to data.
  
%  Display data.

  figure(1)
      imagesc(phi), colorbar
      title('Phase Function')
  figure(2)
      mesh(fftshift(PSF))
      title('Point Spread Function')
  figure(3)
      imagesc(fftshift(log(abs(fft2(PSF))+sqrt(eps)))), colorbar
      title('Log_{10} of Power Spectrum of PSF')
  figure(4)
      imagesc(f_true), colorbar
      title('Object (True Image)')
  figure(5)
      imagesc(dat(1:nfx,1:nfy)), colorbar
      title('Blurred, Noisy Image')

%  Clear nonessential variables from work space.

  clear error_percnt H indx L0m2 pupil r R X Y body body_support 
  clear c_f c_phi f1 f2 f3 f4 f5 h iX iY icen indx1 indx2 ix iy x panel
