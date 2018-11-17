%  gen_data.m
%
%  Generate image data d = conv(s,f) + noise, where conv(.,.) denotes
%  discrete convolution. 

  poisson_flag = 1; %%%input(' Enter 1 for Poisson error; else enter 0: ');
  gauss_var = 2; %%%input(' Gaussian error variance = ');
  iseed = 0; %%%input(' Seed for random number generator = ');
  randn('state',iseed);    % Reset random number generator to initial state.
  
  load pupil
  [nx,ny] = size(pupil);
  n = sum(sum(pupil));

  nfx = nx; nfy = ny;
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
  c_f = 1e4;
  f_true = c_f * f_true;
  f_true = f_true';
  f_true = extend(f_true,Nx,Ny);   %  Extend f_true array to size Nx X Ny.


  [Nx,Ny] = size(f_true);
  pupil = extend(pupil,Nx,Ny);   %  Extend pupil array to size Nx X Ny.
  indx = [0:nx, ceil(nx)-1:-1:1]';
  [indx1,indx2] = meshgrid(indx,indx');
  rsq = indx1.^2 + indx2.^2;
  wn = fft2(randn(Nx,Ny));
  
  C = 200;
  K0sq = 1e-4;
  phase_filter = (rsq + K0sq).^(-11/12);
  spectrum = wn .* phase_filter;
  phi = real(ifft2(spectrum));
  phi = phi - mean(phi(:));
  phi = C * phi;
  
%  Generate simulated data. Convolve object with PSF and add error.

  imath = sqrt(-1);
  H = fftshift( pupil.*exp(imath*phi) );
  s = abs(ifft2(H)).^2;
  S = fft2(s);
  S = S / S(1,1);  %  Scale so mean is preserved.
  F_true = fft2(f_true);
  D = S .* F_true;
  d = real(ifft2(D));

  if poisson_flag == 1    %  Poisson error.
    fprintf(' ... adding Poisson error to image data ...\n');
    d = randpoisson(d);
  end
  if gauss_var > 0        %  Gaussian error.
    d = d + gauss_var *  randn(Nx,Ny);
  end
  
%  Fourier transform data.

  D = fft2(d);

%  Plot data.

  figure(1)
    subplot(221)
      imagesc(extract(f_true,nx,ny)), colorbar
      title('Object f')
    subplot(222)
      imagesc(extract(fftshift(s),nx,ny)), colorbar
      title('Point Spread Function  s')
    subplot(223)
      imagesc(extract(d,nx,ny)), colorbar
      title('Noisy, Blurred Image d')
    
