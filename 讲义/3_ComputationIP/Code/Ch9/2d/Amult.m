  function[Ax] = Amult(A,x);
  
%  Compute (S*S + gam*I)x where S = [S1 ... Sk]' and 
%  Ti*x = real(ifft2(S{i}.*fft2(x))). 

  global TOTAL_FFTS 
  TOTAL_FFTS = TOTAL_FFTS + 2;

  S = A.S;
  gam = A.gam;
  n_frames = A.n_frames;
  
  Ax = zeros(size(x));
  Fourier_x = fft2(x);
  
  for i = 1:n_frames
    Ax = Ax + (abs(S{i}).^2).*Fourier_x;
  end
  
  Ax = (1/n_frames^2)*real(ifft2(Ax)) + gam*x;
