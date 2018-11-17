  function xhat = dft(x);

  [m,n] = size(x);
  xhat = fft2(x) / (sqrt(m)*sqrt(n));
