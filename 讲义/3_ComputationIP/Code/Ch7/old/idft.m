  function xhat = idft(x);

  [m,n] = size(x);
  xhat = ifft2(x) * (sqrt(m)*sqrt(n));
