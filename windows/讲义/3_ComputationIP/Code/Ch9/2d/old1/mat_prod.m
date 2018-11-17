  function g = mat_prod(S,f)

%  g = mat_prod(S,f)

  g = real(ifft2(S.*fft2(f)));
  