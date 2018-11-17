  function F = fourier(n)
  
%  Generate n X n Fourier matrix F.
%  F(i,j) = omega^((i-1)*(j-1)) / sqrt(n).


  imath = sqrt(-1);
  omega = exp(-imath*2*pi/n);
  ivec = [0:1:n-1]';
  F = ones(n,n);
  for j = 2:n
    F(:,j) = (omega^(j-1)).^ivec;
  end
  F = F / sqrt(n);
