  function B = extend(A,Nx,Ny)
%  B = extend(A,Nx,Ny)
%
%  Zero-extend array A to size Nx X Ny.

  [nx,ny] = size(A);
  nxd2 = floor(nx/2);
  nyd2 = floor(ny/2);
  B = zeros(Nx,Ny);
  B(nx-nxd2+1:Nx-nxd2,ny-nyd2+1:Ny-nyd2) = A;