  function kappa = get_kappa(u_mat,beta)
%
%  Compute diffusion coefficient
%    kappa(u) = 1 / sqrt((du/dx)^2 + beta^2)

  [nx,ny] = size(u_mat);
  hx = 1/nx;
  hy = 1/ny;
  Dux = [zeros(1,ny); diff(u_mat)/hx; zeros(1,ny)];
  Dux = (Dux(1:nx,:) + Dux(2:nx+1,:)) / 2;
  Duy = [zeros(1,nx); diff(u_mat')/hy; zeros(1,nx)]';
  Duy = (Duy(:,1:ny) + Duy(:,2:ny+1)) / 2;
  kappa = 1 ./ sqrt(Dux.^2 + Duy.^2 + beta^2);
