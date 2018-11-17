  function [kappa] = get_kappa(u,betaTV)
%
%  [kappa] = get_kappa(u,betaTV)
%
%  Build TV diffusion coefficient,
%
%      kappa = 1 / sqrt(u_x^2 + u_y^2 + beta^2).
%
%  Use finite differences discretization.
%
%      kappa = 1/((del+x(u))^2 + (del+y(u))^2 + beta^2)^{1/2}
%      del+x(u) = (x_{i+1,j}-x_{i,j})/h
%      del_x(u) = (x_{i-1,j}-x_{i,j})/h

  [nx,ny] = size(u);
  hx = 1/(nx-1);
  hy = 1/(ny-1);
  
hx = 1; hy = 1;
   
  u_ext = zeros(nx+2,ny+2);
  u_ext(2:nx+1,2:ny+1) = u;
  uy = diff(u_ext) / hy;
  ux = (diff(u_ext')') / hx;
  kappa = 1./sqrt(ux(1:nx+1,:).^2 + uy(:,1:ny+1).^2 + betaTV^2);
