  function [dudq] = deriv(q_mat,b_vec)
%
%  Compute derivative matrix for mapping
%    q |-> u(q) = A(q)^{-1}b,
%  where 
%    A(q)u = -d/dx(Dx(q) du/dx) - d/dy(Dy(q) du/dy),
%  and Dx(q), Dy(q) are diagonal with entries ...
%  The derivative has the form
%    du/dq = -A(q)^{-1} dA/dq u
%          = +A(q)^{-1} [d/dx dDx/dq dudx + d/dy dDy/dq dudy]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [nx,ny] = size(q_mat);
  N = nx*ny;
  A = get_Amat(q_mat);
  u_vec = A \ b_vec;
  u_mat = reshape(u_vec,nx,ny);

  %  Compute components of gradient.

  dudx = [zeros(1,ny); diff(u_mat); zeros(1,ny)];
  dudy = [2*u_mat(:,1) diff(u_mat')' -2*u_mat(:,ny)];

  %  Evaluate diffusion coefficient. Then take harmonic average 
  %  at x-interfaces to obtain Dx and along y-interfaces to obtain Dy.
  %  First and last columns of Dy correspond to homogeneous Dirichlet
  %  BC along (x,y=0) and (x,y=1), 0<x<1.

  expmq = exp(-q_mat);

  %  Compute dDx/dq * dudx on interior. 

  denomx = 1 ./ (expmq(1:nx-1,:) + expmq(2:nx,:)).^2;
  divx = zeros(N,N);
  for i = 2:nx-1
    dDxdq_dudx_i = 2 * dudx(i,:) .* expmq(i,:) .* denomx(i-1,:);
    dDxdq_dudx_ip1 = 2 * dudx(i+1,:) .* expmq(i,:) .* denomx(i,:);
    divxi_im1 = dDxdq_dudx_i;
    divxi_i = dDxdq_dudx_ip1 - dDxdq_dudx_i;
    divxi_ip1 = -dDxdq_dudx_ip1;
    for j = 1:ny
      indx1 = (j-1)*nx + i;
      divx(indx1-1,indx1) = divxi_im1(j);
      divx(indx1,indx1) = divxi_i(j);
      divx(indx1+1,indx1) = divxi_ip1(j);
    end
  end

  %  Contributions from BC du/dx = 0 at (x=0,y) and (x=1,y), 
  %  0 < y < 1.

  dDxdq_dudx = zeros(nx+1,ny);
  dDxdq_dudx_2 = 2 * dudx(2,:) .* expmq(1,:) .* denomx(1,:);
  divxi_1 = dDxdq_dudx_2;
  divxi_2 = -dDxdq_dudx_2;
  for j = 1:ny
    indx1 = (j-1)*nx + 1;
    divx(indx1,indx1) = divxi_1(j);
    divx(indx1+1,indx1) = divxi_2(j);
  end

  dDxdq_dudx_nx = 2 * dudx(nx,:) .* expmq(nx,:) .* denomx(nx-1,:);
  divxi_nxm1 = dDxdq_dudx_nx;
  divxi_nx = -dDxdq_dudx_nx;
  for j = 1:ny
    indx1 = (j-1)*nx + nx;
    divx(indx1-1,indx1) = divxi_nxm1(j);
    divx(indx1,indx1) = divxi_nx(j);
  end

  %  Compute dDy/dq * dudy on interior. 

  denomy = 1 ./ (expmq(:,1:ny-1) + expmq(:,2:ny)).^2;
  divy = zeros(N,N);
  for i = 2:ny-1
    dDydq_dudy_i = 2 * dudy(:,i) .* expmq(:,i) .* denomy(:,i-1);
    dDydq_dudy_ipq = 2 * dudy(:,i+1) .* expmq(:,i) .* denomy(:,i);
    divyi_im1 = dDydq_dudy_i;
    divyi_i = dDydq_dudy_ipq - dDydq_dudy_i;
    divyi_ip1 = -dDydq_dudy_ipq;
    for j = 1:nx
      indx1 = (i-1)*nx + j;
      divy(indx1-nx,indx1) = divyi_im1(j);
      divy(indx1,indx1) = divyi_i(j);
      divy(indx1+nx,indx1) = divyi_ip1(j);
    end
  end

  %  Contributions from homogeneous Dirichlet BC along (x,y=0) 
  %  and (x,y=1), 0 < x < 1. 

  dDydq_dudy_1 = dudy(:,1) .* exp(q_mat(:,1));
  dDydq_dudy_2 = 2 * dudy(:,2) .* expmq(:,1) .* denomy(:,1);
  divyi_1 = dDydq_dudy_2 - dDydq_dudy_1;
  divyi_2 = -dDydq_dudy_2;
    for j = 1:nx
      indx1 = j;
      divy(indx1,indx1) = divyi_1(j);
      divy(indx1+nx,indx1) = divyi_2(j);
    end

  dDydq_dudy_ny = 2 * dudy(:,ny) .* expmq(:,ny) .* denomy(:,ny-1);
  dDydq_dudy_nyp1 = dudy(:,ny+1) .* exp(q_mat(:,ny));
  divyi_nym1 = dDydq_dudy_ny;
  divyi_ny = dDydq_dudy_nyp1 - dDydq_dudy_ny;
    for j = 1:nx
      indx1 = (ny-1)*nx + j;
      divy(indx1-nx,indx1) = divyi_nym1(j);
      divy(indx1,indx1) = divyi_ny(j);
    end

  %  

  dudq = A \ (divx + divy);
