  function [u_mat,A0,grad_ls] = co_lsgrad(q_mat,b,z,Cindex)
%
%  Costate computation of gradient of least squares functional,
%      grad_ls = <C[du/dq*e_i],resid>, i=1,...,N.

    [nx,ny] = size(q_mat);
    N = nx * ny;
    [A0] = get_Amat(q_mat);
    uvec = A0 \ b;
    u_mat = reshape(uvec,nx,ny);
    expinvq = exp(-q_mat);

%  Compute residual = C*u - z.

    Cu = uvec(Cindex);  
    resid = Cu - z(:);

%  Compute costate = -adjoint(C) * residual.

    Cstar_resid = zeros(N,1);
    Cstar_resid(Cindex) = resid;
    costate = -A0 \ Cstar_resid;
    costate = reshape(costate,nx,ny);

%  Contribution from "x-operators"

  divTx_costate = [costate(1,:); diff(costate); -costate(nx,:)];
  dudx = [zeros(1,ny); diff(u_mat); zeros(1,ny)];
  dDxdq = 2 ./ (expinvq(1:nx-1,:) + expinvq(2:nx,:)).^2;
  dDxdq = [exp(2*q_mat(1,:)); dDxdq; exp(2*q_mat(nx,:))];
  tx = dDxdq .* divTx_costate .* dudx;
  gradx = (tx(1:nx,:) + tx(2:nx+1,:)) .* expinvq;

%  Contribution from "y-operators"

  divTy_costate = [costate(:,1) diff(costate')' -costate(:,ny)];
  dudy = [2*u_mat(:,1) diff(u_mat')' -2*u_mat(:,ny)];
  dDydq = 2 ./ (expinvq(:,1:ny-1) + expinvq(:,2:ny)).^2;
  dDydq = [exp(2*q_mat(:,1)) dDydq exp(2*q_mat(:,ny))];
  ty = dDydq .* divTy_costate .* dudy;
  grady = (ty(:,1:ny) + ty(:,2:ny+1)) .* expinvq;

%  Gradient is sum.

  grad_ls = gradx + grady;
