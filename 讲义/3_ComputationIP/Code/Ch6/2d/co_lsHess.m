  function [Hv] = co_lsHess(v,q0,u0,A0,Cindex)
%
%  Use costate method to evaluate H_ls * v.

  [nx,ny] = size(q0);
  N = nx * ny;
  expinvq = exp(-q0);

%  Compute du/dx, du/dy, dDx/dq, dDy/dq.

  dudx = [zeros(1,ny); diff(u0); zeros(1,ny)];
  dudy = [2*u0(:,1) diff(u0')' -2*u0(:,ny)];

  dDxdq = 2 ./ (expinvq(1:nx-1,:) + expinvq(2:nx,:)).^2;
  dDxdq = [exp(2*q0(1,:)); dDxdq; exp(2*q0(nx,:))];

  dDydq = 2 ./ (expinvq(:,1:ny-1) + expinvq(:,2:ny)).^2;
  dDydq = [exp(2*q0(:,1)) dDydq exp(2*q0(:,ny))];

%  Apply dD/dq to v, and then apply divergence operator to result times
%  grad u.

  ev = expinvq .* v;
  dDxdq_v = dDxdq(2:nx,:) .* (ev(1:nx-1,:) + ev(2:nx,:));
  dDxdq_v = [dDxdq(1,:).*ev(1,:); dDxdq_v; dDxdq(nx+1,:).*ev(nx,:)];
  div_dDxdq_v = diff(dDxdq_v .* dudx);

  dDydq_v = dDydq(:,2:ny) .* (ev(:,1:ny-1) + ev(:,2:ny));
  dDydq_v = [dDydq(:,1).*ev(:,1) dDydq_v dDydq(:,ny+1).*ev(:,ny)];
  div_dDydq_v = diff((dDydq_v .* dudy)')';

  div_dDdq_v = div_dDxdq_v + div_dDydq_v;

%  Apply -A0^{-1} * adjoint(C) * C * A0^{-1}, and then 
%  apply adjoint divergence operator.

  w0 = A0 \ div_dDdq_v(:);
  w = zeros(N,1);
  w(Cindex) = w0(Cindex);
  w = -A0 \ w;
  w = reshape(w,nx,ny);

  divTx_w = [w(1,:); diff(w); -w(nx,:)];
  divTy_w = [w(:,1) diff(w')' -w(:,ny)];

%

  tx = dDxdq .* dudx .* divTx_w;
  Hvx = (tx(1:nx,:) + tx(2:nx+1,:)) .* expinvq;

  ty = dDydq .* dudy .* divTy_w;
  Hvy = (ty(:,1:ny) + ty(:,2:ny+1)) .* expinvq;

  Hv = Hvx + Hvy;
