%  setup.m
%
%  Generate data for 2-d parameter identification problem.

  nx = input(' nx = ');
  ny = input(' ny = ');
  hx = 1/nx;
  hy = 1/ny;
  xmid = [hx/2:hx:1]';
  ymid = [hy/2:hy:1]';

  [YMID,XMID] = meshdom(ymid,xmid);

  q0 = 0;
  q1 = log(.1);
  qmask = (YMID < 1/2);
  q_exact = q0*ones(size(XMID)).*qmask + q1*ones(size(XMID)).*(1-qmask);
  kappa = exp(q_exact);

  A = get_Amat(q_exact);

  N = nx * ny;
  b = zeros(N,1);
  indxx0 = floor((nx+1)/2);
  indxy0 = floor((ny+1)*3/5);
  b((indxy0-1)*nx+indxx0) = 1;   %  Dirac delta as x = 1/2, y = 3/5.

  uvec = A\b;
  u_mat = reshape(uvec,nx,ny);

%  Set up and apply observation operator C. Take every other grid point 
%  as an observation point.

  indx1 = [2:2:nx]';
  indx2 = [2:2:ny]';
  tmpx = zeros(nx,1);
  tmpx(indx1) = ones(size(indx1));
  tmpy = zeros(ny,1);
  tmpy(indx2) = ones(size(indx2));
  tmp = tmpx * tmpy';
  tmp = tmp(:);
  Cindex = (tmp == 1);

  %  Omit observations in neighborhood of point source.

  x0 = (indxx0 - 1/2) * hx;
  y0 = (indxy0 - 1/2) * hy;
  r0 = .1;
  Cindex = ((XMID - x0).^2 + (YMID - y0).^2 > r0^2);
  Cindex = Cindex(:);

%  Cindex = [1:2:N]';
  z = uvec(Cindex);
  Nz = max(size(z));

%  Add error to data z.

  ratio = input(' Noise/Signal ratio = ');
  z = z + ratio * norm(uvec) * randn(Nz,1) / sqrt(Nz);

%  Plot data

  figure(1)
  subplot(221)
    imagesc(kappa)
    xlabel('y index')
    ylabel('x index')
    title('Intensity plot of kappa(x,y)')
    colormap(jet)
    colorbar

  subplot(222)
    imagesc(u_mat)
    xlabel('y index')
    ylabel('x index')
    title('Intensity plot of u(x,y)')
    colormap(hot)
    colorbar

  subplot(223)
    plot(ymid,u_mat(indxx0,:))
    xlabel('y axis')
    ylabel('u(x,y0)')
    title('Cross section u(x=x0,y)')

%  Compute flow field.

  Dux = [zeros(1,ny); diff(u_mat)/hx; zeros(1,ny)];
  Dux = (Dux(1:nx,:) + Dux(2:nx+1,:)) / 2;
  Duy = [zeros(1,nx); diff(u_mat')/hy; zeros(1,nx)]';
  Duy = (Duy(:,1:ny) + Duy(:,2:ny+1)) / 2;
  flowx = -kappa .* Dux;
  flowy = -kappa .* Duy;
  scale = max(max(sqrt(flowx.^2 + flowy.^2)));

  subplot(224)
    contour(u_mat,15), hold on, quiver(flowy,flowx), hold off
    xlabel('y index')
    ylabel('x index')
    title('Fluid Flow Field')
