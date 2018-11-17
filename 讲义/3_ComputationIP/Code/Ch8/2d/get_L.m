  function [L] = get_L(kappa,BCflag)
%
%  [L] = get_L(kappa,BCflag)
%
%      L*v = -div (kappa grad v), in Omega = [0,1]x[0,1]
%  If BCflag = 'D', boundary conditions are homogeneous Dirichlet,
%      v = 0 on boundary of Omega.
%  If BCflag = 'N', boundary conditions are homogeneous Neumann,
%      grad v . normal = 0 on boundary of Omega.

%  Use finite differences discretization.
%  L = del_x(kappa del+x(u)) + del_y(kappa(del+y(u));
%      del+x(u) = (x_{i+1,j}-x_{i,j})/h
%      del_x(u) = (x_{i-1,j}-x_{i,j})/h

  [nxp1,nyp1] = size(kappa);
  nx = nxp1 - 1;
  ny = nyp1 - 1;
  N = nx*ny;
  
  %  Build Lx
  
  Lxmain = zeros(N,1);
  Lxsub = zeros(N,1);
  Lxsuper = zeros(N,1);
  for j = 1:ny
    d = kappa(:,j+1);
    i0 = (j-1)*nx + 1;
    i1 = j*nx;
    Lxmain(i0:i1) = d(1:nx)+d(2:nx+1);
    if BCflag == 'N'
      Lxmain(i0) = d(2);
      Lxmain(i1) = d(nx);
    end
    Lxsub(i0:i1) = [-d(2:nx); 0];
    Lxsuper(i0:i1) = [0; -d(2:nx)];
  end

  Lx = spdiags([Lxsub Lxmain Lxsuper], [-1 0 1], N,N);
  
  %  Build Ly

  Lymain = zeros(N,1);
  Lysub = zeros(N,1);
  Lysuper = zeros(N,1);
  for j = 1:ny
    i0 = (j-1)*nx + 1;
    i1 = j*nx;
    Djm1 = kappa(2:nx+1,j);
    Dj = kappa(2:nx+1,j+1);
    if BCflag == 'N'
      if j == 1
        Djm1 = zeros(nx,1);
      elseif j == ny
	Dj = zeros(nx,1);
      end
    end
    Lymain(i0:i1) = Djm1 + Dj;
    Lysub(i0:i1) = -Dj;
    Lysuper(i1+1:i1+nx) = -Dj;
  end
  Lysuper = Lysuper(1:N);

  Ly = spdiags([Lysub Lymain Lysuper], [-nx 0 nx], N,N);

  L = Lx + Ly;