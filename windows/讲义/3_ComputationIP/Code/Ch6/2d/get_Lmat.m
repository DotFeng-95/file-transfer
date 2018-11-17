  function [L] = get_Lmat(kappa)
%
%  Set up "stiffness" matrix L for CCFD discretization of diffusion operator
%               L(kappa)u = -div( kappa(x) grad u)
%  using MATLAB's sparse matrix format. L is nx X ny. It is block tridiagonal
%  with ny^2 blocks, each of which is nx X nx. The diagonal blocks are each
%  tridiagonal matrices. The sub- and superdiagonal blocks are each diagonal
%  matrices. This corresponds to lexicographical ordering of the unknowns.

%  Initialization.

  [nx,ny] = size(kappa);
  N = nx * ny;
  main = zeros(N,1);
  sub = main;  super = main;
  subsub = main; supersuper = main;

%  Contribution from kappa across x-interfaces.

  kappax = 2./ (1./kappa(1:nx-1,:) + 1./kappa(2:nx,:));  % Harmonic avrg
  for j = 1:ny
    dx = kappax(:,j);
    i0 = (j-1) * nx + 1;
    i1 = j * nx;
    main(i0:i1) = [dx(1); dx(1:nx-2)+dx(2:nx-1); dx(nx-1)];
    sub(i0:i1) = [-dx; 0];
    super(i0:i1) = [0; -dx];
  end
  clear kappax

%  Contribution from kappa across y-interfaces.

  kappay = 2./ (1./kappa(:,1:ny-1) + 1./kappa(:,2:ny));  % Harmonic avrg
  main(1:nx) = main(1:nx) + kappay(:,1);
  for j = 1:ny-2
    i0 = j * nx + 1;
    i1 = (j+1) * nx;
    main(i0:i1) = main(i0:i1) + kappay(:,j) + kappay(:,j+1);
  end
  i0 = (ny-1) * nx + 1;
  main(i0:N) = main(i0:N) + kappay(:,ny-1);

  for j = 1:ny-1
    i0 = j * nx + 1;
    i1 = (j+1) * nx;
    subsub(i0-nx:i1-nx) = -kappay(:,j);
    supersuper(i0:i1) = -kappay(:,j);
  end

%  Set up sparse matrix L.

  L = spdiags([subsub sub main super supersuper], [-nx -1 0 1 nx], N,N);

