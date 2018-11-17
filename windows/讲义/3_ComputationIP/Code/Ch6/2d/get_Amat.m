  function [A] = get_Amat(q)
%
%  Set up "stiffness" matrix A for CCFD discretization of diffusion operator
%               A(q)u = -div( kappa grad u),
%  where kappa = exp(q). Boundary conditions are
%                 u = 0, at y=0, y=1.
%             du/dx = 0, at x=0, x=1.
%  Use MATLAB's sparse matrix format. A is nx X ny. It is block tridiagonal
%  with ny^2 blocks, each of which is nx X nx. The diagonal blocks are each
%  tridiagonal matrices. The sub- and superdiagonal blocks are each diagonal
%  matrices. This corresponds to lexicographical ordering of the unknowns.

%  Initialization.

  [nx,ny] = size(q);
  N = nx * ny;
  Amain = zeros(N,1);
  Asub = Amain;  Asuper = Amain;
  Asubsub = Amain; Asupersuper = Amain;

%  Compute kappax and kappay, values for kappa at cell interfaces, by
%  harmonic average over adjacent cells.

  kappa_inv = exp(-q);
  kappax = 2 ./ (kappa_inv(1:nx-1,:) + kappa_inv(2:nx,:));
  kappay = 2 ./ (kappa_inv(:,1:ny-1) + kappa_inv(:,2:ny));

%  Values of kappa on boundaries y=0 and y=1.

  ky0 = exp(q(:,1));
  ky1 = exp(q(:,ny));

%  Contribution from kappax.

  for j = 1:ny
    dx = kappax(:,j);
    i0 = (j-1) * nx + 1;
    i1 = j * nx;
    Amain(i0:i1) = [dx(1); dx(1:nx-2)+dx(2:nx-1); dx(nx-1)];
    Asub(i0:i1) = [-dx; 0];
    Asuper(i0:i1) = [0; -dx];
  end

%  Contribution from kappay.

  Amain(1:nx) = Amain(1:nx) + kappay(:,1);
  for j = 1:ny-2
    i0 = j * nx + 1;
    i1 = (j+1) * nx;
    Amain(i0:i1) = Amain(i0:i1) + kappay(:,j) + kappay(:,j+1);
  end
  i0 = (ny-1) * nx + 1;
  Amain(i0:N) = Amain(i0:N) + kappay(:,ny-1);

  for j = 1:ny-1
    i0 = j * nx + 1;
    i1 = (j+1) * nx;
    Asubsub(i0-nx:i1-nx) = -kappay(:,j);
    Asupersuper(i0:i1) = -kappay(:,j);
  end

%  Set up sparse matrix A.

  A = spdiags([Asubsub Asub Amain Asuper Asupersuper], [-nx -1 0 1 nx], N,N);

%  Modify to incorporate Dirichlet BC at top and bottom boundaries.

  N = nx * ny;

  for i = 1:nx
    A(i,i) = A(i,i) + 2*ky0(i);
    A(N-nx+i,N-nx+i) = A(N-nx+i,N-nx+i) + 2*ky1(i);
  end
