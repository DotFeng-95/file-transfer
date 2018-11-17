  function C = level1(t)
  
%  Construct level 1 preconditioning matrix for T = BTTB(t).

  [nnx,nny] = size(t);
  if (mod(nnx,2)==0 | mod(nny,2)==0)
    fprintf('\n *** Row and column dimensions of input t must be odd.\n');
    return
  end
  
  nx = ceil(nnx/2);  %  nnx = 2*nx - 1
  ny = ceil(nny/2);  %  nny = 2*ny - 1
  N = nx*ny;         %  T is N X N
  C = zeros(N,N);
  
  for j = 1:ny
    
    %  Replace each Tj = toeplitz(t(:,j)) with its optimal circulant 
    %  approximation.
    
    [cj,Cj] = optimal_circ_toep(t(:,j));
    xindex = [1:j];
    yindex = [ny-j+1:ny];
    for k = 1:j
      i0 = (xindex(k) - 1)*nx;
      j0 = (yindex(k) - 1)*nx;
      C(i0+1:i0+nx,j0+1:j0+nx) = Cj;
    end
  end
  for j = ny+1:nny
    [cj,Cj] = optimal_circ_toep(t(:,j));
    xindex = [j-ny+1:ny];
    yindex = [1:nny-j+1];
    for k = 1:nny-j+1
      i0 = (xindex(k) - 1)*nx;
      j0 = (yindex(k) - 1)*nx;
      C(i0+1:i0+nx,j0+1:j0+nx) = Cj;
    end
  end
    
