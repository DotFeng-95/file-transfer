  function [D,C] = level1_bttb(t)

%  D = level1_bttb(t)
%  [D,C] = level1_bttb(t)
%
%  Construct level 1 block circulant approximation C to T = BTTB(t).
%  D is a 3-D array containing the blocks D_k for which 
%     C = tensor(I,F)' * P' * diag(D_k) * P * tensor(I,F).
%  Here P denotes the permutation matrix corresponding to the switch
%  from column lexicographical ordering to row lexicographical
%  ordering; F denotes the nx X nx Fourier matrix; I denotes the 
%  ny X ny identity matrix; and tensor(.,.) denotes tensor product.
%  If only 1 output argument is specified, then only D is returned.

  [nnx,nny] = size(t);
  if (mod(nnx,2)==0 | mod(nny,2)==0)
    fprintf('\n *** Row and column dimensions of input t must be odd.\n');
    return
  end
  
  nx = ceil(nnx/2);  %  nnx = 2*nx - 1
  ny = ceil(nny/2);  %  nny = 2*ny - 1
  D = zeros(nx,ny,ny);  
  if nargout == 2
    N = nx*ny;       %  C is N X N
    C = zeros(N,N);
  end
  
  for k = 1:ny
    
    %  Compute optimal circulant approximation Ck 
    %  to Tk = toeplitz(t(:,k)).
    
    c_k = optimal_circ_toep(t(:,k));
    Lambda_k = fft(c_k);
    xindex = [1:k];
    yindex = [ny-k+1:ny];
    if nargout == 2
      C_k = circulant(c_k);
    end

    for l = 1:k
      i = xindex(l);
      j = yindex(l);
      D(:,i,j) = Lambda_k;
      if nargout == 2
        i0 = (i - 1)*nx;
        j0 = (j - 1)*nx;
        C(i0+1:i0+nx,j0+1:j0+nx) = C_k;
      end
    end
  end
    
  for k = ny+1:nny
    c_k = optimal_circ_toep(t(:,k));
    Lambda_k = fft(c_k);
    xindex = [k-ny+1:ny];
    yindex = [1:nny-k+1];
    if nargout == 2
      C_k = circulant(c_k);
    end
    for l = 1:nny-k+1
      i = xindex(l);
      j = yindex(l);
      D(:,i,j) = Lambda_k;
      if nargout == 2
        i0 = (i - 1)*nx;
        j0 = (j - 1)*nx;
        C(i0+1:i0+nx,j0+1:j0+nx) = C_k;
      end
    end
  end
