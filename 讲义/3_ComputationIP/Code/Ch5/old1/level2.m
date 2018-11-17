  function [c,C] = level2(t)

%  c = level2(t)
%  [c,C] = level2(t)
%  
%  Construct level 2 preconditioning matrix C for a 
%  block Toeplitz-Toeplitz block matrix T.
%  t is the representer for T, i.e., T = bttb(t).
%  C is block circulant with circulant blocks with representer c,
%  i.e., C = bccb(c).

  [nnx,nny] = size(t);
  if (mod(nnx,2)==0 | mod(nny,2)==0)
    fprintf('\n *** Row and column dimensions of input t must be odd.\n');
    return
  end
  
  nx = ceil(nnx/2);  %  nnx = 2*nx - 1
  ny = ceil(nny/2);  %  nny = 2*ny - 1

  %  Level 1
  
  c1 = zeros(nx,nny);
  for j = 1:nny
    c1(:,j) = optimal_circ_toep(t(:,j));
  end
  
  %  Level 2
  
  c = zeros(nx,ny);
  for i = 1:nx
    c(i,:) = optimal_circ_toep(c1(i,:)).';
  end
  
  if nargout == 2
    C = bccb(c);
  end
  
