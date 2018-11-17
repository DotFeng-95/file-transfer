%%%%%%%%%%%%%%%%%
  nnx = nfx;
  nny = nnx;
  t = extract(fftshift(PSF),2*nnx-1,2*nny-1);
  nargout = 2;
%%%%%%%%%%%%%%%%%  

  alpha = input(' alpha = ');
  [D,C] = level1_bttb(t);  
  M = C'*C + alpha * eye(n);
  
  x_array = randn(nnx,nny);
  y = reshape(M*x_array(:),nnx,nny);
  
  M1 = zeros(nnx,nny,nny);
  for i = 1:nnx
      Di = zeros(nny,nny);
      Di(:,:) = D(i,:,:);
      M1(i,:,:) = Di'*Di + alpha*eye(nnx);
  end
    
  y_array = c1_mult(M1,x_array);
  y_diff = norm(y-y_array,'fro') / norm(y,'fro')
