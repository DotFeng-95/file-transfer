%  Compare.m

  alpha = input(' Regularization parameter alpha = ');
  resetflag = input(' Enter 1 to reset preconditioner: ');    

%  Set up A11 matrix and compute its Choleski factorization.

 if resetflag
  [nx,ny] = size(dat);
  n = nx*ny;
  p = log2(nx)
  p0 = input(' Coarse grid level p0 = ');
  n0 = 2^p0
  N = n0^2;
  b = Kstar_dat;
%  L = laplacian(nx,ny);
L = speye(nx*ny);
  disp(' ... Setting up A11 ...');
  A11 = zeros(N,N);
  G = zeros(N,N);
  for i = 1:N
      phi_i = zeros(N,1);
      phi_i(i) = 1;
      phi_i = prolongs(reshape(phi_i,n0,n0), p-p0);
      APhi_i = amult(phi_i,k_hat,alpha,L);
      Phistar_APhi_i = restricts(APhi_i,p-p0);
      A11(:,i) = Phistar_APhi_i(:);
      LPhi_i = reshape(L*phi_i(:),nx,ny); 
      Phistar_LPhi_i = restricts(LPhi_i,p-p0);
      G(:,i) = Phistar_LPhi_i(:);
  end
  AR = chol(A11);  %  A11 = AR'*AR.
  LR = chol(L);    %  L = LR'*LR.
  GR = chol(G);    %  G = GR'*GR.
  disp(' ... A11 setup complete ...')
  inv_hat = 1./(abs(k_hat).^2 + alpha);
 end
  
%  Generate highly accurate approximate solution u_alpha. Use this 
%  as the "truth".

  cg_maxiter = 50;
  [u_alpha] = circ_cg(k_hat,b,L,alpha,cg_maxiter, ...
      inv_hat,f_true);

%  Generate Schur complement CG approximation.

  fprintf('      Schur Complement CG\n')
  cg_maxiter = 15;
  [u,schur_enormvec] = schur_cg(k_hat,b,L,alpha,cg_maxiter,AR,LR, ...
      p,p0,f_true,u_alpha);

%  Generate symmetric multiplicative Schwarz PCG approximation.

  fprintf('      SMS PCG\n')
  cg_maxiter = 15;
  [u,sms_enormvec] = sms_pcg(k_hat,b,L,alpha,cg_maxiter,AR,LR,GR, ...
      p,p0,f_true,u_alpha);

%  Generate additive Schwarz PCG approximation.

  fprintf('      Additive Schwarz PCG\n')
  cg_maxiter = 25;
  [u,as_enormvec] = as_pcg(k_hat,b,L,alpha,cg_maxiter,AR,LR,GR, ...
      p,p0,f_true,u_alpha);

%  Generate block circulant extension PCG approximation.

  fprintf('      Block Circulant PCG\n')
  cg_maxiter = 25;
  [uh,circ_enormvec] = circ_cg(k_hat,b,L,alpha,cg_maxiter, ...
      inv_hat,f_true,u_alpha);

%  Generate unpreconditioned CG approximation.

  fprintf('      Unpreconditioned CG\n')
  cg_maxiter = 25;
  [uh,cg_enormvec] = cg(k_hat,b,L,alpha,cg_maxiter,f_true,u_alpha);

%  Display results

 figure(3)
  schur_indx = [1:max(size(schur_enormvec))];
  bce_indx = [1:max(size(circ_enormvec))];
  cgindx = [1:max(size(cg_enormvec))];
  asindx = [1:max(size(as_enormvec))];
  smsindx = [1:max(size(sms_enormvec))];
  semilogy(cgindx,cg_enormvec,'*-', schur_indx,schur_enormvec,'o-', bce_indx,circ_enormvec,'+-', asindx,as_enormvec,'x-', smsindx,sms_enormvec,'v-')
  xlabel('CG iteration')
  title('Soln Error Norm for CG(*), Schur Compl(o), Circ(+), AS(x), SMS(v)')
