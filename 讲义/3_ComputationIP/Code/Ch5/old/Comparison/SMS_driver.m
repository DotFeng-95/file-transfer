%  SMS_driver.m
%
%  Use symmetric multiplicative Schwarz PCG iteration to solve Ah*u=b, 
%  where
%      Ah = K'*K + alpha*L
%      b = K'*dat

  alpha = input(' Regularization parameter alpha = ');
  cg_maxiter = input(' Max no. CG iterations = ');
  resetflag = input(' Reset preconditioner? (0 for No, 1 for Yes): ');
  if resetflag == 1
    [nx,ny] = size(dat);
    n = nx*ny;
    hx = 1/nx;
    hy = 1/ny;
    p = log2(nx)
    p0 = input(' Level p0 for coarse subspace = ');
    n0 = 2^p0
    N = n0^2;
    b = Kstar_dat;

    %  Set up regularization operator.
    
    L = speye(nx*ny); %%%laplacian(nx,ny);
    
    %  Set up matrices A11, G and compute their Choleski factorizations.
    
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
  end
  
%  CG iteration.

%  u = zeros(nx,ny);
  Phistar_b = restricts(b,p-p0);
  u = prolongs(reshape(AR\(AR'\Phistar_b(:)),n0,n0),p-p0);
  g = amult(u,k_hat,alpha,L) - b;
  residnormvec = norm(g);
  stepnormvec  = [];
  outflag = 1;
  
  for nu = 1:cg_maxiter
    
    %  Apply SMS preconditioner. Solve M_{SMS}*z = g.

    g1 = restricts(g,p-p0);
    v = prolongs(reshape( AR \ (AR' \ g1(:)), n0,n0), p-p0);
    e = g - amult(v,k_hat,0,L);
    w = reshape( LR\(LR'\e(:)) ,nx,ny);
    e1 = restricts(e,p-p0);
    e1 = prolongs(reshape( GR \ (GR' \ e1(:)), n0,n0), p-p0);
    uQ = (w - e1) / alpha;
    y = g - amult(uQ,k_hat,0,L);
    y1 = restricts(y,p-p0);
    uP = prolongs(reshape( AR \ (AR' \ y1(:)), n0,n0), p-p0);
    z = uP + uQ;

    delta = g(:)'*z(:);
    if nu == 1,
       d = -z; 
    else
       beta_nu = delta / delta_last;
       d = -z + beta_nu * d;
    end

    Ad = amult(d,k_hat,alpha,L);
    alpha_nu = delta / (d(:)'*Ad(:));
    u = u + alpha_nu * d;
    g = g + alpha_nu * Ad;
    delta_last = delta;
    
    residnormvec = [residnormvec; norm(g(:))];
    stepnormvec  = [stepnormvec; abs(alpha_nu)*norm(d(:))/norm(u(:))];
    if exist('u_alpha_true')
      cg_enormvec = [cg_enormvec;  norm(u(:) - u_alpha_true(:))];
    end
    
    fprintf('  CG iter%3.0f, ||resid||=%6.4e, ||step||=%6.4e \n', ... 
         nu, residnormvec(nu), stepnormvec(nu));
    
    if outflag == 1
      figure(2)
        subplot(221)
	imagesc(f_true); colormap(jet), colorbar
	title('True Image')
	subplot(222)
	imagesc(reshape(u,nx,ny)); colormap(jet), colorbar
	title('Approx. soln')
        subplot(223)
        semilogy(residnormvec/residnormvec(1),'o')
        xlabel('CG iteration')
        title('CG residual norm')
        subplot(224)
        semilogy(stepnormvec,'o')
        xlabel('CG iteration')
        title('CG relative step norm')
      drawnow
    end
    
  end



    
