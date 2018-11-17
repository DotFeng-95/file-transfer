%  schur1.m
%
%  Solve system A*u=b using Schur complement CG iteration.

  alpha = input(' Regularization parameter alpha = ');
  cg_maxiter = input(' Max no. CG iterations = ');
  resetflag = input(' Reset preconditioner? (0 for No, 1 for Yes): ');
  if resetflag == 1
    [nx,ny] = size(z);
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
    
    %  Set up A11 matrix and compute its Choleski factorization.
    
    disp(' ... Setting up A11 ...');
    A11 = zeros(N,N);
    for i = 1:N
      phi_i = zeros(N,1);
      phi_i(i) = 1;
      phi_i = prolongs(reshape(phi_i,n0,n0), p-p0);
      APhi_i = amult(phi_i,k_hat,alpha,L);
      Phistar_APhi_i = restricts(APhi_i,p-p0);
      A11(:,i) = Phistar_APhi_i(:);
    end
    AR = chol(A11);
    LR = chol(L);
    disp(' ... A11 setup complete ...')
  end
  
  residnormvec = [];
  stepnormvec  = [];
  schur_enormvec = [];

  %  CG Initialization
  
  Phistar_b = restricts(b,p-p0);
  uh = prolongs(reshape(AR\(AR'\Phistar_b(:)),n0,n0),p-p0);
  Au = amult(uh,k_hat,alpha,L);
  resid = b - Au;
  yh = reshape(LR\(LR'\resid(:)),nx,ny);
  dh = resid;
  zh = yh;
  delta0 = resid(:)'*yh(:);
  
  %  Iteration
  
  for cg_iter = 1:cg_maxiter
    Az = amult(zh,k_hat,alpha,L);
    Phistar_Az = restricts(Az,p-p0);
    v = zh-prolongs(reshape(AR\(AR'\Phistar_Az(:)),n0,n0),p-p0);
    w = amult(v,k_hat,alpha,L);
    g = reshape(LR\(LR'\w(:)),nx,ny);
    tau = delta0 / (dh(:)'*g(:));
    
    uh = uh + tau*v;
    resid = resid - tau*w;
    yh = yh - tau*g;
    
    delta1 = resid(:)'*yh(:);
    beta = delta1/delta0;
    delta0 = delta1;
    
    dh = resid + beta*dh;
    zh = yh + beta*zh;
    
    %  Compute error indicators and display results.
    
    residnormvec = [residnormvec; norm(resid(:))];
    stepnormvec  = [stepnormvec; abs(tau)*norm(v(:))/norm(uh(:))];
    if exist('u_alpha_true')
      schur_enormvec = [schur_enormvec; norm(uh(:) - u_alpha_true(:))];
    end
    fprintf('  CG iter%3.0f, ||resid||=%6.4e, ||step||=%6.4e \n', ... 
         cg_iter, residnormvec(cg_iter), stepnormvec(cg_iter));
      figure(2)
        subplot(221)
        imagesc(f_true); colormap(jet), colorbar
        title('True Image')
        subplot(222)
        imagesc(reshape(uh,nx,ny)); colormap(jet), colorbar
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
