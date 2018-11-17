%  CG_driver.m
%
%  Use CG iteration without preconditioning to solve Ah*u=b, where
%      Ah = K'*K + alpha*L
%      b = K'*dat

  alpha = input(' Regularization parameter alpha = ');
  cg_maxiter = input(' Max no. CG iterations = ');
  bh = Kstar_dat;

  %  Set up regularization operator.
    
  L = speye(nx*ny); %%%  laplacian(nx,ny);
  
%  CG iteration.

  uh = zeros(nx,ny);
  g = -bh;
  residnormvec = norm(bh);
  stepnormvec  = [];
  outflag = 1;
  
  for nu = 1:cg_maxiter

    delta = norm(g(:))^2;
    if nu == 1,
       d = -g; 
    else
       beta_nu = delta / delta_last;
       d = -g + beta_nu * d;
    end

    Ad = amult(d,k_hat,alpha,L);
    alpha_nu = delta / (d(:)'*Ad(:));
    uh = uh + alpha_nu * d;
    g = g + alpha_nu * Ad;
    delta_last = delta;
    
    residnormvec = [residnormvec; norm(g(:))];
    stepnormvec  = [stepnormvec; abs(alpha_nu)*norm(d(:))/norm(uh(:))];
    if exist('u_alpha_true')
      cg_enormvec = [cg_enormvec;  norm(uh(:) - u_alpha_true(:))];
    end
    
    fprintf('  CG iter%3.0f, ||resid||=%6.4e, ||step||=%6.4e \n', ... 
         nu, residnormvec(nu), stepnormvec(nu));
    
    if outflag == 1
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
    
  end



    
