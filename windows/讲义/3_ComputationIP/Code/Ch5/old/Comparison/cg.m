  function [uh,cg_enormvec] = ... 
    conjgrad(k_hat,bh,L,alpha,cg_maxiter,u_exact,u_alpha_true)
%
%  Use Conjugate Gradient iteration 
%  to solve Ah*u=b, where
%      Ah = K'*K + alpha*L
%      b = K'*z

  [nx,ny] = size(u_exact);  
  outflag = 1;
  
%  CG iteration.

  uh = zeros(nx,ny);
  residnormvec = norm(bh);
  stepnormvec  = [];
  cg_enormvec = [];
  
  for cg_iter = 1:cg_maxiter
    
    %  Form product Ah*uh.

    Ahuh = amult(uh,k_hat,alpha,L);
    
    resid = bh - Ahuh;
    d = resid;
    rd = resid(:)'*d(:);
    if cg_iter == 1,
       pvec = d; 
    else
       betak = rd / rdlast;
       pvec = d + betak * pvec;
    end

%  Compute Ah*p

    Ap = amult(pvec,k_hat,alpha,L);

    alphak = rd / (pvec(:)'*Ap(:));
    uh = uh + alphak*pvec;
    resid = resid - alphak*Ap;
    rdlast = rd;
    residnormvec = [residnormvec; norm(resid(:))];
    stepnormvec  = [stepnormvec; abs(alphak)*norm(pvec(:))/norm(uh(:))];
    if exist('u_alpha_true')
      cg_enormvec = [cg_enormvec;  norm(uh(:) - u_alpha_true(:))];
    end
    
    fprintf('  CG iter%3.0f, ||resid||=%6.4e, ||step||=%6.4e \n', ... 
         cg_iter, residnormvec(cg_iter), stepnormvec(cg_iter));
    
    if outflag == 1
      figure(2)
        subplot(221)
	imagesc(u_exact); colormap(jet), colorbar
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



    
