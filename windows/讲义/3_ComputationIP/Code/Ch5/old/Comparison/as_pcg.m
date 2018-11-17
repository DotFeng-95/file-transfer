  function [u,as_enormvec] = as_pcg(k_hat,b,L,alpha,cg_maxiter,AR,LR,GR, ...
      p,p0,f_true,u_alpha_true)
%
%  Use additive Schwarz PCG iteration to solve Ah*u=b, where
%      Ah = K'*K + alpha*L
%      b = K'*dat

  [nx,ny] = size(b);
  n0 = 2^p0;
  
%  CG iteration.

  u = zeros(nx,ny);
  g = -b;
  residnormvec = norm(b);
  stepnormvec  = [];
  as_enormvec = [];
  outflag = 1;
  
  for nu = 1:cg_maxiter
    
    %  Apply additive Schwarz preconditioner. Solve M_{AS}*z = g.

    g1 = restricts(g,p-p0);
    v = AR \ (AR' \ g1(:));
    w = (GR \ (GR' \ g1(:)))/alpha;
    x = prolongs(reshape(v - w,n0,n0), p-p0);
    y = reshape(LR\(LR'\g(:)),nx,ny)/alpha;
    z = x + y;

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
      as_enormvec = [as_enormvec;  norm(u(:) - u_alpha_true(:))];
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



    
