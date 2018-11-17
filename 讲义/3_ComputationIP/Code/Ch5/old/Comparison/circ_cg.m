  function [uh,circ_enormvec] = circ_cg(k_hat,bh,L,alpha,cg_maxiter, ...
      inv_hat,u_exact,u_alpha_true)

%
%  Solve system A*u=b using PCG iteration with block circulant extension
%  preconditioner.

  residnormvec = norm(bh);
  stepnormvec  = [];
  circ_enormvec = [];
  [nx,ny] = size(u_exact);

  %  CG Initialization
  
  uh = zeros(nx,ny);
  g = -bh;
  
  %  Iteration
  
  for cg_iter = 1:cg_maxiter
    
    %  Apply preconditioner.
    
    z = convolve_2d(inv_hat,g);
    z = z(1:nx,1:ny);
    
    delta = g(:)'*z(:);
    if cg_iter == 1
      d = -z;
    else
      beta = delta / delta0;
      d = -z + beta*d;
    end

    Hd = amult(d,k_hat,alpha,L);
    tau = delta / (d(:)'*Hd(:));
    uh = uh + tau*d;
    g = g + tau*Hd;
    delta0 = delta;
    
    %  Compute error indicators and display results.
    
    residnormvec = [residnormvec; norm(g(:))];
    stepnormvec  = [stepnormvec; abs(tau)*norm(d(:))/norm(uh(:))];
    if exist('u_alpha_true')
      circ_enormvec = [circ_enormvec; norm(uh(:) - u_alpha_true(:))];
    end
    fprintf('  CG iter%3.0f, ||resid||=%6.4e, ||step||=%6.4e \n', ... 
         cg_iter, residnormvec(cg_iter), stepnormvec(cg_iter));
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
