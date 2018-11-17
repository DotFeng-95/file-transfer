  function [u,schur_enormvec] = schurcg(k_hat,b,L,alpha,cg_maxiter,AR,LR, ...
    p,p0,u_exact,u_alpha_true)

%
%  Solve system A*u=b using Schur complement CG iteration.

  residnormvec = norm(b);
  stepnormvec  = [];
  schur_enormvec = [];
  [nx,ny] = size(u_exact);
  n0 = 2^p0;

  %  CG Initialization
  
  Phistar_b = restricts(b,p-p0);
  u = prolongs(reshape(AR\(AR'\Phistar_b(:)),n0,n0),p-p0);
  Au = amult(u,k_hat,alpha,L);
  g = b - Au;
  y = reshape(LR\(LR'\g(:)),nx,ny);
  d = g;
  z = y;
  delta0 = g(:)'*y(:);
  
  %  Iteration
  
  for cg_iter = 1:cg_maxiter
    Phistar_Az = restricts( amult(z,k_hat,alpha,L), p-p0);
    v = z - prolongs(reshape(AR\(AR'\Phistar_Az(:)),n0,n0),p-p0);
    w = amult(v,k_hat,alpha,L);
    h = reshape(LR\(LR'\w(:)),nx,ny);
    tau = delta0 / (d(:)'*h(:));
    
    u = u + tau*v;
    g = g - tau*w;
    y = y - tau*h;
    
    delta1 = g(:)'*y(:);
    beta = delta1/delta0;
    
    d = g + beta*d;
    z = y + beta*z;
    delta0 = delta1;
    
    %  Compute error indicators and display results.
    
    residnormvec = [residnormvec; norm(g(:))];
    stepnormvec  = [stepnormvec; abs(tau)*norm(v(:))/norm(u(:))];
    if exist('u_alpha_true')
      schur_enormvec = [schur_enormvec; norm(u(:) - u_alpha_true(:))];
    end
    fprintf('  CG iter%3.0f, ||resid||=%6.4e, ||step||=%6.4e \n', ... 
         cg_iter, residnormvec(cg_iter), stepnormvec(cg_iter));
      figure(2)
        subplot(221)
        imagesc(u_exact); colormap(jet), colorbar
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
