  function [f,cgiter] = pcg_pls(S,b,alpha,L,Active,...
      maxiter,cg_fig,steptol,residtol,cg_tab_flag)
  
%
%

  [Nx,Ny] = size(b);
  f = zeros(size(b));
  bbar = (1-Active) .* b;
  resid = bbar;
  residnormvec = norm(resid(:));
  stepnormvec = [];
  cgiter = 0;
  stop_flag = 0;
  cgiter = 0;
  stop_flag = 0;
  
%  CG iterations.
  
  while stop_flag == 0

    cgiter = cgiter + 1;
    dh = resid;
    rd = resid(:)'*dh(:);
    
    %  Compute conjugate direction ph.

    if cgiter == 1,
      ph = dh; 
    else
      betak = rd / rdlast;
      ph = dh + betak * ph;
    end

    %  Form product Ah*ph.
    
    Iph = (1-Active).*ph;
    LIph = reshape(L*Iph(:),Nx,Ny);
    Ahph = (1-Active).*(real(ifft2(abs(S).^2 .* fft2(Iph))) + alpha*LIph);
 
    %  Update f and residual.
    
    alphak = rd / (ph(:)'*Ahph(:));
    f = f + alphak*ph;
    resid = resid - alphak*Ahph;
    rdlast = rd;

    residnorm = norm(resid(:));
    stepnorm = abs(alphak)*norm(ph(:))/norm(f(:));
    residnormvec = [residnormvec; residnorm];
    stepnormvec = [stepnormvec; stepnorm];
    
    %  Check stopping criteria.
    
    if cgiter >= maxiter
      stop_flag = 1;
    elseif stepnorm < steptol
      stop_flag = 2;
    elseif residnorm / residnormvec(1) < residtol
      stop_flag = 3;
    end
    
    %  Output convergence information.

    if cg_tab_flag == 1
      fprintf('  CG iter%3.0f, ||resid||=%6.4e, ||step||=%6.4e \n', ... 
         cgiter, residnormvec(cgiter), stepnormvec(cgiter));
    end
    if cg_fig > 0
      figure(cg_fig)
        subplot(221)
          semilogy(residnormvec/residnormvec(1),'o')
          xlabel('CG iteration')
          title('CG residual norm')
        subplot(222)
          semilogy(stepnormvec,'o')
          xlabel('CG iteration')
          title('CG relative step norm')
      drawnow
    end

  end

  f = f + Active.*b; 
