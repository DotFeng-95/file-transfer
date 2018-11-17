%  Lhd_gradproj.m
%
%  Compute minimizer of penalized log likelihood functional
%     J(f) =  sum[Kf+sigsq - (d+sigsq)*log(Kf+sigsq) + alpha/2*f'*L*f
%  subject to [f]_i \geq 0 for each i.
%  Use the Gradient Projection method.
%  gradient:             g(f) = K'*diag(1./(Kf+sigsq))*(Kf-d) + alpha*L*f
%  Hessian:              H = K'*diag[(d+sigsq)./(Kf+sigsq).^2]*K + alpha*L

  alpha = 1e-7; %%%input(' Regularization parameter alpha = ');
  max_iter = 500; %%%input(' Max. no. Gradient Projection iterations = ');
  tol = 1e-5; %%%input(' Projected gradient relative stopping tolerance = ');
  sigsq = 9; %%%input(' sigma^2 = ');
  iterhist = [];
  beta = .5;

%  Set up penalty operator L.
  
  [Nx,Ny] = size(d);
  Lflag = 1; %%%input(' Enter 0 for L^2 reg.; else enter 1 for H^1 reg: ');
  if Lflag == 0
    L = speye(Nx*Ny);      %  Identity
  elseif Lflag ==1
    L = laplacian(Nx,Ny);  %  Negative Laplacian
  end

%  Initial guess.

  f_tik = real(ifft2(conj(S).*fft2(d)./(abs(S).^2+5e-3)));
  f = max(f_tik,0);
  dps = d + sigsq;
  [J,Kfps,g] = eval_lhd(f,S,dps,alpha,L,sigsq);
  D = dps ./ Kfps.^2;
  Pg = f - max(f-g,zeros(size(f)));
  normPg = norm(Pg(:));
  iter = 0;
  
  while (normPg>rel_tol & iter <= max_iter)
    iter = iter + 1;
    Active = ((0 < f & f <= epsilon) | (f==0 & g>0));
    n_Active = sum(Active(:));
    s = -g;
    lambda = 1;
    f_new = max(f+lambda*s,zeros(size(f)));
    delta_f = f - f_new;
    J_new = eval_lhd(f_new,S,dps,alpha,L,sigsq);
    J_goal = J - 1e-4*(g(:)'*delta_f(:));
    ls_iter = 0;
    
    %  Simple line search.
    
    while (J_new > J_goal)
      lambda = lambda*beta;
      ls_iter = ls_iter + 1;
      f_new = max(f+lambda*s,zeros(size(f)));
      delta_f = f - f_new;
      J_new = eval_lhd(f_new,S,dps,alpha,L,sigsq);
      if (ls_iter > 10)
	disp('*** Error in line search ***');
	break
      end
      J_goal = J - 1e-4*(g(:)'*delta_f(:));
    end
    f = f_new;
    [J,Kfps,g] = eval_lhd(f,S,dps,alpha,L,sigsq);
    D = dps ./ Kfps.^2;
    Pg = f - max(f-g,zeros(size(f)));
    normPg = norm(Pg(:));
    epsilon = min(epsilon0,normPg);

    %  Record and display data.
    
    iterhist(iter,1) = normPg;
    iterhist(iter,2) = ls_iter;
    iterhist(iter,3) = n_Active;
    iterhist(iter,4) = cgiter;

   if mod(iter,max_iter/20) == 0
   fprintf(' iter=%3.0f ||Pg||=%6.4e #LS=%3.0f #Active=%5.0f\n',...
       iter, normPg, ls_iter, n_Active);
    figure(3)
      subplot(221)
        imagesc(extract(f_true,nx,ny)), colorbar
	title('True Object')
      subplot(222)
        imagesc(extract(f,nx,ny)), colorbar
	title('Nonnegative LHD Reconstruction')
      subplot(223)
        yy = iterhist(:,3)/prod(size(f));
	xx = [1:max(size(yy))]'; 
        plot(xx,yy,'o', xx,yy,'-')
        xlabel('Iteration')
        title('Fraction of Constraints Active')
      subplot(224)
        yy = iterhist(:,1);
	xx = [1:max(size(yy))]';
        semilogy(xx,yy,'o', xx,yy,'-')
	xlabel('Iteration')
	title('Projected Gradient Norm')
    colormap(1-hot)
    drawnow
   end
  end

  rel_error = norm(f-f_true,'fro')/norm(f_true,'fro');
  
  fprintf('\n  Relative soln error = %7.4e\n', rel_error);
