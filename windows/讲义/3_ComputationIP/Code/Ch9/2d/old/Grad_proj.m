%  Grad_proj.m
%
%  Use gradient projection to minimize the penalized least squares cost 
%  function
%    J(f) = 1/2*(||S*f-d||^2 + alpha*f'*L*f)
%  subject to the constraint u > or = 0. 

  alpha = 1e-3; %%%input(' Regularization parameter alpha = ');
  max_iter = 500; %%%input(' Max. no. gradient projection iterations = ');
  tol = 1e-5; %%%input(' Gradient norm relative stopping tolerance = ');
  iterhist = [];
  beta = .1;

%  Initial guess.

  f_tikh = real(ifft2(conj(S).*fft2(d)./(abs(S).^2 + 1e-3)));
  f = max(f_tikh,0);

%  Set up regularization operator.

  Lflag = 0; %%%input(' Enter 0 for L^2 reg; 1 for H^1 reg: ');
  if Lflag == 0      %  Identity
    L = speye(Nx*Ny);
  elseif Lflag == 1  %  Negative Laplacian
    L = laplacian(Nx,Ny);
  end
  
  [J,g] = eval_pls(f,S,d,L,alpha);
  Pg = f - max(f-g,zeros(size(f)));
  normPg = norm(Pg(:));
  rel_tol = tol * normPg;
  epsilon0 = max(f_true(:)) / 10;
  epsilon = min(epsilon0,normPg);
  iter = 0;
  
  while (normPg > rel_tol & iter <= max_iter)

    Active = (f==0 & g>0);
    n_Active = sum(Active(:));

    s = -g;
    lambda = 1;
    iter = iter + 1;
    f_new = max(f+lambda*s,zeros(size(f)));
    delta_f = f - f_new;
    J_new = eval_pls(f_new,S,d,L,alpha);
    J_goal = J - 1e-4*(g(:)'*delta_f(:));
    ls_iter = 0;
    
    %  Simple line search.
    
    while (J_new > J_goal)
      lambda = lambda*beta;
      ls_iter = ls_iter + 1;
      f_new = max(f+lambda*s,zeros(size(f)));
      delta_f = f - f_new;
      J_new = eval_pls(f_new,S,d,L,alpha);
      if (ls_iter > 10)
	disp('*** Error in line search ***');
	break
      end
      J_goal = J - 1e-4*(g(:)'*delta_f(:));
    end
    f = f_new;
    [J,g] = eval_pls(f,S,d,L,alpha);
    Pg = f - max(f-g,zeros(size(f)));
    normPg = norm(Pg(:));
    epsilon = min(epsilon0,normPg);

    %  Record and display data.
    
    iterhist(iter,1) = normPg;
    iterhist(iter,2) = ls_iter;
    iterhist(iter,3) = n_Active;
    iterhist(iter,4) = 1;
   if mod(iter,max_iter/20) == 0
   fprintf(' iter=%3.0f ||Pg||=%6.4e #LS=%3.0f #Active=%5.0f\n',...
       iter, normPg, ls_iter, n_Active);
    figure(3)
      subplot(221)
        imagesc(extract(f_true,nx,ny)), colorbar
	title('True Object')
      subplot(222)
        imagesc(extract(f,nx,ny)), colorbar
	title('Nonnegative Reconstruction')
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
    colormap(1-gray)
    drawnow
   end

  end
    
  rel_error = norm(f-f_true,'fro')/norm(f_true,'fro');
  
  fprintf('\n  Relative soln error = %7.4e\n', rel_error);
  
