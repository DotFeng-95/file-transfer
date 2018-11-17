%  GP_ls
%
%  Gradient Projection algorithm to solve nonnegatively
%  constrained least squares minimization problem. Minimize 
%    J(f) = 0.5*||K*f-d||^2 + alpha/2 * f'*L*f
%         = 0.5 * f'*H*f - f'*b
%  subject to f >= 0, where
%    H = K'*K + alpha*L
%    b = K'*d.

  alpha = input(' Regularization parameter alpha = ');
  max_GP_iter = input(' Max. no. of gradient projection iterations = ');
  reg_flag = 0;
  pgrad_tol = 1e-10;

  if reg_flag == 0      %  Identity regularization operator
    L = eye(n);
  elseif reg_flag == 1  %  Negative Laplacian regularization operator
    L = spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], [-1 0 1], n,n);
    L(1,1) = 1;  L(n,n) = 1;
  end
  
  H = K'*K + alpha*L;
  b = K'*d;
  f0 = zeros(n,1);
  termcode = 0;
  k = 0;
  f = f0;
  g = H*f-b;
  norm_pgrad = norm_projgrad(g,f);
  pgnorm_vec = norm_pgrad;
  if exist('f_ls','var')
    enorm_vec = norm(f - f_ls) / norm(f_ls);
  end
  stepnormvec = [];
  fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
  
  %  Gradient Projection iteration.
  
  while termcode == 0
    k = k + 1;
    Active = (f==0);
    p = -g.*((1-Active)+Active.*(g<0));
    
    %  Line search.
    
    tau0 = -g'*p / (p'*H*p);
    J = 0.5*f'*H*f - f'*b;
    [f,J,tau,termcode_ls] = line_srch_ls(f,J,tau0,g,p,H,b);

    %  Compute norm of projected gradient.
    
    g = H*f - b;
    norm_pgrad = norm_projgrad(g,f);
    pgnorm_vec = [pgnorm_vec; norm_pgrad];
    stepnorm = norm(f-f0);
    stepnormvec = [stepnormvec; stepnorm];
    if exist('f_ls','var')
      enorm_vec = [enorm_vec; norm(f-f_ls)/norm(f_ls)];
    end
    
    f0 = f;
    
    if k >= max_GP_iter
      termcode = 1;
    elseif norm_pgrad < pgrad_tol
      termcode = 2;
    elseif termcode_ls == 4
      termcode = 3;
    end
    
    %  Display numerical performance info and reconstruction.
    
    fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
    figure(3)
      plot(x,f, x,f_true,'--')
      xlabel('x axis')
      title('Regularized Nonnegative Least Squares Image')
  end
  
  figure(4)
    indx_pg = [0:k]';
    subplot(221)
      semilogy(indx_pg,pgnorm_vec,'o-')
      xlabel('GP Iteration')
      title('Projected Gradient Norm')
    subplot(222)
      semilogy(stepnormvec,'o-')
      xlabel('GP Iteration')
      title('Step Norm')
    if exist('f_ls','var')
      subplot(223)
        semilogy(enorm_vec,'o-')
	xlabel('GP Iteration')
	title('Rel Iterative Soln Error Norm')
    end
    
  rel_error = norm(f-f_true) / norm(f_true)
  