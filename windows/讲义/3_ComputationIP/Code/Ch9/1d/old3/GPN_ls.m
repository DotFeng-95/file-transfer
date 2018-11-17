%  GPN_ls
%
%  Gradient Projection / Newton algorithm to solve nonnegatively
%  constrained least squares minimization problem. Minimize 
%    J(f) = 0.5*||K*f-d||^2 + alpha/2 * f'*L*f
%         = 0.5 * f'*H*f - f'*b
%  subject to f >= 0, where
%    H = K'*K + alpha*L
%    b = K'*d.

  alpha = input(' Regularization parameter alpha = ');
  max_GPCG_iter = input(' Max. no. of GPN iterations = ');
  max_iter_GP = input(' Max. no. of gradient projection iters / GPN = ');
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
  g = H*f0-b;
  pgrad = proj_grad(g,f0);
  norm_pgrad = norm(pgrad);
  pgnorm_vec = norm_pgrad;
  stepnormvec = [];
  fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
  
  %  Gradient Projection/Newton iteration.
  
  while termcode == 0
    k = k + 1;
    
    %  Projected gradient iterations.
    
    [f,k_GP,termcode_GP] = grad_proj_ls(f0,H,b,max_iter_GP);
    
    %  Solve H_reduced * p = D_I * g.
    
    g = H*f - b;
    Active = (f==0);
    D_A = diag(Active);
    D_I = eye(n) - D_A;
    H_reduced = D_I * H * D_I + D_A;
    p = -H_reduced \ (D_I*g);
    
    %  Line search.
    
    tau0 = -g'*p / (p'*H*p);
    J = 0.5*f'*H*f - f'*b;
    [tau,f,termcode_ls] = line_srch_ls(tau0,J,f,g,p,H,b);

    %  Compute norm of projected gradient.
    
    pgrad = proj_grad(g,f);
    norm_pgrad = norm(pgrad);
    pgnorm_vec = [pgnorm_vec; norm_pgrad];
    stepnorm = norm(f-f0);
    stepnormvec = [stepnormvec; stepnorm];
    
    f0 = f;
    
    if k >= max_GPCG_iter
      termcode = 1;
    elseif norm_pgrad < pgrad_tol
      termcode = 2;
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
      xlabel('GPCG Iteration')
      title('Projected Gradient Norm')
    subplot(222)
      semilogy(stepnormvec,'o-')
      xlabel('GPCG Iteration')
      title('Step Norm')
      
  rel_error = norm(f-f_true) / norm(f_true)
  