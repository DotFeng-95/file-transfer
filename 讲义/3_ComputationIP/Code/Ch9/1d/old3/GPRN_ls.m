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
  max_GPN_iter = input(' Max. no. of GPN iterations = ');
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
  f0_ls = zeros(n,1);
  termcode = 0;
  k = 0;
  g = H*f0_ls-b;
  pgrad = proj_grad(g,f0_ls);
  norm_pgrad = norm(pgrad);
  pgnorm_vec = norm_pgrad;
  stepnormvec = [];
  fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
  
  %  Gradient Projection/Newton iteration.
  
  while termcode == 0
    k = k + 1;
    
    %  Projected gradient iterations.
    
    [f_ls,k_GP,termcode_GP] = grad_proj_ls(f0_ls,H,b,max_iter_GP);
    
    %  Solve H_reduced * p = D_I * g.
    
    g = H*f_ls - b;
    Active = (f_ls==0);
    D_A = diag(Active);
    D_I = eye(n) - D_A;
    H_reduced = D_I * H * D_I + D_A;
    p = -H_reduced \ (D_I*g);
    
    %  Line search.
    
    tau0 = -g'*p / (p'*H*p);
    J = 0.5*f_ls'*H*f_ls - f_ls'*b;
    [tau,f_ls,termcode_ls] = line_srch_ls(tau0,J,f_ls,g,p,H,b);

    %  Compute norm of projected gradient.
    
    pgrad = proj_grad(g,f_ls);
    norm_pgrad = norm(pgrad);
    pgnorm_vec = [pgnorm_vec; norm_pgrad];
    stepnorm = norm(f_ls-f0_ls);
    stepnormvec = [stepnormvec; stepnorm];
    
    f0_ls = f_ls;
    
    if k >= max_GPN_iter
      termcode = 1;
    elseif norm_pgrad < pgrad_tol
      termcode = 2;
    end
    
    %  Display numerical performance info and reconstruction.
    
    fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
    figure(3)
      plot(x,f_ls, x,f_true,'--')
      xlabel('x axis')
      title('Regularized Nonnegative Least Squares Image')
  end
  
  figure(4)
    indx_pg = [0:k]';
    subplot(221)
      semilogy(indx_pg,pgnorm_vec,'o-')
      xlabel('GPN Iteration')
      title('Projected Gradient Norm')
    subplot(222)
      semilogy(stepnormvec,'o-')
      xlabel('GPN Iteration')
      title('Step Norm')
    subplot(223)
      hist(log10(g(Active)))
      xlabel('log_{10}(g_{active})')
      ylabel('No. of Occurences')
      title('Hist of active grad components')
      
  rel_error_ls = norm(f_ls-f_true) / norm(f_true)
  