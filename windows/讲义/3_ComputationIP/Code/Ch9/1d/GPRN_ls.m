%  GPRN_ls
%
%  Gradient Projection / reduced Newton algorithm to solve nonnegatively
%  constrained regularized least squares minimization problem. Minimize 
%    J(f) = 0.5*||K*f-d||^2 + alpha/2 * f'*L*f
%         = 0.5 * f'*H*f - f'*b
%  subject to f >= 0, where
%    H = K'*K + alpha*L
%    b = K'*d.

  alpha = 5e-3; %%%input(' Regularization parameter alpha = ');
  max_GPRNiter = 20; %%%input(' Max. no. of GPRN iterations = ');
  max_GPiter = 5;%input(' Max. no. of Stage I (gradient projection) iters = ');
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
  f_ls = zeros(n,1);
  J = 0.5 * f_ls'*H*f_ls - f_ls'*b;
  g = H*f_ls-b;
  norm_pgrad = norm_projgrad(g,f_ls);
  pgnorm_vec = norm_pgrad;
  stepnormvec = [];
  termcode = 0;
  k = 0;
  fprintf(' outer iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
  
  %  Gradient Projection / Reduced Newton iteration.
  
  while termcode == 0
    k = k + 1;
    
    %  Stage I: Perform projected gradient iterations.
    
    %  Initialization.
    
    sigma = 0.25;
    J0_GP = J;
    max_delJ = 0;
    Active0 = (f_ls==0);  %  Current active set of indices.
    GPiter = 0;
    GPtermcode = 0;
    f_GP = f_ls;
  
    while GPtermcode == 0
      GPiter = GPiter + 1;    %  GP iteration count.
      p = -g.*((1-Active0)+Active0.*(g<0));
      tau0 = -(p'*g) / (p'*H*p);
      [f_GP,J_GP,tau,ls_termcode] = line_srch_ls(f_GP,J0_GP,tau0,g,p,H,b);
      g = H*f_GP - b;
      delta_J = J0_GP - J_GP;
      max_delJ = max(delta_J,max_delJ);
      Active = (f_GP==0);
    
      %  Check GP stopping criteria
    
      if delta_J < sigma * max_delJ    %  Insufficient decrease in J
        GPtermcode = 1;
      elseif norm(Active-Active0) == 0 %  Active set unchanged.
        GPtermcode = 2; 
      elseif GPiter >= max_GPiter      %  Max iteration count exceeded.
        GPtermcode = 3;
      elseif ls_termcode == 4          %  Line search failed.
        GPtermcode = -ls_termcode;
      end
      J0_GP = J_GP;
      Active0 = Active;
      
      %  Output Stage I info.
      
      norm_pgrad = norm_projgrad(g,f_GP);
      fprintf('   GP it %d delJ=%6.4e |Projected grad|= %6.4e LStrmcd=%d\n',...
	   GPiter, delta_J, norm_pgrad, ls_termcode);
    
    end
    fprintf('     End GP iterations. GP termcode = %d\n', GPtermcode);
    
    %  Stage II: Solve H_reduced * p = D_I * g and then perform
    %  projected line search.
    
    f_ls = f_GP;
    J = J_GP;
    Active = (f_ls==0);
    D_A = diag(Active);
    D_I = eye(n) - D_A;
    H_reduced = D_I * H * D_I + D_A;
    p = -H_reduced \ (D_I*g);
    
    %  Line search.
    
    tau0 = -g'*p / (p'*H*p);
    [f_ls,J,tau,termcode_ls] = line_srch_ls(f_ls,J,tau0,g,p,H,b);

    %  Compute gradient and norm of projected gradient.
    
    g = H*f_ls - b;
    norm_pgrad = norm_projgrad(g,f_ls);
    pgnorm_vec = [pgnorm_vec; norm_pgrad];
    stepnorm = norm(f_ls-f_GP);
    stepnormvec = [stepnormvec; stepnorm];
    
    if k >= max_GPRNiter
      termcode = 1;
    elseif norm_pgrad < pgrad_tol
      termcode = 2;
    end
    
    %  Display numerical performance info and reconstruction.
    
    fprintf(' outer iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
    figure(3)
      plot(x,f_ls, x,f_true,'--')
      xlabel('x axis')
      title('Regularized Nonnegative Least Squares Image')
  end
  
  figure(4)
    indx_pg = [0:k]';
    subplot(221)
      semilogy(indx_pg,pgnorm_vec,'o-')
      xlabel('GPRN Iteration')
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
