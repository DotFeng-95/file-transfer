%  GPRN_lhd
%
%  Gradient Projection / Reduced Newton algorithm to solve the
% nonnegatively constrained likelihood minimization problem. Minimize 
%    J(f) = sum(K*f+sigsq - (d+sigsq).*log(K*f+sigsq)) + alpha/2*f'*L*f
%  subject to f >= 0.

  alpha = input(' Regularization parameter alpha = ');
  sigsq = input(' Parameter sigma^2 = ');
  max_GPN_iter = input(' Max. no. of GPN iterations = ');
  max_iter_GP = input(' Max. no. of gradient projection iters / GPN = ');
  reg_flag = 0;
  pgrad_tol = 1e-9;

  if reg_flag == 0      %  Identity regularization operator
    L = eye(n);
  elseif reg_flag == 1  %  Negative Laplacian regularization operator
    L = spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], [-1 0 1], n,n);
    L(1,1) = 1;  L(n,n) = 1;
  end

  dn = max(d,0);
  f0_lhd = zeros(n,1);  %  Initial guess.
  termcode = 0;
  k = 0;
  f_lhd = f0_lhd;
  Kf = K*f_lhd;
  g = K'*diag(1./(Kf+sigsq))*(Kf-dn) + alpha*L*f_lhd;
  pgrad = proj_grad(g,f_lhd);      %  Compute norm of projected gradient.
  norm_pgrad = norm(pgrad);
  pgnorm_vec = norm_pgrad;
  stepnormvec = [];
  fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
  
  while termcode == 0
    k = k + 1;
    
    %  Gradient Projection iterations.
    
    [f_lhd,k_GP,termcode_GP] = ...
	grad_proj_lhd(f_lhd,K,dn,L,alpha,sigsq,max_iter_GP);
    
    %  Solve H_reduced * p = D_I * g.
    
    Kf = K*f_lhd;
    g = K'*diag(1./(Kf+sigsq))*(Kf-dn) + alpha*L*f_lhd;
    H = K'*diag((dn+sigsq)./(Kf+sigsq).^2)*K + alpha*L;
    Active = (f_lhd==0);
    D_A = diag(Active);
    D_I = eye(n) - D_A;
    H_reduced = D_I * H * D_I + D_A;
    g_reduced = D_I*g;
    p = -H_reduced \ g_reduced;
    
    %  Line search.
    
    tau0 = -g'*p / (p'*H*p);
    J0 = sum(Kf+sigsq - (dn+sigsq).*log(Kf+sigsq)) + alpha/2*f_lhd'*L*f_lhd;
    [tau,f_lhd,termcode_ls] = ...
	line_srch_lhd(tau0,J0,f_lhd,g,p, K,dn,L,alpha,sigsq);

    %  Compute norm of projected gradient.
    
    Kf = K*f_lhd;
    g = K'*diag(1./(Kf+sigsq))*(Kf-dn) + alpha*L*f_lhd;
    pgrad = proj_grad(g,f_lhd);
    norm_pgrad = norm(pgrad);
    pgnorm_vec = [pgnorm_vec; norm_pgrad];
    stepnorm = norm(f_lhd-f0_lhd);
    stepnormvec = [stepnormvec; stepnorm];

    f0_lhd = f_lhd;
    
    if k >= max_GPN_iter
      termcode = 1;
    elseif norm_pgrad < pgrad_tol
      termcode = 2;
    elseif termcode_GP < 0 | termcode_GP > 3
      termcode = 3;
    elseif termcode_ls > 0
      termcode = 4;
    end

    %  Display numerical performance info and reconstruction.
    
    fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
    figure(3)
      plot(x,f_lhd, x,f_true,'--')
      xlabel('x axis')
      title('Regularized Nonnegative Likelihood Image')
  end
  
  figure(4)
    indx = [0:k]';
    subplot(221)
      semilogy(indx,pgnorm_vec,'o-')
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
    
  rel_error_lhd = norm(f_lhd-f_true) / norm(f_true)