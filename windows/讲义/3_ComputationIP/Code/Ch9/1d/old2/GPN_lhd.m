%  GPN_lhd
%
%  Gradient Projection / Newton algorithm to solve the nonnegatively
%  constrained likelihood minimization problem. Minimize 
%    J(f) = sum(K*f+sigsq - (d+sigsq).*log(K*f+sigsq)) + alpha/2*f'*L*f
%  subject to f >= 0.

  alpha = input(' Regularization parameter alpha = ');
  sigsq = input(' Parameter sigma^2 = ');
  max_GPCG_iter = input(' Max. no. of GPCG iterations = ');
  max_iter_GP = 10;     %  Max. no. of gradient projection iterations.
  reg_flag = 0;

  if reg_flag == 0      %  Identity regularization operator
    L = eye(n);
  elseif reg_flag == 1  %  Negative Laplacian regularization operator
    L = spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], [-1 0 1], n,n);
    L(1,1) = 1;  L(n,n) = 1;
  end

  dn = max(d,0);
  f = zeros(n,1);  %  Initial guess.
%f = zeros(n,1);
  termcode = 0;
  k = 0;
  Kf = K*f;
  g = K'*diag(1./(Kf+sigsq))*(Kf-dn) + alpha*L*f;
  pgrad = f - proj(f-g);      %  Compute norm of projected gradient.
  norm_pgrad = norm(pgrad);
  pgnorm_vec = norm_pgrad;
  fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
  
  while termcode == 0
    k = k + 1;
    
    %  Projected gradient iterations.
    
    [f,k_GP,termcode_GP] = grad_proj_lhd(f,K,dn,L,alpha,sigsq,max_iter_GP);
    
    %  Solve H_reduced * p = D_I * g.
    
    Kf = K*f;
    g = K'*diag(1./(Kf+sigsq))*(Kf-dn) + alpha*L*f;
    H = K'*diag((dn+sigsq)./(Kf+sigsq).^2)*K + alpha*L;
    Active = (f==0);
    D_A = diag(Active);
    D_I = eye(n) - D_A;
    H_reduced = D_I * H * D_I + D_A;
%    p = cg(zeros(n,1),H_reduced,g_reduced,max_iter_CG);
    p = -H_reduced \ (D_I*g);
    
    %  Line search.
    
    tau0 = -g'*p / (p'*H*p);
    J0 = sum(Kf+sigsq - (dn+sigsq).*log(Kf+sigsq)) + alpha/2*f'*L*f;
    [tau,f] = line_srch_lhd(tau0,J0,f,g,p, K,dn,L,alpha,sigsq);

    if k >= max_GPCG_iter
      termcode = 1;
    end

    %  Compute norm of projected gradient.
    
    Kf = K*f;
    g = K'*diag(1./(Kf+sigsq))*(Kf-dn) + alpha*L*f;
    pgrad = f - proj(f-g);
    norm_pgrad = norm(pgrad);
    pgnorm_vec = [pgnorm_vec; norm_pgrad];
    
    %  Display numerical performance info and reconstruction.
    
    fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
    figure(3)
      plot(x,f, x,f_true,'--')
      xlabel('x axis')
      title('Regularized Nonnegative Least Squares Image')
  end
  
  figure(4)
    semilogy(pgnorm_vec,'o-')
    xlabel('GPCG Iteration')
    title('Projected Gradient Norm vs Iteration Count')
    
  rel_error = norm(f-f_true) / norm(f_true)