%  GP_lhd
%
%  Gradient Projection algorithm to solve the nonnegatively
%  constrained likelihood minimization problem. Minimize 
%    J(f) = sum(K*f+sigsq - (d+sigsq).*log(K*f+sigsq)) + alpha/2*f'*L*f
%  subject to f >= 0.

  alpha = input(' Regularization parameter alpha = ');
  sigsq = input(' Parameter sigma^2 = ');
  max_GP_iter = input(' Max. no. of gradient projection iterations = ');
  reg_flag = 0;
  pgrad_tol = 1e-9;

  if reg_flag == 0      %  Identity regularization operator
    L = eye(n);
  elseif reg_flag == 1  %  Negative Laplacian regularization operator
    L = spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], [-1 0 1], n,n);
    L(1,1) = 1;  L(n,n) = 1;
  end

  dn = max(d,0);
  f0 = zeros(n,1);  %  Initial guess.
  termcode = 0;
  k = 0;
  f = f0;
  Kf = K*f;
  g = K'*diag(1./(Kf+sigsq))*(Kf-dn) + alpha*L*f;
  pgrad = proj_grad(g,f);      %  Compute norm of projected gradient.
  norm_pgrad = norm(pgrad);
  pgnorm_vec = norm_pgrad;
  stepnormvec = [];
  if exist('f_lhd','var')
    enorm_vec = norm(f - f_lhd) / norm(f_lhd);
  end
  
  fprintf(' iter%4.0f, |Projected grad| = %6.4e\n', k, norm_pgrad);
  
  %  Gradient Projection iterations.
  
  while termcode == 0
    k = k + 1;
    Active = (f==0);
    p = -g.*((1-Active)+Active.*(g<0));
    
    %  Line search.
    
    H = K'*diag((dn+sigsq)./(Kf+sigsq).^2)*K + alpha*L;
    tau0 = -g'*p / (p'*H*p);
    J0 = sum(Kf+sigsq - (dn+sigsq).*log(Kf+sigsq)) + alpha/2*f'*L*f;
    [tau,f,termcode_ls] = line_srch_lhd(tau0,J0,f,g,p, K,dn,L,alpha,sigsq);

    %  Compute norm of projected gradient.
    
    Kf = K*f;
    g = K'*diag(1./(Kf+sigsq))*(Kf-dn) + alpha*L*f;
    pgrad = proj_grad(g,f);
    norm_pgrad = norm(pgrad);
    pgnorm_vec = [pgnorm_vec; norm_pgrad];
    stepnorm = norm(f-f0);
    stepnormvec = [stepnormvec; stepnorm];
    if exist('f_lhd','var')
      enorm_vec = [enorm_vec; norm(f-f_lhd)/norm(f_lhd)];
    end
    
    f0 = f;
    
    if k >= max_GP_iter
      termcode = 1;
    elseif norm_pgrad < pgrad_tol
      termcode = 2;
    elseif termcode_ls > 2
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
    subplot(221)
      semilogy(pgnorm_vec,'o-')
      xlabel('GP Iteration')
      title('Projected Gradient Norm')
    subplot(222)
      semilogy(stepnormvec,'o-')
      xlabel('GP Iteration')
      title('Step Norm')
    if exist('f_lhd','var')
      subplot(223)
        indx = [0:k]';
        semilogy(indx,enorm_vec,'o-')
        xlabel('GP Iteration')
        title('Rel Iterative Soln Error Norm')
    end
    
  rel_error = norm(f-f_true) / norm(f_true)
  