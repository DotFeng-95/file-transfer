%  MRNSD.m
%
%  Use the Modified Residual Norm Steepest Descent (MRNSD) algorithm
%  to approximately solve 
%    min 0.5*||Kf - d||^2  subject to  f>=0

  max_iter = 100; %%%input(' Max. no. MRNSD iterations = ');

%  Initialization

  f_MRNSD = ones(size(f_true));
  Kf = mat_prod(S{1},f_MRNSD);
  g_k = mat_prod(conj(S{1}),Kf-dat);
  norm_f = norm(f_true(:));
  MRNSD_error = norm(f_MRNSD(:)-f_true(:))/norm_f;

%  Begin iteration

  for k=1:max_iter
    
    p_k = -f_MRNSD .* g_k;        %  Descent direction.
    gam = -g_k(:)'*p_k(:);
    u = mat_prod(S{1},p_k);
    
    %  Compute tau_bndry
    
    neg_indx = (p_k < 0);
    if sum(neg_indx(:)) == 0
      tau_bndry = inf;
    else
      tau_bndry = min( -f_MRNSD(neg_indx) ./ p_k(neg_indx) );
    end
    
    %  Compute tau_uncnstr
    
    tau_uncnstr = gam / norm(u(:))^2;
    
    tau = min(tau_bndry,tau_uncnstr);
    f_MRNSD = f_MRNSD + tau*p_k;
    g_k = g_k + tau*mat_prod(conj(S{1}),u);
    ek = norm(f_MRNSD(:)-f_true(:))/norm_f;
    MRNSD_error = [MRNSD_error; ek];
    
    %  Display MRNSD reconstruction
    
    if mod(k,max_iter/20) == 0
      fprintf('   Iteration %i, relative error %6.4e\n',k,ek);
      figure(1)
        subplot(221)
          imagesc(f_true(1:nfx,1:nfy)), colorbar
          title('True Object')
        subplot(222)
          imagesc(f_MRNSD(1:nfx,1:nfy)), colorbar
          title('MRNSD Reconstruction')
        subplot(223)
          semilogy(MRNSD_error)
          xlabel('MRNSD Iteration')
          title('MRNSD Relative Solution Error')
      colormap(1-gray)
      drawnow
    end
    
  end
    
  [mval,mindex] = min(MRNSD_error);
  fprintf(' Min relative soln error of %5.2f percent occurred at MRNSD iteration %d\n', mval*100,mindex);






