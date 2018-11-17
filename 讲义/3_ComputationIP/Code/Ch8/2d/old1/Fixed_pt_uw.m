%  Fixed_pt_uw.m    (Upwind difference gradient approximation used)
%
%  Use "lagged diffusivity" fixed point iteration to minimize
%      T(u) = ||K*u - d||^2/2 + alpha*J(u),
%  where K is a discretized integral operator, d is discrete data, 
%  ||.|| denotes the l^2 norm, alpha is a positive regularization 
%  parameter, and J is a smooth approximation to the 
%  Total Variation functional.
%      J(u) = sum_i psi(|[D*u]_i|^2,beta) * Delta_x,
%  where D is a discretization of the first derivative operator and
%  beta is a positive smoothing parameter.
%
%  At each iteration, replace u by u + Delta_u, where Delta_u solves
%    (K'*K + alpha*L(u)) * Delta_u = K'(K*u-d) + alpha*L(u)*u,
%  where 
%    L(u)*v = D'* diag(psi'(|[D*u]_i|^2,beta) * D * Delta_x. 

  if ~exist('max_fp_iter')
    alpha = input(' Regularization parameter alpha = ');
    beta = input(' TV smoothing parameter beta = ');
    max_fp_iter = input(' Max. no. of fixed point iterations = ');
    max_cg_iter = input(' Max. no. of CG iterations = ');
    cg_steptol = 1e-4;
    cg_residtol = 1e-4;
    cg_out_flag = 0;  %  If flag = 1, output CG convergence info.
  end

  %  Discretize gradient operator using upwind differencing.
  
  n = nfx;
  nsq = n^2;
  Delta_x = 1 / n;
  Delta_y = Delta_x;
  D = spdiags([-ones(n-1,1) ones(n-1,1)], [0 1], n-1,n) / Delta_x;
  I_trunc = spdiags(ones(n-1,1), 0, n-1,n);
  Dx = kron(D,I_trunc);
  Dy = kron(I_trunc,D);

  %  Initialization.

  fp_gradnorm = [];
  snorm_vec = [];
  k_hat_sq = abs(k_hat).^2;
  Kstar_d = integral_op(dat,conj(k_hat),n,n);   %  Compute K'*dat.
  U_fp = zeros(n,n);
  
  for fp_iter = 1:max_fp_iter
      
    %  Set up regularization operator L.

    u_fp = U_fp(:);
    psi_1 = psi_prime((Dx*u_fp).^2 + (Dy*u_fp).^2, beta);
    diff_coef = spdiags(psi_1, 0, (n-1)^2,(n-1)^2);
    L = Dx' * diff_coef * Dx + Dy' * diff_coef * Dy;
    L = L * Delta_x * Delta_y / 2;
    KstarKu = integral_op(U_fp,k_hat_sq,n,n);    
    G = KstarKu - Kstar_d + alpha* reshape(L*u_fp,n,n);
    gradnorm = norm(G(:));
    fp_gradnorm = [fp_gradnorm; gradnorm];
    
    %  CG initialization.
  
    Delta_U = zeros(n,n);
    resid = -G;
    residnormvec = norm(resid(:));
    stepnormvec = [];
    cgiter = 0;
    stop_flag = 0;

    %  CG iteration.

    while stop_flag == 0
      cgiter = cgiter + 1;

      dh = resid;

      %  Compute conjugate direction p.
 
      rd = resid(:)'*dh(:);
      if cgiter == 1,
         ph = dh; 
       else
         betak = rd / rdlast;
         ph = dh + betak * ph;
      end

      %  Form product Ah*ph.

      KstarKp = integral_op(ph,k_hat_sq,n,n);
      Ahph = KstarKp + alpha * reshape(L*ph(:),n,n);

      %  Update Delta_U and residual.
    
      alphak = rd / (ph(:)'*Ahph(:));
      Delta_U = Delta_U + alphak*ph;
      resid = resid - alphak*Ahph;
      rdlast = rd;
      
      residnorm = norm(resid(:));
      stepnorm = abs(alphak)*norm(ph(:))/norm(Delta_U(:));
      residnormvec = [residnormvec; residnorm];
      stepnormvec = [stepnormvec; stepnorm];

      %  Check stopping criteria.
    
      if cgiter >= max_cg_iter
        stop_flag = 1;
      elseif stepnorm < cg_steptol
        stop_flag = 2;
      elseif residnorm / residnormvec(1) < cg_residtol
        stop_flag = 3;
      end

      %  Display CG convergence information.
      
      if cg_out_flag == 1
          fprintf('   CG iter%3.0f, ||resid||=%6.4e, ||step||=%6.4e \n', ... 
            cgiter, residnorm, stepnorm)
        figure(1)
          subplot(221)
            semilogy(residnormvec/residnormvec(1),'o')
            title('CG Residual Norm')
          subplot(222)
            semilogy(stepnormvec,'o')
            title('CG Step Norm')
      else
	if cgiter == 1
	  fprintf(' ... computing Delta_u via CG iterations ...\n');
	end
      end
      
    end %end for CG

    %  Update TV regularized reconstruction.
    
    U_fp = U_fp + Delta_U;
    snorm = norm(Delta_U(:));
    snorm_vec = [snorm_vec; snorm];
    
    %   Output fixed point convergence information.
   
    fprintf(' FP iter=%3.0f, ||grad||=%6.4e, ||step||=%6.4e, nCG=%3.0f\n', ...
       fp_iter, gradnorm, snorm, cgiter);
   
    figure(2)
      subplot(221)
        semilogy(residnormvec/residnormvec(1),'o')
        xlabel('CG iteration')
        title('CG Relative Residual Norm')
      subplot(222)
        semilogy(stepnormvec,'o')
        xlabel('CG iteration')
        title('CG Relative Step Norm')
      subplot(223)
        semilogy([1:fp_iter],fp_gradnorm,'o-')
        xlabel('Fixed Point Iteration')
        title('Norm of FP Gradient')
      subplot(224)
        semilogy([1:fp_iter],snorm_vec,'o-')
        xlabel('Fixed Point Iteration')
        title('Norm of FP Step')
    figure(3)
      umax = max(U_fp(:));
      imagesc(umax-U_fp), colorbar
      title('Reconstruction')
    figure(4)
      plot([1:nfx]',U_fp(ceil(nfx/2),:), [1:nfx]',f_true(ceil(nfx/2),:))
      title('Cross Section of Reconstruction')
    drawnow
       
  end %for fp

  rel_soln_error = norm(U_fp(:)-f_true(:))/norm(f_true(:))
  clear max_fp_iter;
