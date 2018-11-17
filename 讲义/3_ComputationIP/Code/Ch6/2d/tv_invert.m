%  tv_invert.m   
%
%  (Tikhonov) regularized least squares solution of parameter estimation
%  problem u(q) = z. Compute q to minimize
%    ||u(q)-z||^2/2 + alpha*J_TV(q).
%  This requires solution to a sequence of quadratic problems
%    min_dq ||K*dq + r(q)||^2/2 + alpha * <L(q)(q+dq),q+dq>,
%  where r(q) = u(q)-z.
%  The matrix L is a discretization of a diffusion operator with no-flux
%  boundary conditions. It has a 1-dimensional null space, spanned by vector 
%  v0 = ones(n,1). Set dq = a*v0 + dq_perp, where dq_perp is the
%  projection of dq onto the orthogonal complement of null(L). This yields
%  the block system
%
%     (v0'*w0) * a  +       w0' * dq_perp       = -v0'*g
%        w0 * a     +  (H_ls + alpha*L)*dq_perp = -g
%
%  where the gradient g = dudq'*r(q), the (Gauss-Newton approximation to
%  the) least squares Hessian H_ls = dudq'*dudq, and 
%  w0 = H_ls * v0. Convert to Schur complement form
%
%     (v0'*w0) * a  +       w0' * dq_perp       = -v0'*g
%     [H_ls + alpha*L - w0*w0'/(v0'*w0)]*dq_perp = -g + (v0'*g)/(v0'*w0)*w0

%  Run "setup.m" first.

  disp(' *** Be sure to run setup.m first! ***');
  alpha = input(' Regularization parameter alpha = ');
  beta = 1; %%%input(' TV regularity parameter beta = ');
  Newton_max = input(' Max no. of QuasiNewton iterations = ');
  pcg_max = N/2; %input(' Max no. of PCG iterations = ');
  pcg_tol = 1e-3; %%%input(' PCG relative error stopping tolerance = ');
  pcgioflag = 1;  %  Output info inside function pcg.m
  pcgflag = input(' PCG_flag (1 if PCG; else 0): ');
  if pcgflag == 1,
    Linv_iter = 20;
    Linv_tol = 1e-6;
    nu = 2;           %  No. of Multigrid "smoother" sweeps.
    level = log2(nx);  %  No. of Multigrid levels.
  end

%  Initialization.

  reset = input(' Enter 1 to reset initial guess; else enter 0: ');
  if reset == 1
    q_vec = zeros(N,1);
  end
  q_mat = reshape(q_vec,nx,ny);
  gradnormvec = [];
  PCG_conv_factor = [];
  v0 = ones(N,1) / sqrt(N);

%  Perform Quasi-Newton iterations. Keep track of the norm of the gradient 
  
  for iter = 1:Newton_max

    %  Set up diffusion coefficient for operator alpha*L(u).
  
    kappa = get_kappa(u_mat,beta);
    alphaL = alpha * get_Lmat(kappa);

    %  Compute (approximate) Hessian and gradient of objective functional.

    [u_mat,A0,grad_ls] = co_lsgrad(q_mat,b,z,Cindex);
    g = grad_ls(:) + alphaL * q_vec;
    gradnormvec(iter) = norm(g);

    %  Compute projections onto null(L)^{perp}.

    w0 = co_lsHess(reshape(v0,nx,ny),q_mat,u_mat,A0,Cindex);
    c = v0'*w0(:);
    gperp = -g + (v0'*g/c) * w0(:);

    %  Solve Shur system using PCG. The preconditioner takes the form
    %     P_perp = alpha * L_perp,
    %  where L_perp is the restriction of L to null(L)^perp.

    %  PCG initialization.

    residnormvec = [];
    pcgiter = 0;
    residrat = 1;
    dqperp = zeros(N,1);
    resid = gperp;

    while pcgiter < pcg_max & residrat > pcg_tol
      pcgiter = pcgiter + 1;

      residnormvec(pcgiter) = norm(resid);

      %  Solve P_r * d = resid, where the preconditioning matrix 
      %  P_r = alpha*L restricted to null(L)-perp.
      
      if pcgflag == 1,
        [d,tmp] = Linvert(zeros(N,1),alpha*kappa, ...
           resid,Linv_iter,Linv_tol,level,nu);
      else   
        d = resid; %%% No preconditioning
      end
      rd = resid'*d;
      if pcgiter == 1,
         p = d; 
       else
         betak = rd / rdlast;
         p = d + betak * p;
      end

      Hls_p = co_lsHess(reshape(p,nx,ny),q_mat,u_mat,A0,Cindex);
      Hp = Hls_p(:) + alphaL * p - ((w0(:)'*p)/c) * w0(:);
      alphak = rd / (p'*Hp);
      dqperp = dqperp + alphak*p;
      resid = resid - alphak*Hp;
      rdlast = rd;
      residrat = residnormvec(pcgiter)/residnormvec(1);

      if pcgioflag == 1,
        fprintf('  PCG iter %3.0f, ||resid|| = %6.4e \n', ...
          pcgiter, residnormvec(pcgiter));
      end

    end  %  end of PCG iteration.

    a = (-v0'*g - w0(:)'*dqperp)/c;
    dq = a*v0 + dqperp;

    q_vec = q_vec + dq;
    q_mat = reshape(q_vec,nx,ny);

    %  Output PCG information and approximate minimizer.

    PCG_conv_factor(iter) = -mean(diff(log(residnormvec)));
    npcg = max(size(residnormvec));

    fprintf('\n Q-N iter %2.0f, ||grad||= %6.4e, PCG iter= %3.0f,', ...
        iter, gradnormvec(iter), npcg);
    fprintf(' PCG conv factor= %6.4e\n\n', PCG_conv_factor(iter));

    figure(2)
      clg
      subplot(221)
      imagesc(exp(q_mat)), colorbar
      title('diffusivity kappa')
    subplot(222)
      semilogy(residnormvec,'o')
      xlabel('PCG iteration')
      title('Norm of PCG residual')
    drawnow

  end  %  end of Quasi-Newton iteration.

%  Quasi-Newton information.

  figure(3)
    clg
    subplot(221)
      imagesc(exp(q_mat)), colorbar
      title('diffusivity kappa')
    subplot(222)
      semilogy(gradnormvec,'o')
      xlabel('Quasi-Newton iteration')
      title('Norm of gradient')
    subplot(223)
      plot(PCG_conv_factor,'o')
      xlabel('Quasi-Newton iteration')
      title('PCG convergence factor')
    drawnow

%  Compute and plot spectrum.

if 0==1
  disp(' ... computing spectra of operators ...');

  dudq = deriv(q_mat,b);
  H_ls = dudq'*dudq;
  alphaL = alpha * get_Lmat(kappa);
  H = H_ls + alphaL;
  spec_H = sort(real(eig(H)));
  Hperp = H - w0(:)*w0(:)'/c;
  CinvHperp = pinv(full(alphaL)) * Hperp;
  spec_PCG = sort(real(eig(CinvHperp)));
  
  figure(4)
  subplot(221)
    semilogy(spec_H,'o')
    xlabel('index i')
    ylabel('lambda_i')
    title('Spectrum of H')
  subplot(223)
    semilogy(spec_PCG(2:N),'o')
    xlabel('index i')
    ylabel('lambda_i')
    title('Spectrum of H_PCG')

end
