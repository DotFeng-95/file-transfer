%  Primal_dual.m
%
%  Use primal-dual Newton's method to minimize the functional
%      T(u) = ||K*u - d||^2/2 + alpha*J(u),
%  where K is a discretized integral operator, d is discrete data, 
%  ||.|| denotes the l^2 norm, alpha is a positive regularization 
%  parameter, and J is a smooth approximation to the 
%  Total Variation functional.
%      J(u) = .5 * sum_i psi(|[D*u]_i|^2,beta) * Delta_x * Delta_y,
%  where D =[Dx; Dy] is a discretization of the gradient operator 
%  and beta is a positive smoothing parameter. The primal-dual system
%      G1  = K'*(K*u-d) + alpha*(Dx'*Vx + Dy'*Vy) = 0
%      G2x = Dx*u - psi'(|[D*u]_i|^2,beta)*Vx     = 0
%      G2y = Dy*u - psi'(|[D*u]_i|^2,beta)*Vy     = 0
%      max_i |v_i| <= 1
%  is solved using Newton's method. At each iteration, the size of 
%  the V-step is controlled with a line search to ensure that the 
%  constraint max_i sqrt(Vx(i)^2 + Vy(i)^2) < 1 is maintained.

    alpha = input(' Regularization parameter alpha = ');
    beta = input(' TV smoothing parameter beta = ');
    max_pd_iter = input(' No. of primal-dual Newton iterations = ');
    max_cg_iter = input(' No. of CG iterations = ');
    cg_steptol = 1e-7;
    cg_residtol = 1e-4;
    cg_out_flag = 0;  %  If flag = 1, output CG convergence info.

  %  Discretize first derivative operators.
  
  n = nfx;
  nsq = n^2;
  Delta_x = 1 / n;
  Delta_y = Delta_x;
  Delta_xy = Delta_x * Delta_y;
  D = spdiags([-ones(n-1,1) ones(n-1,1)], [0 1], n-1,n) / Delta_x;
  I_trunc = spdiags(ones(n-1,1), 0, n-1,n);
  Dx = kron(D,I_trunc);
  Dy = kron(I_trunc,D);
  
  %  Initialization.

  pd_gradnorm = [];
  snorm_vec = [];
  k_hat_sq = abs(k_hat).^2;
  Kstar_d = integral_op(dat,conj(k_hat),n,n);   %  Compute K'*d.
  f_pd = zeros(n,n);
  u = zeros(n-1,n-1);
  v = zeros(n-1,n-1);
  I = speye((n-1)^2);
  
  for pd_iter = 1:max_pd_iter

    Dxf = Dx * f_pd(:);
    Dyf = Dy * f_pd(:);
    Df_squared = Dxf.^2 + Dyf.^2;
    psi_1 = psi_prime(Df_squared,beta);
    psi_2 = psi_doubleprime(Df_squared,beta);
    Binv = spdiags(psi_1, 0, (n-1)^2,(n-1)^2);

    %  Construct matrix Lbar.
    
%  uu = Binv*Dxf; vv = Binv*Dyf;
  uu = u(:); vv = v(:);
    w = 2 * psi_2 ./ psi_1.^2;
    E11 = spdiags(w.*uu.*Dxf, 0, (n-1)^2,(n-1)^2);
    E12 = spdiags(w.*uu.*Dyf, 0, (n-1)^2,(n-1)^2);
    E21 = spdiags(w.*vv.*Dxf, 0, (n-1)^2,(n-1)^2);
    E22 = spdiags(w.*vv.*Dyf, 0, (n-1)^2,(n-1)^2);
    Lbar = Dx'*Binv*(I + E11)*Dx + Dy' * Binv * (I + E22) * Dy ...
 	+ Dx'*Binv*E12*Dy + Dy'*Binv*E21*Dx;
    Lbar = (Lbar + Lbar')/2;

    %  Solve primal-dual system.
    
    KstarKf = integral_op(f_pd,k_hat_sq,n,n);
    r = Kstar_d(:) - KstarKf(:) - alpha*Dx'*Binv*Dxf - alpha*Dy'*Binv*Dyf;
    
    gradnorm = norm(r);
    pd_gradnorm = [pd_gradnorm; gradnorm];
    
    %  Use PCG iteration to solve linear system
    % (K'*K + alpha*Lbar)*Delta_f = r
    
    Delta_f = zeros(n,n);
    resid = reshape(r,n,n);
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

      %  Form product (K'*K + alpha*Lbar)*ph.

      KstarKp = integral_op(ph,k_hat_sq,n,n);
      Ahph = KstarKp + alpha * reshape(Lbar*ph(:),nfx,nfy);

      %  Update Delta_f and residual.
    
      alphak = rd / (ph(:)'*Ahph(:));
      Delta_f = Delta_f + alphak*ph;
      resid = resid - alphak*Ahph;
      rdlast = rd;
      
      residnorm = norm(resid(:));
      stepnorm = abs(alphak)*norm(ph(:))/norm(Delta_f(:));
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
           cgiter, residnormvec(cgiter), stepnormvec(cgiter))
        figure(1)
          subplot(221)
            semilogy(residnormvec/residnormvec(1),'o')
            title('CG Residual Norm')
          subplot(222)
            semilogy(stepnormvec,'o')
            title('CG Step Norm')
      else
        if cgiter == 1
          fprintf(' ... computing Delta_f via CG iterations ...\n');
        end
      end
    end %end for CG
keyboard
    
    %  Compute dual Newton steps.
    
    uu = u(:); vv = v(:);
    E11 = spdiags(w.*uu.*Dxf, 0, (n-1)^2,(n-1)^2);
    E12 = spdiags(w.*uu.*Dyf, 0, (n-1)^2,(n-1)^2);
    E21 = spdiags(w.*vv.*Dxf, 0, (n-1)^2,(n-1)^2);
    E22 = spdiags(w.*vv.*Dyf, 0, (n-1)^2,(n-1)^2);
    Delf = Delta_f(:);
    Delta_u = -u(:) + Binv*(Dx*(f_pd(:)+Delf) + E11*Dx*Delf + E12*Dy*Delf);
    Delta_u = reshape(Delta_u,n-1,n-1);
    Delta_v = -v(:) + Binv*(Dy*(f_pd(:)+Delf) + E21*Dx*Delf + E22*Dy*Delf);
    Delta_v = reshape(Delta_v,n-1,n-1);
    
    %  Update primal variable.
    
    f_pd = f_pd + Delta_f;
    
    %  Perform line search in the dual variables u,v.  Calculate
    %  minimum rho in (0,1] for which |V(i,j) + rho*Delta_V(i,j)| = 1 
    %  for all (i,j). Here |.| denotes Euclidean norm on R^2.
      
    VdotdV = u.*Delta_u + v.*Delta_v;
    absVsq = u.^2 + v.^2;
    absdVsq = Delta_u.^2 + Delta_v.^2;
    rho = (sqrt(VdotdV.^2 + absdVsq.*(1-absVsq)) - VdotdV) ./ (absdVsq + eps);
    rho_min = min(min(rho(:)),1);
    if rho_min < 1
      rho_min = .9*rho_min;
    end
    u = u + rho_min*Delta_u;
    v = v + rho_min*Delta_v;
      
    %  Display results.
    
    snorm = sqrt(norm(Delta_f(:))^2 + norm(Delta_u(:))^2 + norm(Delta_v(:))^2);
    snorm_vec = [snorm_vec; snorm];
    
    %   Output primal-dual Newton convergence information.
   
    fprintf(' PD iter=%3.0f, ||grad||=%6.4e, ||step||=%6.4e, nCG=%3.0f\n', ...
       pd_iter, gradnorm, snorm, cgiter);
   
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
        semilogy([1:pd_iter],pd_gradnorm,'o-')
        xlabel('Primal-Dual Newton Iteration')
        title('Norm of PD Gradient')
      subplot(224)
        semilogy([1:pd_iter],snorm_vec,'o-')
        xlabel('Primal-Dual Newton Iteration')
        title('Norm of PD Step')
       
     %  Leave out corner entry f_pd(n,n), since it is coupled to other
     %  entries only through K'*K, and not through Lbar.
     
     U_PD = f_pd(1:n-1,1:n-1);
     figure(3)
       umax = max(U_PD(:));
       imagesc(umax-U_PD), colorbar
       title('Reconstruction')
     figure(4)
       plot([1:n-1]',U_PD(ceil(nfx/2),:), [1:n-1]',f_true(ceil(nfx/2),1:n-1))
       title('Cross Section of Reconstruction')
     drawnow
       
  end %for pd iteration

  rel_soln_error = norm(f_pd(:)-f_true(:))/norm(f_true(:))
  clear max_pd_iter;







