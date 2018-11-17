%  Primal_dual.m
%
%  Use primal-dual Newton's method to minimize the functional
%      T(u) = ||K*u - d||^2/2 + alpha*J(u),
%  where K is a discretized integral operator, d is discrete data, 
%  ||.|| denotes the l^2 norm, alpha is a positive regularization 
%  parameter, and J is a smooth approximation to the 
%  Total Variation functional.
%      J(u) = sum_i 2*psi(|[D*u]_i|^2,beta) * Delta_x,
%  where D is a discretization of the first derivative operator and
%  beta is a positive smoothing parameter.
%  The primal-dual system
%      g1 = K'*(Ku-d) + alpha*D'*v         = 0
%      g2 = D*u - psi'(|[D*u]_i|^2,beta)*v = 0
%      max_i |v_i| <= 1
%  is solved using Newton's method. At each iteration, the size of 
%  the v-step is controlled with a line search to ensure that the 
%  constraint max_i |v_i| <= 1 is maintained.

  if ~exist('max_pd_iter')
    alpha = input(' Regularization parameter alpha = ');
    beta = input(' TV smoothing parameter beta = ');
    max_pd_iter = input(' No. of primal-dual Newton iterations = ');
    max_cg_iter = input(' No. of CG iterations = ');
    cg_steptol = 1e-4;
    cg_residtol = 1e-4;
  end

  %  Set up discretization of first derivative operators.
  
  n = nfx;
  nsq = n^2;
  Delta_x = 1 / n;
  Delta_y = Delta_x;
  D = spdiags([-ones(n-1,1) ones(n-1,1)], [0 1], n-1,n) / Delta_x;
  I_trunc = spdiags(ones(n-1,1), 0, n-1,n);
  Dx = kron(D,I_trunc);
  Dy = kron(I_trunc,D);
  
%  I_trunc2 = spdiags(ones(n-1,1), 1, n-1,n);
%  Dx2 = kron(D,I_trunc2);
%  Dy2 = kron(I_trunc2,D);

  %  Initialization.

  pd_gradnorm = [];
  snorm_vec = [];
  k_hat_sq = abs(k_hat).^2;
  b = integral_op(d1,conj(k_hat));   %  Compute b = K'*d.
  b  = b(1:nfx,1:nfy);
  U_pd = zeros(n,n);
  Vx = zeros(n-1,n-1);
  Vy = zeros(n-1,n-1);
  
  for pd_iter = 1:max_pd_iter

    Dxu = Dx * U_pd(:);
    Dyu = Dy * U_pd(:);
    psi_1 = 1./ sqrt((Dxu).^2 + (Dyu).^2 + beta^2);
    psi_2 = -.5 * ((Dxu).^2 + (Dyu).^2 + beta^2).^(-1.5);
    B = spdiags(psi_1, 0, (n-1)^2,(n-1)^2);
    
    %  Compute gradient.
    
    KstarKu = integral_op(U_pd,k_hat_sq);
    Div_v = reshape(Dx' * Vx(:) + Dy' * Vy(:), n,n);
    G1 = KstarKu(1:n,1:n) - b + alpha * Div_v;
    G2x = reshape(Dxu,n-1,n-1) - Vx ./ reshape(psi_1,n-1,n-1);
    G2y = reshape(Dyu,n-1,n-1) - Vy ./ reshape(psi_1,n-1,n-1);
    gradnorm = sqrt(norm(G1(:))^2 + norm(G2x(:))^2 + norm(G2y(:))^2);
    pd_gradnorm = [pd_gradnorm; gradnorm];
    
    %  Set up regularization operator L. First compute symmetrized M.
    
    n2 = (n-1)^2;
    B = spdiags(psi_1, 0, n2,n2);
    fpsi = 2*psi_2 ./ psi_1.^2;
    M11 = spdiags(1 + fpsi .* Dxu .* Vx(:), 0, n2,n2);
    M22 = spdiags(1 + fpsi .* Dyu .* Vy(:), 0, n2,n2);
    m12 = fpsi .* Dyu .* Vx(:);
    m21 = fpsi .* Dxu .* Vy(:);
    m_average = (m12 + m21) / 2;
    Msym = spdiags(m_average, 0, n2,n2);
    M12 = spdiags(m12, 0, n2,n2);
    M21 = spdiags(m21, 0, n2,n2);
    
    L = Dx'*B*(M11*Dx + Msym*Dy) + Dy'*B*(Msym*Dx + M22*Dy);
  
    %  Compute Delta_U using PCG iteration.
    
%  CG initialization.
  
    Delta_U = zeros(n,n);
    resid = -G1 - alpha * reshape(Dx' * B * G2x(:) + Dy' * B * G2y(:),n,n);
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

      KstarKp = integral_op(ph,k_hat_sq);
      Ahph = KstarKp(1:nfx,1:nfy) + alpha * reshape(L*ph(:),nfx,nfy);

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
      
      fprintf('   CG iter%3.0f, ||resid||=%6.4e, ||step||=%6.4e \n', ... 
         cgiter, residnormvec(cgiter), stepnormvec(cgiter))
      figure(1)
        subplot(221)
          semilogy(residnormvec/residnormvec(1),'o')
          title('CG Residual Norm')
        subplot(222)
          semilogy(stepnormvec,'o')
          title('CG Step Norm')
    end %end for CG

    Delta_Vx = B * (G2x(:) + M11*Dx*Delta_U(:) + M12*Dy*Delta_U(:));
    Delta_Vy = B * (G2y(:) + M21*Dx*Delta_U(:) + M22*Dy*Delta_U(:));
    Delta_Vx = reshape(Delta_Vx,n-1,n-1);
    Delta_Vy = reshape(Delta_Vy,n-1,n-1);
    
    %  Update primal variable
    
    U_pd = U_pd + Delta_U;
    
    %  Perform line search in the dual variable v.  Calculate
    %  minimum rho in (0,1] for which |v(i) + rho*Delta_v(i)| = 1.
      
    VdotdV = Vx.*Delta_Vx + Vy.*Delta_Vy;
    absVsq = Vx.^2 + Vy.^2;
    absdVsq = Delta_Vx.^2 + Delta_Vy.^2;
    rho = (sqrt(VdotdV.^2 + absdVsq.*(1-absVsq)) - VdotdV) ./ (absdVsq + eps);
    rho_min = min(min(rho(:)),1);
    if rho_min < 1
      rho_min = .9*rho_min;
    end
    Vx = Vx + rho_min*Delta_Vx;
    Vy = Vy + rho_min*Delta_Vy;
      
    %  Display results.
    
    snorm = sqrt(norm(Delta_U)^2 + norm(Delta_Vx)^2 + norm(Delta_Vy)^2);
    snorm_vec = [snorm_vec; snorm];
    
    %   Output primal-dual Newton convergence information.
   
    fprintf(' PD iter = %3.0f, ||gradient|| = %6.4e,  ||step|| = %6.4e.\n', ...
       pd_iter, gradnorm, snorm);
   
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
     drawnow
       
     %  Leave out corner entry U_pd(n,n), since it is coupled to other
     %  entries only through K'*K, and not through L.
     
     U_PD = U_pd(1:n-1,1:n-1);
     figure(3)
       umax = max(U_PD(:));
       imagesc(umax-U_PD), colorbar
       title('Reconstruction')
     figure(4)
       plot([1:n-1]',U_PD(ceil(nfx/2),:), [1:n-1]',f_true(ceil(nfx/2),1:n-1))
       title('Cross Section of Reconstruction')
       
  end %for pd iteration

  rel_soln_error = norm(U_pd(:)-f_true(:))/norm(f_true(:))
  clear max_pd_iter;







