%  Primal_dual.m
%
%  Use primal-dual Newton's method to minimize the functional
%      T(u) = ||K*u - d||^2/2 + alpha*J(u),
%  where K is a discretized integral operator, d is discrete data, 
%  ||.|| denotes the l^2 norm, alpha is a positive regularization 
%  parameter, and J is a smooth approximation to the 
%  Total Variation functional.
%      J(u) = sum_i 2*psi(|[D*u]_i|^2,beta) * Delta_x * Delta_y,
%  where D =[Dx; Dy] is a discretization of the gradient operator 
%  and beta is a positive smoothing parameter. The primal-dual system
%      G1  = K'*(K*u-d) + alpha*(Dx'*Vx + Dy'*Vy) = 0
%      G2x = Dx*u - psi'(|[D*u]_i|^2,beta)*Vx     = 0
%      G2y = Dy*u - psi'(|[D*u]_i|^2,beta)*Vy     = 0
%      max_i |v_i| <= 1
%  is solved using Newton's method. At each iteration, the size of 
%  the V-step is controlled with a line search to ensure that the 
%  constraint max_i sqrt(Vx(i)^2 + Vy(i)^2) < 1 is maintained.

  if ~exist('max_pd_iter')
    alpha = input(' Regularization parameter alpha = ');
    beta = input(' TV smoothing parameter beta = ');
    max_pd_iter = input(' No. of primal-dual Newton iterations = ');
    max_cg_iter = input(' No. of CG iterations = ');
    cg_steptol = 1e-7;
    cg_residtol = 1e-4;
    cg_out_flag = 0;  %  If flag = 1, output CG convergence info.
  end

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
  
%  I_trunc2 = spdiags(ones(n-1,1), 1, n-1,n);
%  Dx2 = kron(D,I_trunc2);
%  Dy2 = kron(I_trunc2,D);

  %  Initialization.

  pd_gradnorm = [];
  G1norm = [];
  snorm_vec = [];
  k_hat_sq = abs(k_hat).^2;
  Kstar_d = integral_op(dat,conj(k_hat),n,n);   %  Compute K'*d.
  U_pd = zeros(n,n);
  Vx = zeros(n-1,n-1);
  Vy = zeros(n-1,n-1);
  
  for pd_iter = 1:max_pd_iter

    Dxu = Dx * U_pd(:);
    Dyu = Dy * U_pd(:);
    Du_squared = Dxu.^2 + Dyu.^2;
    psi_1 = psi_prime(Du_squared,beta);
    psi_2 = psi_doubleprime(Du_squared,beta);
    B = spdiags(psi_1, 0, (n-1)^2,(n-1)^2);
    
    %  Compute components G1, G2 of primal-dual system.
    
    KstarKu = integral_op(U_pd,k_hat_sq,n,n);
    Div_V = reshape(Dx' * Vx(:) + Dy' * Vy(:), n,n);
    G1 = KstarKu - Kstar_d + alpha * Delta_xy * Div_V;
    G2x = reshape(Dxu,n-1,n-1) - Vx ./ reshape(psi_1,n-1,n-1);
    G2y = reshape(Dyu,n-1,n-1) - Vy ./ reshape(psi_1,n-1,n-1);
    gradnorm = sqrt(norm(G1(:))^2 + norm(G2x(:))^2 + norm(G2y(:))^2);
    pd_gradnorm = [pd_gradnorm; gradnorm];
    G1norm = [G1norm; norm(G1(:))];
    
    %  Set up regularization operator L. First compute blocks of 
    %     [1 + c*Dxu*Vx    c*Dyu*Vx  ]
    %     [  c*Dxu*Vy    1 + c*Dyu*Vy]
    %  where c = 2*psi''/(psi')^2.
    
    n2 = (n-1)^2;
    B = spdiags(psi_1, 0, n2,n2);
    c_psi = 2*psi_2 ./ psi_1.^2;
    M11 = spdiags(1 + c_psi .* Dxu .* Vx(:), 0, n2,n2);
    M22 = spdiags(1 + c_psi .* Dyu .* Vy(:), 0, n2,n2);
    m12 = c_psi .* Dyu .* Vx(:);
    m21 = c_psi .* Dxu .* Vy(:);
    m_average = (m12 + m21) / 2;
    Msym = spdiags(m_average, 0, n2,n2);
    M12 = spdiags(m12, 0, n2,n2);
    M21 = spdiags(m21, 0, n2,n2);
    L = Dx'*B*(M11*Dx + Msym*Dy) + Dy'*B*(Msym*Dx + M22*Dy);
    L = L * Delta_xy;
  
    %  Use PCG iteration to solve linear system
    % (K'*K + alpha*L)*Delta_u = -K'*(K*u-d) - alpha*(Dx'*B*Dxu + Dy'*B*Dyu)
    
    Delta_U = zeros(n,n);
    resid = Kstar_d - KstarKu - alpha * Delta_xy * reshape(Dx'*B*Dxu(:) ...
	+ Dy'*B*Dyu(:),n,n);
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

      %  Form product (K'*K + alpha*L)*ph.

      KstarKp = integral_op(ph,k_hat_sq,n,n);
      Ahph = KstarKp + alpha * reshape(L*ph(:),nfx,nfy);

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
          fprintf(' ... computing Delta_u via CG iterations ...\n');
        end
      end

    end %end for CG

    Delta_Vx = -Vx(:) + B * (Dxu(:) + M11*Dx*Delta_U(:) + M12*Dy*Delta_U(:));
    Delta_Vy = -Vy(:) + B * (Dyu(:) + M21*Dx*Delta_U(:) + M22*Dy*Delta_U(:));
    Delta_Vx = reshape(Delta_Vx,n-1,n-1);
    Delta_Vy = reshape(Delta_Vy,n-1,n-1);
    
    %  Update primal variable
    
    U_pd = U_pd + Delta_U;
    
    %  Perform line search in the dual variables Vx,Vy.  Calculate
    %  minimum rho in (0,1] for which |V(i,j) + rho*Delta_V(i,j)| = 1 
    %  for all (i,j). Here |.| denotes Euclidean norm on R^2.
      
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
     drawnow
       
  end %for pd iteration

  rel_soln_error = norm(U_pd(:)-f_true(:))/norm(f_true(:))
  clear max_pd_iter;







