%  Figs.m
%
%  Generate plots for Figs 5.1 and 5.2. Fig. 5.1 compares the performance 
%  of the steepest descent method vs. primal Newton's method on a 1-D 
%  test problem. Fig. 5.2 compares lagged diffusivity fixed point iteration
%  vs. primal-dual Newton's method.

  disp(' Be sure to run Setup.m and Gen_ualpha.m first.')
  load u_alpha;   %  Get "true" reconstruction

%  Regularization parameters.

  alpha = 1e-4;
  beta = .1;
  
%  Fig. 4.1

  SD_max = 100;    %  Max iteration count for steepest descent
  Newt_max = 100;  %  Max iteration count for primal Newton's method
  Newt_grad_tol = 1e-10;
  Steepest_d;
  Newton;
  
%  Compute solution, or "estimation", errors from convergence histories.

  [dum,n_iter] = size(SD_conv_history);
  SD_enorm = [];
  for k = 1:n_iter
    uk = SD_conv_history(:,k);
    SD_enorm = [SD_enorm; norm(uk - u_alpha)];
  end
  
  [dum,n_iter] = size(Newton_conv_history);
  Newton_enorm = [];
  for k = 1:n_iter
    uk = Newton_conv_history(:,k);
    Newton_enorm = [Newton_enorm; norm(uk - u_alpha)];
  end
  
  figure(3)
    ind1 = [1:max(size(sd_gradnorm))]';
    ind2 = [1:max(size(newt_gradnorm))]';
    semilogy(ind1,sd_gradnorm,'o-', ind2,newt_gradnorm,'*-')
    xlabel('Iteration')
    title('Norm of Gradient for Steepest Descent (o) and Newton Method (*)')

  figure(4)
    ind1 = [1:max(size(SD_enorm))]';
    ind2 = [1:max(size(Newton_enorm))]';
    semilogy(ind1,SD_enorm,'o-', ind2,Newton_enorm,'*-')
    xlabel('Iteration')
    title('Norm of Soln Error for Steepest Descent (o) and Primal Newton (*)')
      
%  Fig. 4.2

  fp_iter = 20;
  grad_tol = 1e-2;
  PDnewt_iter = 20;
  Fixed_pt;
  Primal_dual;

%  Compute solution, or "estimation", errors from convergence histories.

  [dum,n_iter] = size(FPconv_history);
  FPenorm = [];
  for k = 1:n_iter
    uk = FPconv_history(:,k);
    FPenorm = [FPenorm; norm(uk - u_alpha)];
  end
  [dum,n_iter] = size(PDconv_history);
  PDenorm = [];
  for k = 1:n_iter
    uk = PDconv_history(:,k);
    PDenorm = [PDenorm; norm(uk - u_alpha)];
  end
  
  figure(5)
    ind1 = [1:max(size(fp_gradnorm))]';
    ind2 = [1:max(size(pd_gradnorm))]';
    semilogy(ind1,fp_gradnorm,'o-', ind2,pd_gradnorm,'*-')
    xlabel('Iteration')
    title('Norm of Gradient for Fixed Point (o) and PD Newton Method (*)')

  figure(6)
    ind1 = [1:max(size(FPenorm))]';
    ind2 = [1:max(size(PDenorm))]';
    semilogy(ind1,FPenorm,'o-', ind2,PDenorm,'*-')
    xlabel('Iteration')
    title('Norm of Soln Error for Fixed Point (o) and PD Newton Method (*)')





