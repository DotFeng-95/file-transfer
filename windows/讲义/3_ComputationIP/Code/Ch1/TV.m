%  TV.m
%
%  Compute TV regularized solutions for the 1-D test problem.

  n_alpha = 40;
  alpha_vec = logspace(-7,-2,n_alpha);
  beta = .01;
  enorm = zeros(n_alpha,1);
  f_mat = [];
  dx = 1/n; 
  
  for ii = 1:n_alpha
    alpha = alpha_vec(ii);
    [f_alpha,gnorm] = primal_dual(alpha,beta,20,K,d,0);
  
    %  Display solution and gradient norm
 
    figure(4)
      plot(x,f_true,'--', x,f_alpha,'-')
      xlabel('x axis')
      title('True and TV Regularized Solutions')
    figure(5)
      semilogy(gnorm,'o-')
      xlabel('Quasi-Newton Iteration')
      title('G Norm')
    drawnow
      
    f_mat = [f_mat f_alpha];
    enorm(ii) = norm(f_alpha - f_true);
  end
  
  [alpha_min,indx] = min(enorm);

 figure(1)
  subplot(221)
    semilogx(alpha_vec,enorm,'o-')
    xlabel('\alpha')
    ylabel('||e_{\alpha}||')
    title('Solution Error Norm')
    xlim([alpha_vec(1) alpha_vec(n_alpha)])
  subplot(222)
    fmin = f_mat(:,indx(1));
    plot(x,fmin,'-', x,f_true,'--')
    xlabel('x axis')
    ylabel('f_{\alpha}(x)')
    ttl = ['\alpha = ', num2str(alpha_vec(indx(1)))];
    title(ttl)
  subplot(223)
    f1 = f_mat(:,1);
    plot(x,f1,'-', x,f_true,'--')
    xlabel('x axis')
    ylabel('f_{\alpha}(x)')
    ttl = ['\alpha = ', num2str(alpha_vec(1))];
    title(ttl)
  subplot(224)
    f2 = f_mat(:,n_alpha);
    plot(x,f2,'-', x,f_true,'--')
    xlabel('x axis')
    ylabel('f_{\alpha}(x)')
    ttl = ['\alpha = ', num2str(alpha_vec(n_alpha))];
    title(ttl)
    
