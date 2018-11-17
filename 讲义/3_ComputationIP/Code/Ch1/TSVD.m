%  TSVD.m
%
%  Compute TSVD solution to equation Kf=d.

%  Run "Gen_data.m" first.


  alpha_vec = logspace(-3,0,40);
  n_alpha = max(size(alpha_vec));
  error_norm=zeros(n_alpha,1);
  
  for i=1:n_alpha
    alpha = alpha_vec(i);
    w = (svals > alpha);
    f_alpha = V*((V'*d).*(w./svals)); %%%V*diag(w./svals)*(V'*d);
  
  %  Display solution

    figure(4)
      plot(x,f_true,'--', x,f_alpha,'-')
      xlabel('x axis')
    header = ['True and Reconstructed Solutions for \alpha = ',num2str(alpha)];
      title(header)
    disp(' Hit any key to continue.'); pause
    
    error_norm(i) = norm(f_alpha-f_true);
  end
  
  [e_min,index] = min(error_norm);
  alpha_min = alpha_vec(index(1));
  w = (svals > alpha_min);
  f_min = V*((V'*d).*(w./svals));
  n1 = 1;
  w = (svals > alpha_vec(n1));
  f_1 = V*((V'*d).*(w./svals));
  n2 = n_alpha-1;
  w = (svals > alpha_vec(n2));
  f_2 = V*((V'*d).*(w./svals));
  
  subplot(221)
    semilogx(alpha_vec,error_norm,'-')
    xlabel('\alpha')
    ylabel('||e_{\alpha}||')
    title('Norm of TSVD Solution Error')
    xlim([alpha_vec(1) alpha_vec(n_alpha)]);
  subplot(222)
    plot(x,f_true,'--', x,f_min,'-')
    xlabel('x axis')
    ylabel('f_{\alpha}(x)');
    ttl = ['\alpha = ', num2str(alpha_min)];
    title(ttl)
  subplot(223)
    plot(x,f_true,'--', x,f_1,'-')
    xlabel('x axis')
    ylabel('f_{\alpha}(x)');
    ttl = ['\alpha = ', num2str(alpha_vec(n1))];
    title(ttl)
  subplot(224)
    plot(x,f_true,'--', x,f_2,'-')
    xlabel('x axis')
    ylabel('f_{\alpha}(x)');
    ttl = ['\alpha = ', num2str(alpha_vec(n2))];
    title(ttl)
    
