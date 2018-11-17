%  Landweber.m
%
%  Compute Landweber solution to equation Kf=d.

%  Run "Gen_data.m" first.

  n_iter = 10000; %%%input(' Max. no. of Landweber iterations = ');
  tau = .75; %%%input(' Landweber parameter tau = ');
  error_norm=[];
  it_vec = [];
  Vd = V'*d;
  
  for nu=1:n_iter
    w = 1 - (1-tau*svals.^2).^nu;
    f_nu = V*(Vd.*(w./svals));
  
    %  Compute error and display solution

    power = max(floor(log10(nu))-1,0);
    if mod(nu,10^power) == 0
      error = norm(f_nu-f_true);
      error_norm = [error_norm; error];
      it_vec = [it_vec; nu];
      figure(4)
        plot(x,f_true,'--', x,f_nu,'-')
        xlabel('x axis')
        title('True and Reconstructed Solutions')
      drawnow
    end
    
  end
  
  [e_min,index] = min(error_norm);
  nu_min = it_vec(index(1));
  w = 1 - (1-tau*svals.^2).^nu_min;
  f_min = V*(Vd.*(w./svals));
  nu1 = 3;
  w = 1 - (1-tau*svals.^2).^nu1;
  f_1 = V*(Vd.*(w./svals));
  nu2 = n_iter;
  w = 1 - (1-tau*svals.^2).^nu2;
  f_2 = V*(Vd.*(w./svals));
  
  subplot(221)
    semilogx(it_vec,error_norm,'-o')
    xlabel('Iteration \nu')
    ylabel('||e_{\nu}||')
    title('Landweber Solution Error Norm')
  subplot(222)
    plot(x,f_true,'--', x,f_min,'-')
    xlabel('x axis')
    ylabel('f_{\nu}(x)');
    ttl = ['\nu = ', num2str(nu_min)];
    title(ttl)
  subplot(223)
    plot(x,f_true,'--', x,f_2,'-')
    xlabel('x axis')
    ylabel('f_{\nu}(x)');
    ttl = ['\nu = ', num2str(nu2)];
    title(ttl)
  subplot(224)
    plot(x,f_true,'--', x,f_1,'-')
    xlabel('x axis')
    ylabel('f_{\nu}(x)');
    ttl = ['\nu = ', num2str(nu1)];
    title(ttl)
    