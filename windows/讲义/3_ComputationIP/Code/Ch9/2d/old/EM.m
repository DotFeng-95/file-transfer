%  EM.m
%
%  Find approximate solution to Sf = conv(s,f) = d using EM iteration. 
%  EM seeks to minimize the Poisson negative log likelihood function
%    J(f) = sum_i {[Sf]_i - (d_i + sigma^2)*log([Sf]_i + sigma^2)}.

  max_iter = 1000; %%%input(' Max. no. EM iterations = ');
  sigsq = 9; %%%input(' sigma^2 = ');

  dps = d + sigsq;  
  min_d = min(dps(:)+sigsq);
  if (min_d<0)
    fprintf(1,'\n *** Warning: observed image has negative values \n')
    fprintf(1,'              min(d(:) + sigma^2) = %d \n', min_d)
    fprintf(1,'     Compensating... \n\n')
    dps = max(dps,0);
  end;
  f_em = dps + sqrt(eps);
  norm_f = norm(f_true(:));
  em_error = norm(f_em(:)-f_true(:))/norm_f;
  
  
  for k=0:max_iter-1,
    Kfps = mat_prod(S,f_em) + sigsq;
    f_em = f_em.*mat_prod( conj(S),dps./Kfps);
    ek = norm(f_em-f_true)/norm_f;
    em_error = [em_error;ek];
    
    %  Display EM reconstruction
    
    if mod(k,max_iter/20) == 0
      fprintf('   Iteration %i, relative error %6.4e\n',k,ek);
      figure(1)
        subplot(221)
          imagesc(extract(f_true,nx,ny)), colorbar
	  title('True Object')
        subplot(222)
          imagesc(extract(f_em,nx,ny)), colorbar
	  title('EM Reconstruction')
        subplot(223)
          semilogy(em_error)
          xlabel('EM Iteration')
	  title('EM Relative Solution Error')
      colormap(1-gray)
      drawnow
    end
  end;
  
  [mval,mindex] = min(em_error)


