% Min_nn_qp.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initial guess                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  f_0 = ones(size(f_true));        % Initial guess
  %f_0 = f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Cost Function and necessary parameters                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  cost_fn                  = 'ls_function_new';        % Objective fn
  cost_params.cost_fn      = cost_fn;
  cost_params.gam          = input('regularization parameter = ');
  cost_params.S            = S;                                      
  cost_params.d            = d;
  cost_params.n_frames     = n_frames;
  cost_params.Amult        = 'Amult'; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Necessary optimization parameters                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  opt_params.max_iter    = input('Max number of iterations = ');
  opt_params.step_tol    = 1e-15;       % termination criterion norm(step)
  opt_params.grad_tol    = 1e-10; 	% termination criterion norm(grad)
  opt_params.b           = Tmult(cost_params,d);

  out_params             = [];
  
  global TOTAL_COST_EVALS TOTAL_FFTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function call                                                           %
%   Output: x = solution                                                   %
%         histout = iteration history                                      %
%             Each row of histout is                                       %
%   [iteration count, norm(projgrad), J, projgrad norm, step length,       %
%            relative size of active set, linesearch iters, cg iters,      %
%            total cost function evals, total fft's and ifft's]            %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [f,histout] = GPCG_MT(f_0,opt_params,cost_params,out_params);
    %[f,histout] = GPCG_FM(f_0,opt_params,cost_params,out_params);
    %[f,histout] = GPBB(f_0,opt_params,cost_params,out_params);
    %[f,histout] = GP(f_0,opt_params,cost_params,out_params);
    %[f,histout] = PCG(f_0,opt_params,cost_params,out_params);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Display data.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     feval(out_params.output_fn,f,histout);
	
    figure(5)
      imagesc(f), colorbar
      title('Approximate Image')

    figure(6)
      subplot(221)
        xx = histout(1,:);
	yy = histout(3,:); 
        semilogy(xx,yy,'o', xx,yy,'-')
        xlabel('Iteration')
        title('Semilog plot of projected gradient norm')
      subplot(222)
        xx = histout(1,:);
	yy = histout(2,:); 
        plot(xx,yy,'o', xx,yy,'-')        
	title('Function Value')
      subplot(223)
        xx = histout(1,:);
	yy = histout(4,:); 
        semilogy(xx,yy,'o', xx,yy,'-')
        xlabel('Iteration')
        title('Step length reductions')
      subplot(224)
        xx = histout(1,:);
	yy = histout(5,:);
        semilogy(xx,yy,'o', xx,yy,'-')
	title('Relative size of active set')

     fprintf('Total cost function evaluations = %d\n',TOTAL_COST_EVALS)
     fprintf('Total fast fourier transforms = %d\n',TOTAL_FFTS)

  rel_error = norm(f(:)-f_true(:))/norm(f_true(:));
  fprintf(' Relative solution error = %5.2f percent.\n', rel_error*100);
