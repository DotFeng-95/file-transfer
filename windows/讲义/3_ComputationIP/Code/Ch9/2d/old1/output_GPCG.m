  function [iter_hist,flag] = ...
      output_GPCG(iter,ls_iter,cg_iter,x_old,x,g,J,Active,iter_hist,params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Output and storage of iteration history.                                 %
%  Check stopping criteria.                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  global TOTAL_COST_EVALS TOTAL_FFTS

  s = x - x_old;
  stepnorm = norm(s(:));
  pg = x - max(x - g,0);
  gradnorm = norm(pg(:));
  ia = sum(Active(:));
  ndim = length(Active(:));

  iter_hist = [iter_hist [iter; J; gradnorm; stepnorm; ia/ndim; ...
                      TOTAL_COST_EVALS;TOTAL_FFTS]];

  fprintf('It=%d cost=%6.3e |s|=%6.3e |pg|=%6.3e ls_it=%d cg_it=%d #Active=%d\n',iter, J, stepnorm, gradnorm, ls_iter, cg_iter, ia); 
  
        %  Check stopping criteria.
    
    if iter >= params.max_iter
        flag = 1;
	fprintf('   *** Max iterations exceeded ***\n');
    elseif iter_hist(4,iter+1) < params.step_tol%/ norm(x(:)) < params.step_tol
        flag = 2;
	fprintf('   *** Min step norm tolerance met ***\n');
    elseif iter_hist(3,iter+1) / iter_hist(3,1) < params.grad_tol
        flag = 3;
	fprintf('   *** Min gradient norm tolerance met ***\n');
    else
        flag = 0;
    end

