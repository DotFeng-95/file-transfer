  function [x_c,iter_hist] = ...
      GPCG_MT(x_c,opt_params,cost_params,out_params);

%  This m-file was written by John Bardsley, Dept. of Mathematical Sciences,
%  Montana State University, in September 2001.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Numerical parameters and function names.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  max_iter     = opt_params.max_iter;    %  Max. no. iterations.
  cost_fn      = cost_params.cost_fn;    %  m-file name for cost fn eval.
  iter = 0;
  total_iters = 0;
  flag = 0;
  global TOTAL_COST_EVALS TOTAL_FFTS
  TOTAL_FFTS = 0;
  TOTAL_COST_EVALS = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Cost function and gradient evaluation. Active set determination.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  x_c = max(x_c,0);               % Put initial guess in feasible set
  Active = (x_c == 0);            % Compute active set
  ia = sum(Active(:)); 			% # of active indices
  [J_c,g_c] = feval(cost_fn,x_c,cost_params); 
  pg_c = proj_grad(g_c,x_c);  % Projected gradient
  npg_c = norm(pg_c(:)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Storage vector for important numerical information                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %  Column k of iter_hist contains the following information:
  %    iter_hist(1,k) = Outer iteration count
  %    iter_hist(2,k) = cost function value
  %    iter_hist(3,k) = projected gradient norm
  %    iter_hist(4,k) = step norm
  %    iter_hist(5,k) = relative size of active set
  %    iter_hist(6,k) = total cost function evals 
  %    iter_hist(7,k) = total fft's and ifft's

  iter_hist = [iter; J_c; npg_c; 0; ia/length(x_c(:)); ...
      TOTAL_COST_EVALS;TOTAL_FFTS];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Outer iteration.                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  while flag == 0;

    current_cost_evals = TOTAL_COST_EVALS;
    ls_flag = 0;
    ls_iter = 0;
    J_diff_max = 0;
    iter = iter + 1;
    y_c = x_c;
    x_p = x_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Projected gradient iteration. Identify constraints!!!!                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  while ls_flag == 0
    
    ls_iter = ls_iter + 1;

    d_c = -g_c.*((1 - Active) + Active.*(g_c < 0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  Linesearch                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [y_new,J_new,g_new,ls_param,opt_flag] = ...
          linesearch_MT(y_c,d_c,g_c,J_c,cost_fn,cost_params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check stopping criteria and update information.                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      J_diff = J_c - J_new;
      J_diff_max = max(J_diff,J_diff_max);
      if J_diff <= .25*J_diff_max         %  Insufficient decrease in J.
	ls_flag = 1;
      else
        Active_new = (y_new == 0);
        if  Active_new == Active          %  Active set unchanged.
	  ls_flag = 1;
        end
      end
      
      Active = Active_new;
      y_c = y_new;
      J_c = J_new;
      g_c = g_new;
      
  end % while ls_flag == 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Information output and storage.                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    opt_params.choice = 1; 
%    [iter_hist,flag] = ...
%     output_GPCG(iter,ls_iter,x_c,y_c,g_c,J_c,Active,iter_hist,opt_params);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Subspace Minimization                                                   %
%  Use CG to approx solve A*d_c = -g  with A,g = projected Hess, grad.     %
%  Initializions.                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    x_old = y_c; g_old = g_c; J_old = J_c; Active_old = Active;
    b = -g_c;
    J = J_c;
    f = zeros(size(b));
    resid = (1-Active) .* b;
    J_diff_max = 0;
    cgiter = 0;
    cg_flag = 0;
    cgiter0 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CG iterations.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
    while cg_flag == 0

      cgiter = cgiter + 1;
      dh = resid;
      rd = resid(:)'*dh(:);
    
      %  Compute conjugate direction ph.

      if cgiter == 1,
        ph = dh; 
      else
        betak = rd / rdlast;
        ph = dh + betak * ph;
      end

      %  Form product Ah*ph.
    
      Iph = (1-Active).*ph;
      Ahph = (1-Active).*feval(cost_params.Amult,cost_params,Iph);

      %  Update f and residual.
    
      alphak = rd / (ph(:)'*Ahph(:));
      f = f + alphak*ph;
      resid = resid - alphak*Ahph;
      rdlast = rd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Check Stopping criteria                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      x_i = x_old + f;
      [J_i,g_i] = feval(cost_fn,x_i,cost_params);
      J_diff = J - J_i;
      J = J_i;
      J_diff_max = max(J_diff,J_diff_max);
      if J_diff <= .1*J_diff_max
       % Insufficient decrease.
        [x_c,J_c,g_c,ls_param,cg_flag] = ...
          linesearch_MT(x_old,f,g_old,J_old,cost_fn,cost_params);

        Active = (x_c == 0);
	Bactive = Active.*(g_c>0);        %  Binding set.
        Same = (Active == Active_old);
	AB_diff = (Active == Bactive);

        if min(Same(:)) == 0 | min(AB_diff) == 0
          cg_flag = 1;
        else
         % Active = Active_old. Return to CG.
	  cg_flag = 0;
	  J_diff_max = 0;
        end  
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Information output and storage.                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [iter_hist,flag] = ...
            output_GPCG(iter,ls_iter,cgiter-cgiter0,x_p,x_c,g_c,...
  	       J_c,Active,iter_hist,opt_params);
	if flag ~= 0
	  cg_flag = flag;
	elseif cg_flag == 0
	  x_p = x_c;
          ls_iter = 0;
          cgiter0 = cgiter;
	  iter = iter + 1;
	end
	
      end  % if J_diff ...
    end % CG iteration
    
    total_iters = total_iters + ls_iter + cgiter;
    
  end % while flag == 0

  fprintf('TOTAL CG + LINESEARCH ITERATIONS = %d\n',total_iters); 

  epsilon = .1;
  if flag == 2 | flag == 3 % Gradient or step tolerance reached.
    Active = find(x_c == 0 & abs(g_c)<epsilon);
    if ~isempty(Active)
      g_Active = g_c(Active); 
      figure(4)
        hist(log10(g_Active))
        xlabel('log_{10}(g_{active})')
        ylabel('No. of Occurences')
        title('Hist of active grad components')

      fprintf('# i such that |g(Active)_i| < %0.1e is %d\n',epsilon,length(Active(:)))

      Active = find(x_c == 0 & g_c == 0);
      fprintf('# i such that |g(Active)_i| = 0 is %d\n',length(Active(:)))
    else
      %fprintf('|g(Active)_i|>=%0.1e for all i\n',epsilon)
    end 
  end
  














