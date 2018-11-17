  function[x_new,phi_alpha,phi_p_alpha,alpha,opt_flag] = ... 
                          linesearch_MT(x,p,g,q,q_fn,cost_params)
 
%  Initialize parameters and vectors.
  
  phi_0 = q;
  phi_p_0 = g(:)'*p(:);  
  Ap = feval(cost_params.Amult,cost_params,p);
  phi_pp_0 = p(:)'*Ap(:);
  alpha = -phi_p_0/phi_pp_0;
 
  I = find(p < 0);    

  if isempty(I) 	   % Step towards the interior.
    x_new = x + alpha*p;
    [phi_alpha,phi_p_alpha] = feval(q_fn,x_new,cost_params);  
    opt_flag = 1;
    return
  end
  
  t = -x(I)./p(I);
  J = find(t > 0);
  if isempty(J)
    beta_1 = 0;
  else
    t = t(J);
    beta_1 = min(t(:));
  end  

  if alpha < beta_1
    x_new = x + alpha*p;
    [phi_alpha,phi_p_alpha] = feval(q_fn,x_new,cost_params);  
    opt_flag = 2;
    return
  end

  x_new = max(x + alpha*p,0);
  [phi_alpha,phi_p_alpha] = feval(q_fn,x_new,cost_params);  
  opt_flag = 0;
  i = 0;
  max_iter = 20;
  mu = .1;

%  Iteration.
  
  while opt_flag == 0

    i = i + 1;
       
   % Check sufficient decrease condition.
    
    x_new_minus_x = x_new - x;
    if phi_alpha < phi_0 + mu*g(:)'*x_new_minus_x(:) %mu*norm(x_new_minus_x(:))^2 
      opt_flag = 3;
      return
   end

   % Minimize the quadratic which interpolates phi_0, phi_p_0 and phi_alpha. 
    
    alpha_new = -.5*phi_p_0*(alpha^2/(phi_alpha - alpha*phi_p_0 - phi_0));

   % Determine new alpha.

    m = median([.01*alpha,alpha_new,.5*alpha]);
    alpha = max(m,beta_1);

   % Evaluate phi(alpha).

    x_new = max(x + alpha*p,0);
    [phi_alpha,phi_p_alpha] = feval(q_fn,x_new,cost_params);  

    if i >= max_iter
      disp('*************** Linesearch Error ***************');
      opt_flag = 4;
    end

  end % while flag == 0    


