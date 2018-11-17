%  test.m
%

  nq = max(size(q_true));
  q0 = .5*ones(nq,1);

  fd_grad;
  
%  costate;
  [g_costate,J,A,u] = eval_gradls(q0,b,d);