%  test0.m
%
%  Test cg.m.

  n = input(' System size n = ');
  randn('state',0);  %  Initialize random number generator.
  A = randn(n,n);
  A = A'*A + eye(n);
  b = randn(n,1);
  x_true = A\b;
  x0 = zeros(n,1);
  
  params.matrix = A;
  params.max_cg_iter = input(' Max CG iterations = ');
  params.cg_step_tol = 0;
  params.grad_tol_factor = 1e-6;
  params.cg_io_flag = 1;
  params.cg_figure_no = 1;
  
  [x,term_code] = cg(x0,b,params,'a_mult0','m_inv0');