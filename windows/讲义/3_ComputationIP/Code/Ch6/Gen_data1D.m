%  Gen_data1D.m
%
%  First set up FEM discretization of 1-d steady-state diffusion equations,
%    -d/dx(q(x) du_l/dx) = f_l(x),  0 < x < 1,  l = 1,2,
%    u_l(0) = u_l(1) = 0.
%  The forcing functions f_l(x) correpond to point sources at x=1/3
%  and x=2/3,
%       f_1(x) = delta(x-1/3),
%       f_2(x) = delta(x-2/3).
%  Then solve for u_l and generate data
%    d_l = u_l + noise.
%  The index l corresponds to the lth experiment.

%  Set up FEM mesh.

  n_nodes = 50; %%% input(' No. of interior FEM nodes = ');
  percent_error = 1; %%% input(' Noise level percentage = ');
  h = 1 / (n_nodes +1);
  x = [0:h:1]';
  x0 = [h:h:1-h]';
  x_mid = [h/2:h:1-h/2]';
  q_true = 1 - .25*exp(-(x_mid - .45).^2/.01);

%  Set up left and right BC.

  u_left = 0;
  u_right = 0;

%  Compute load vectors b(:,l), l=1,2. For interior nodes i,
%  b(i,l) = phi_i(l/3), where phi_i(x) is the ith piecewise linear
%  "hat" basis function.

  b = zeros(n_nodes,2);
  b(floor(n_nodes/3),1) = 1;
  b(floor(2*n_nodes/3),2) = 1;
  
  %  Add on contributions from boundary conditions.
  
  b(1,:) = b(1,:) + u_left * q_true(1) / h * ones(1,2);
  b(n_nodes,:) = b(n_nodes,:) + u_right * q_true(n_nodes+1) / h * ones(1,2);

%  Compute stiffness matrix A(q).

  Adiag = (q_true(1:n_nodes) + q_true(2:n_nodes+1)) / h;
  Asub = -q_true(2:n_nodes) / h;
  Asuper = Asub;
  A = spdiags([[Asub;0] Adiag [0;Asuper]], [-1 0 1], n_nodes,n_nodes);

%  Solve systems A*u_l = b_l and append BC.

  u0 = A\b;
  u = [u_left*ones(1,2); u0; u_right*ones(1,2)];

%  Generate noisy data.

  randn('state',0);  %  Reset random number generator to initial state.
  noise = percent_error / 100 * norm(u) * randn(n_nodes,2) / sqrt(n_nodes);
  d = u0 + noise;

%  Plot solution u.

  figure(1)
    plot(x,u,'-', x0,d,'o')
    xlabel('x axis')
    title('Exact (-) and Noisy (o) Data')
