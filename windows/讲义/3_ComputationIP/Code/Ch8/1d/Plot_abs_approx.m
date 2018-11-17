%  Plot_abs_approx.m

%  Plot approximations to the absolute value function.

  beta = input(' Parameter beta = ');
  epsilon = input(' Parameter epsilon = ');
  t = [-1.5:.002:1.5]';
  t2 = t.^2;
  
  phi = abs(t);
  phi_beta = sqrt(t2 + beta^2);
  eps2 = epsilon^2;
  i1 = (t2 < eps2);
  i2 = 1 - i1;
  phi_eps = i1.*t2/epsilon/2 + i2.*(sqrt(t2)-epsilon/2);
  
  plot(t,phi,'-.', t,phi_beta,'--', t,phi_eps,'-')
  xlabel('t')
  