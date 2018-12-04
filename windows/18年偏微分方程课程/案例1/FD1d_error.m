function [e0,e1,emax] = FD1d_error(x,U,u_exact)
N = length(x);
h = (x(end)-x(1))/(N-1);
ue = u_exact(x);
ee = ue-U;
e0 = h*sum(ee.^2);
e1 = sum((ee(2:end)-ee(1:end-1)).^2)/h;
e0 = sqrt(e0);
e1 = sqrt(e1);
emax = max(abs(ue-U));
end
