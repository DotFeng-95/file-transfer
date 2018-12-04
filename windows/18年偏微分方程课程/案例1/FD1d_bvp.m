function [x,U] = FD1d_bvp(N,f,a,b,u)
h = (b-a)/(N-1);
x = (a:h:b)';
c1 = -1/h/h;
c2 = 2/h/h;
g = [c1*ones(1,N-2),0];
c = [0, c1*ones(1,N-2)];
d = [1,c2*ones(1,N-2),1];
A = diag(g,-1)+diag(d)+diag(c,1);
rhs = f(x);
rhs(1) = u(x(1));
rhs(N) = u(x(N));
U = A \ rhs;
end