N = [6,11,21,41,81];
L = 0;
R = 1;
emax = zeros(5,1);
e0 = zeros(5,1);
e1 = zeros(5,1);

for i = 1:5
    [x,U] = FD1d_bvp(N(i),@f,L,R,@u);
    [e0(i),e1(i),emax(i)] = FD1d_error(x,U,@u);
    X{i} = x;
    UN{i} = U;
end
ue = u(X{5});
plot (X{5} , ue , '-k*', X{1} , UN {1} , '-ro ', X {2} ,...
UN {2} , '-gs ', X{3} , UN {3} ,'-bd ' ,...
X{4} , UN {4} , '-ch',X{5} , UN {5} , '-mx');
title ('The   solution   plot ');
xlabel ('x'); ylabel ('u');
legend ('exact','N=6','N =11','N =21','N =41','N =81');

format short e
disp('     emax           e0           e1');
disp([emax,e0,e1]);