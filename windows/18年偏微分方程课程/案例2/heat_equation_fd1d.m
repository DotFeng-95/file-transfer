function [X,T,U] = heat_equation_fd1d(NS,NT,pde,method)
[X,h] = pde.space_grid(NS);
[T,tau] = pde.time_grid(NT);
N = length(X);
M = length(T);
r = pde.a()*tau/h/h;
if r >= 0.5 && ismember(method,{'F','f','forward'})
    error('时间空间离散不满足向前差分的稳定条件！')
end
U = zeros(N,M);
U(:,1) = pde.u_initial(X);
U(1,:) = pde.u_left(T);
U(end,:) = pde.u_right(T);
switch(method)
    case{'F','f','forward'}
        forward();
    case{'B','b','backward'}
        backward();
    case{'CN','cn','crank-nicholson','Crank-Nicholson'}
        crank_nicholson();
    otherwise
        disp(['Sorry, I do not know your ',method]);
end
%%向前差分
    function forward()
        d = 1-2*ones(N-2,1)*r;
        c = ones(N-3,1)*r;
        A = diag(c,-1)+diag(c,1)+diag(d);
        for i = 2:M
            RHS = tau*pde.f(X,T(i));
            RHS(2) = RHS(2) + r*U(1,i-1);
            RHS(end-1) = RHS(end-1) + r*U(end,i-1);
            U(2:end-1,i) = A*U(2:end-1,i-1) + RHS(2:end-1);
        end
    end
%%向后差分
   function backward()
       d = 1 + 2*ones(N-2,1)*r;
       c = -ones(N-3,1)*r;
       A = diag(c,1) + diag(c,-1) + diag(d);
       for i = 2:M
           RHS = tau*pde.f(X,T(i));
           RHS(2) = RHS(2) + r*U(1,i);
           RHS(end-1) = RHS(end-1) + r*U(end,i);
           U(2:end-1,i) = A\(U(2:end-1,i-1) + RHS(2:end-1));
       end
   end
%%六点对称
    function crank_nicholson()
        d1 = 1 + ones(N-2,1)*r;
        d2 = 1 - ones(N-2,1)*r;
        c = 0.5*ones(N-3,1)*r;
        A1 = diag(-c,-1) + diag(-c,1) + diag(d1);
        A0 = diag(c,-1) + diag(c,1) + diag(d2);
        for i = 2:M;
            RHS = tau*pde.f(X,T(i));
            RHS(2) = RHS(2) + 0.5*r*(U(1,i)+U(1,i-1));
            RHS(end-1) = RHS(end-1) + 0.5*r*(U(end,i)+U(end,i-1));
            U(2:end-1,i) = A\(A0*U(2:end-1,i-1)+RHS(2:end-1));
        end
    end    
end