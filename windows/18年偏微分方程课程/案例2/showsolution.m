function showsolution(X,T,U)
[x,t] = meshgrid(X,T);
mesh(x,t,U');
xlabel('X');
ylabel('T');
zlabel('U(X,T)');
end