pde = model_data();
[X,T,U] = heat_equation_fd1d (100 ,10000 , pde ,'forward');
showvarysolution (X,T,U);
showsolution (X,T,U); 

[X,T,U] = heat_equation_fd1d (100 ,100 , pde ,'backward');
showvarysolution (X,T,U);
showsolution (X,T,U);

[X,T,U] = heat_equation_fd1d (100 ,100 , pde ,'crank - nicholson');
showvarysolution (X,T,U);
showsolution (X,T,U);