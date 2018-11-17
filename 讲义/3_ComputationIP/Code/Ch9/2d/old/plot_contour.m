%  Plot_contour.m
%

 figure(1)
  h = .05;
  xx = [-0.5:h:1]';
  [x,y] = meshgrid(xx,xx);
  z = sqrt((x-0.5).^2 + (y+0.5).^2/9);
  h1 = .1;
  x1 = [-0.5:h1:1]';
  [xx,yy] = meshgrid(x1,x1);
  z1 = sqrt((xx-0.5).^2 + (yy+0.5).^2/9);
  [px,py] = gradient(z1,h1,h1);
  contour(x,y,z)
  hold on
  quiver(xx,yy,px,py);
  plot([0,0],[1,0])
  plot([0,1],[0,0])
  xlabel('x axis')
  ylabel('y axis')
  hold off

 figure(2)
  h = .05;
  xx = [-0.5:h:1]';
  [x,y] = meshgrid(xx,xx);
  z = sqrt((x-0.5).^2 + (y+0.5).^2/9);
  contour(x,y,z)
  hold on
  plot([0,0],[1,0])
  plot([0,1],[0,0])
  xlabel('x axis')
  ylabel('y axis')
  hold off
