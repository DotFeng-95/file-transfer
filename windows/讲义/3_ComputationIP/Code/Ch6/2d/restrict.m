  function [v] = restrict(u)
%
%  [v] = restrict(u)
%
%  2-D restriction operator based on piecewise constant basis functions.

  [nx,ny] = size(u);
  if nx == ny
    umat = u;
    n = nx;
    sflag = 1;
  else
    nsq = max(nx,ny);
    n = sqrt(nsq);
    umat = reshape(u,n,n);
    sflag = 0;
  end
  nd2 = n/2;
  vmat = zeros(nd2,nd2);
  for i = 1:nd2
    vi = sum([umat(2*i-1,:); umat(2*i,:)]);
    vmat(i,:) = (vi(1:2:n-1) + vi(2:2:n))/4;
  end

  if sflag == 1
    v = vmat;
  else
    v = vmat(:);
  end