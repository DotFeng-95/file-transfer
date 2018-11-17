  function H_GN = eval_HGN(q0,b,d)
  
%  Evaluate Gauss-Newton approximation to Hessian of least squares 
%  functional.  
  
  nq = max(size(q0));
  h = 1 / nq;
  [n,m] = size(b);
  H_GN = zeros(nq,nq);
  [J0,A0,u0] = eval_Jls(q0,b,d);
  Du = diff(u0)/h;
  for i = 1:nq
    Aiu = zeros(n,m);
    if i == 1
      Aiu(1,:) = u0(1,:)/h;
    elseif i == nq
      Aiu(n,:) = u0(n,:)/h;
    else
      Aiu(i-1,:) = -Du(i-1,:);
      Aiu(i,:)   = Du(i-1,:);
    end
    si = A0 \ (A0 \ Aiu);

    H_GN(1,i) = sum(u0(1,:).*si(1,:)) / h;
    for j = 2:n
      H_GN(j,i) = sum(-Du(j-1,:).*si(j-1,:) + Du(j-1,:).*si(j,:));
    end
    H_GN(n+1,i) = sum(u0(n,:).*si(n,:)) / h;
  end
  