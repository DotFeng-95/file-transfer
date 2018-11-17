  function [uh] = vcycle(uh,kappah,bh,level,nu)
%
%  [uh] = vcycle(uh,kappah,bh,level,nu)
%
%  Multigrid V-cycle

  n = max(size(uh));
  Ah =  get_Lmat(kappah);
  L = tril(Ah,-1); 
  D = spdiags(diag(Ah),0,n,n); 
  U = triu(Ah,1);
  unull = ones(n,1)/sqrt(n);  % Basis for null(L).
  
%  Perform nu steps of Gauss-Seidel presmoothing.

  for gsiter = 1:nu
    uh = uh + (D + L)\(bh - Ah*uh);
    uh = uh - (uh'*unull) * unull;  %  Orthogonalize
  end

%  Recursive call to vcycle.

  if level ~= 1,
    level = level-1;
    rh = bh - Ah*uh;
    [r2h] = restrict(rh);
    kappa2h = restrict(kappah)/4;
    u2h = zeros(size(r2h));
    [u2h] = vcycle(u2h,kappa2h,r2h,level,nu);
    uh = uh + prolong(u2h);
  else
    uh = pinv(full(Ah)) * bh;
  end
  
%  Perform nu steps of Gauss-Seidel postsmoothing.  Note that pre- and 
%  post-smoothing operators are adjoints.

  for gsiter = 1:nu
    uh = uh + (D + U)\(bh - Ah*uh);
    uh = uh - (uh'*unull) * unull;  %  Orthogonalize
  end
