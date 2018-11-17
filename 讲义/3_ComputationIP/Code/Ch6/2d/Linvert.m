  function [u,residnormvec] = Linvert(u,kappa,b,pcg_iter,pcg_tol,level,nu)
%
%  Compute the (pseudoinverse) solution to (the CCFD discretization of)
%           L[u] = -div( kappa(x,y) grad u) = f(x), 0<x,y<1,
%           du/dn = 0,                              on the boundary.
%  Note that null(L) consists of constant vectors.
%  The method of solution is preconditioned conjugate gradient with
%  a multigrid preconditioner.

  n = max(size(u));
  L = get_Lmat(kappa);
  npcg = 0;
  residrat = 1;
  resid = b - L*u;
  residnormvec = norm(resid);
  
  while npcg < pcg_iter & residrat > pcg_tol
    npcg = npcg + 1;

    %  Apply multigrid preconditioner.

    d = zeros(n,1);
    d = vcycle(d,kappa,resid,level,nu);

    rd = resid'*d;
    if npcg == 1,
      p = d; 
    else
      betak = rd / rdlast;
      p = d + betak * p;
    end
    Lp = L*p;
    alphak = rd / (p'*Lp);
    u = u + alphak*p;
    resid = resid - alphak*Lp;
    residnormvec = [residnormvec; norm(resid)];
    rdlast = rd;
    residrat = residnormvec(npcg)/residnormvec(1);
  end
