  function [J,Kfps,g] = eval_lhd(f,S,dps,alpha,L,sigsq)
  
%  [J,g] = eval_lhd(f,S,dps,alpha,sigsq)
%
%  Evaluate penalized likelihood cost functional
%      J = sum(Kfps - dps.*log(Kfps)) + alpha/2*norm(f)^2
%  and its gradient
%      g = K'*((Kfps-dps)./Kfps) + alpha * f,
%  where Kfps = K*f + sigsq.

  Kfps = mat_prod(S,f) + sigsq;
  Lf = L*f(:);
  J =  sum(sum(Kfps - dps.*log(Kfps))) + alpha/2*f(:)'*Lf;
  
  if nargout==1, return, end
  
  [Nx,Ny] = size(f);
  g = mat_prod(conj(S),(Kfps-dps)./Kfps) + alpha*reshape(Lf,Nx,Ny);
  
