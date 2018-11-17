  function [J,g] = eval_pls(f,S,d,L,alpha)
  
%  [J,g] = eval_pls(f,S,d,L,alpha)
%
%  Evaluate penalized least squares cost functional
%      J = .5*(||Kf-d||^2 + alpha*f'*L*f)
%  and its gradient
%      g = K'*(Kf-d) + alpha*L*f.

  resid = real(ifft2(S .* fft2(f))) - d;
  Lf = L*f(:);
  J = .5 * (norm(resid(:))^2 + alpha*f(:)'*Lf);
  
  if nargout==1, return, end
  
  [nx,ny] = size(f);
  Kstar_r = real(ifft2(conj(S).*fft2(resid)));
  g = Kstar_r + alpha*reshape(Lf,nx,ny);
  
