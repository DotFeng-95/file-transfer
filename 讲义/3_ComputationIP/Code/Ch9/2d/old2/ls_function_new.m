  function [J,g] = ls_function(f,J_params)
  
%  Evaluate penalized least squares cost functional
%      J = .5*(||Tf-d||^2 + alpha*||f||^2),  T = [T_1,..,T_k]'
%  where the T_i's are matrices, and the gradient 
%  g = T*(Tf-d)+alpha*f.
 
  global TOTAL_COST_EVALS TOTAL_FFTS
  TOTAL_COST_EVALS = TOTAL_COST_EVALS + 1;

  gam = J_params.gam;
  S = J_params.S;
  d = J_params.d;
  n_frames = J_params.n_frames;
  
  J = 0;
  g_temp = zeros(size(f));
  Fourier_f = fft2(f);
  TOTAL_FFTS = TOTAL_FFTS + 1;
  for i = 1:n_frames
    temp_i = (1/n_frames)*S{i} .* Fourier_f;
    resid_i = real(ifft2(temp_i)) - d{1,i};
    TOTAL_FFTS = TOTAL_FFTS + 1;
    J = J + .5*norm(resid_i(:))^2;
    if nargout == 2                       % Compute grad(J).
      g_temp = g_temp + 1/n_frames*conj(S{i}).*(temp_i - fft2(d{1,i}));
      TOTAL_FFTS = TOTAL_FFTS + 1;
    end
  end

  J = J + .5*gam*norm(f(:))^2;
  
  if nargout == 2
    g = real(ifft2(g_temp))+gam*f;
    TOTAL_FFTS = TOTAL_FFTS + 1;
  end
