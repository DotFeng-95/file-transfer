%  Compute and plot spectrum.

  disp(' ... computing spectra of operators ...');

  dudq = deriv(q_mat,b);
  H_ls = dudq'*dudq;
  alphaL = alpha * get_Lmat(kappa);
  H = H_ls + alphaL;
  spec_H = sort(real(eig(H)));
  Hperp = H - w0(:)*w0(:)'/c;
  CinvHperp = pinv(full(alphaL)) * Hperp;
  spec_PCG = sort(real(eig(CinvHperp)));
  
  figure(4)
  subplot(221)
    semilogy(spec_H,'o')
    xlabel('index i')
    ylabel('lambda_i')
    title('Spectrum of H')
  subplot(223)
    semilogy(spec_PCG(2:N),'o')
    xlabel('index i')
    ylabel('lambda_i')
    title('Spectrum of H_PCG')

