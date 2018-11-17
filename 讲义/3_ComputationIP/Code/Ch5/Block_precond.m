%  Block_precond.m
%
%  Apply PCG to solve regularized BTTB system.

  alpha = input(' Regularization parameter alpha = ');
  precond_flg=input(' Enter 0 for no precond; 1 for level 1; 2 for level 2; 3 for circ ext: ');

%  Get data.

  nnx = nfx;
  nny = nnx;
  n = nnx*nny;
  t = extract(fftshift(PSF),2*nnx-1,2*nny-1);
  f_array = extract(f_true,nnx,nny);
  d_array = dat(1:nnx,1:nny);
  
%  Compute the Fourier representer of the block cirulant extension of
%  T = BTTB(t). This is used to compute matrix-vector products T*x
%  using FFT's.

  nx2 = 2*nnx;
  ny2 = 2*nny;
  t_ext = zeros(nx2,ny2);
  t_ext(2:nx2,2:ny2) = t;
  t_ext_hat = fft2(fftshift(t_ext));
  
%  Store information needed to compute A*x, where A = T'*T + alpha*L.

  params.array_size = [nnx,nny];
  params.t_extension_hat = t_ext_hat;
  params.penalty_matrix = speye(n);
  params.reg_param = alpha;
  mat_mult_fn = 'a_mult';

%  Compute b_array = array(T'*d), where T = BTTB(t) and d is 
%  the data vector.
  
  d_ext = zeros(nx2,ny2);
  d_ext(1:nnx,1:nny) = d_array;
  b_array = real(ifft2( conj(t_ext_hat) .* fft2(d_ext) ));
  b_array = b_array(1:nnx,1:nny);

%  Construct preconditioner.
    
  if precond_flg == 1
    
    %  Construct level 1 preconditioner, M1 = C1'*C1 + alpha*I,
    %  where C1 is the level 1 BCCB approximation to T.
  
    D = level1_bttb(t);
    M1 = zeros(nnx,nny,nny);
    for i = 1:nnx
      Di = zeros(nny,nny);
      Di(:,:) = D(i,:,:);
      M1(i,:,:) = Di'*Di + alpha*eye(nny);
    end
      
    params.precond_representer = M1;
    mat_inv_fn = 'm1_inverse';
  
  elseif precond_flg == 2
    
    %  Construct level 2 preconditioner, M2 = C2'*C2 + alpha*I,
    %  where C2 is the level 2 BCCB approximation to T.
  
    [c2] = level2(t);
    m2_hat = abs(fft2(c2)).^2 + alpha;
    params.precond_representer = m2_hat;
    mat_inv_fn = 'm2_inverse';
    
  elseif precond_flg == 3
    
    %  Construct block circulant preconditioner.
  
    params.precond_representer = 1 ./ (abs(t_ext_hat).^2 + alpha);
    mat_inv_fn = 'circ_ext_inverse';
    
  end
  
  %  CG iteration and output information.
    
  params.max_cg_iter = input(' Max CG iterations = ');
  params.cg_step_tol = 0;
  params.grad_tol_factor = input(' CG gradient stopping tolerance = ');
  params.cg_io_flag = 1;
  params.cg_figure_no = 1;

  f0 = zeros(nnx,nny);   %  Initial guess.
  if precond_flg == 0
    [f,iter_hist,term_code] = cg(f0,b_array,params,mat_mult_fn);
  else
    [f,iter_hist,term_code] = cg(f0,b_array,params,mat_mult_fn,mat_inv_fn);
  end

  figure(2)
    imagesc(reshape(f,nnx,nny)), colorbar
    title('Tikhonov/CG Reconstruction')
