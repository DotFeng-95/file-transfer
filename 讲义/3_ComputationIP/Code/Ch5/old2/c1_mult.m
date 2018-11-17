  function y_array = c1_mult(D,x_array);

%  Compute array(C*vec(x_array)), where C is a block Circulant matrix.
%  D is a 3-D array containing the blocks D_i for which 
%     C = tensor(I,F)' * P' * diag(D_i) * P * tensor(I,F).
%  See comments in level1_bttb for details.
  
  %  Extract and check array sizes.
  
  if ndims(D) ~= 3
    fprintf('***** First argument to c1_mult must be a 3-D array.\n');
    return
  end
  if ndims(x_array) ~= 2
    fprintf('***** Second argument to c1_mult must be a 2-D array.\n');
    return
  end
  [nx,ny] = size(x_array);
  [n1,n2,n3] = size(D);
  if nx ~= n1
    fprintf('***** 1st dimension of D must equal nx.\n');
    return
  end
   if (ny ~= n2 | n2 ~= n3)
    fprintf('***** 2nd and 3rd dimension of D must equal ny.\n');
    return
  end
   
%  Compute array(C*vec(x)).
  
  Z = zeros(ny,nx);
  xhat = fft(x_array);  %  Apply 1-D FFT's to columns.
  xhat = xhat.';        %  Permute rows and columns.
  for i = 1:nx
    vi = xhat(:,i);     %  ith column of xhat.
    Di = zeros(ny,ny);
    Di(:,:) = D(i,:,:); %  ith block; dimension is ny X ny.
    Z(:,i) = Di*vi;
  end
  
  y_array = real(ifft(Z.'));
  
  