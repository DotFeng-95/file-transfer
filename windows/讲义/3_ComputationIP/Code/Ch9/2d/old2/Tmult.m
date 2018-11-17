  function[Tx] = Tmult(T,x);
  
%  Compute ctranspose(T)*b where T = [T1 ... Tk]' and 
%  Ti*x = real(ifft2(S{i}.*fft2(x))). 

  S = T.S;
  n_frames = T.n_frames;
  
  Tx = zeros(size(x{1,1}));
  
  for i = 1:n_frames
    Tx = Tx + conj(S{i}).*fft2(x{1,i});
  end
  
  Tx = (1/n_frames)*real(ifft2(Tx));
