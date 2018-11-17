%  Demo_Tmult.m
%
%  Show that BTTB matrix-vector multiplication can be carried out
%  using zero extension combined with 2-D fft's.

%  Get data.

nnx = nfx/2;
nny = nnx;
t = extract(fftshift(PSF),2*nnx-1,2*nny-1);
f_array = extract(f_true,nnx,nny);

%  Compute T*f using matrix-vector multiplication.

fvec = f_array(:);
T = bttb(t);
Tfvec = T*fvec;
Tf_array = reshape(Tfvec,nnx,nny);
figure(1), imagesc(Tf_array), colorbar

%  Use fft's to compute T*f.

nx2 = 2*nnx;
ny2 = 2*nny;
c_ext = zeros(nx2,ny2);
c_ext(2:nx2,2:ny2) = t;
c_ext = fftshift(c_ext);
f_ext = zeros(nx2,ny2);
f_ext(1:nnx,1:nny) = f_array;
g_ext = real(ifft2(fft2(c_ext).*fft2(f_ext)));
g_array = g_ext(1:nnx,1:nny);
figure(2), imagesc(g_array), colorbar
