%  Spectrum.m
%
%  Compute spectra of preconditioned systems corresponding to 
%  level 1 and level 2 preconditioners.

%  Generate t array and T = bttb(t).

nnx = nfx/2;
nny = nnx;
t = extract(fftshift(PSF),2*nnx-1,2*nny-1);
T = bttb(t);

%  Generate matrix A = T'*T + alpha*I.

alph = input(' alpha = ');
N = nnx*nny;
A = T'*T + alph * eye(N);
eA = sort(eig(A));

%  Compute level 1 and level 2 preconditioning matrices.

[C1] = level1(t);
[c2,C2] = level2(t);
M1 = C1'*C1 + alph * eye(N);
M2 = C2'*C2 + alph * eye(N);


%  Compute eigenvalues for preconditioned systems.

eP1 = sort(eig(M1\A));
eP2 = sort(eig(M2\A));

indx = [1:N]';
semilogy(indx,eA,'o-', indx,eP1,'*-', indx,eP2,'x-')

