function Eout = FreePropagation(Ein,z,lambda,kx,ky)
% This function will propagate a wave using Fresnel propagation within the paraxial regime.
% Input:
% Ein - initial 2D electrical field distribution at z=0 - E(x,y,z=0)([N1 x N2] complex)
% z - distance of propagation
% lam0 - the central wavelength
% kx, ky - the x, y components of wave number (at K space) ([N1 x N2] double)
% Output:
% Eout - the field after the propagation - E(x,y,z).
k0 = 2*pi/lambda;
propagator = exp(1i*k0.*z.*sqrt(1-(kx.^2+ky.^2)/k0.^2));
Eout = ifft2(ifftshift(propagator.*fftshift(fft2(Ein))));
end