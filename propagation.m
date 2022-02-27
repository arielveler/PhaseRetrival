function PsiZ = propagation(Psi0,steps,z,k0,x,kx,n0,dn)
% This function will propagate a wave acccording to what we have learn in
% class, when at every iteration we propagate the wave at the z axis using
% Fourier transform and Fresnel propagator within the paraxial regime, and
% the propagation of a "potential" of the refracting index.
% Input:
% Psi0 - initial wave at z=0
% steps - number of steps in z propagation
% z - distance of propagation
% ko - the wave number
% x - x axis
% kx - the x component of wave number
% n0 - main refracting index
% dn - small change in refracting index
% Output:
% PsiZ - the image of the propagation along Z axis.
dz = z/steps;
PsiZ = zeros([steps,numel(x)]);
PsiZ(1,:) = Psi0;
PsiZ(2,:) = exp(-1i.*k0./n0.*dn.*dz).*ifft(fftshift(exp(1i.*kx.^2.*(dz/2)/(2.*k0))).*fft(Psi0));
for ii=3:steps-1
    PsiZ(ii,:) = exp(-1i.*k0./n0.*dn.*dz).*ifft(fftshift(exp(1i.*kx.^2.*(dz)/(2.*k0))).*fft(PsiZ(ii-1,:)));
end
PsiZ(steps,:) = ifft(fftshift(exp(1i.*kx.^2.*(dz/2)/(2.*k0))).*fft(PsiZ(ii,:)));