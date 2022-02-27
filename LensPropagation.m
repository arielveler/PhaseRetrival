function Eout = LensPropagation(Ein,f,lambda,x,y)
% This function will propagate a wave througth a thin lens using Fresnel propagation within the paraxial regime.
% Input:
% Ein - initial 2D electrical field distribution at z=0 ([N1 x N2] complex)
% f - focal length of the lens
% lam0 - the central wavelength
% x, y - the x, y coordinates of the field ([N1 x N2] double)
% Output:
% Eout - the field after the thin lens.

% add restrictions on lens focal length - Nyquist of the lens (2*pi foldings) 
% dx = max(abs(x(1,2)-x(1,1)),abs(x(2,1)-x(1,1)));
% fn = max(size(x))*dx^2/lambda;
% if f < fn
%     error(['focal length of the lens is too small. minimum for this system is ',num2str(fn), '[m]']);
% end

k0 = 2*pi/lambda;
propagator = exp(-1i*k0/(2*f)*(x.^2+y.^2));
%try:
a = zeros(size(propagator));
a(size(propagator,1)/2-size(propagator,1)/4:size(propagator,1)/2+size(propagator,1)/4,size(propagator,2)/2-size(propagator,2)/4:size(propagator,2)/2+size(propagator,2)/4)=1;
Eout = propagator.*Ein;
end