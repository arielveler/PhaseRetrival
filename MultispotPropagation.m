function Eout = MultispotPropagation(Ein,theta,n,lambda,x,y)
% This function will propagate a wave througth a Multi-spot component yelding
% n spots in theta degrees in propagation between them using Fresnel propagation within the paraxial regime.
% Input:
% Ein - initial 2D electrical field distribution at z=0 ([N1 x N2] complex)
% theta - anlge of deviation between adjacent spots (double).
% n - number of spots in each axis ([1 x 2] integer)
% lambda - the central wavelength (double)
% x, y - the x, y coordinates of the field ([N1 x N2] double)
% Output:
% Eout - the field after the thin lens ([N1 x N2] complex).
k0 = 2*pi/lambda;
dk = sind(theta).*k0;
k_mat = zeros(size(Ein));

% option b:
kx = -dk(1)*(n(2)-1)/2:dk(1):dk(1)*(n(2)-1)/2;
ky = -dk(2)*(n(1)-1)/2:dk(2):dk(2)*(n(1)-1)/2;
for i=1:n(2)
    for j=1:n(1)
        k_mat = k_mat + exp(1i.*(kx(i).*x+ky(j).*y));
    end
end
mask = k_mat;

Eout = mask.*Ein;
end