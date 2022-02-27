function [Eout,centers,sub_size] = MicroLensPropagation(Ein,n,b,pix_size,f,lambda,x,y)
% This function will propagate a wave througth a thin micro-lens array using Fresnel propagation within the paraxial regime.
% Input:
% Ein - initial 2D electrical field distribution at z=0 ([N1 x N2] complex)
% f - focal length of the lenses in the micro-lens array.
% lam0 - the central wavelength
% x, y - the x, y coordinates of the field ([N1 x N2] double)
% n - number of micro-lenses (spots) in the micro-lens array. 
% Output:
% Eout - the field after the thin micro-lens array.
k0 = 2*pi/lambda;
Eout = zeros(size(Ein));
propagator = zeros(size(Ein));
s = flip(b.*n); % physical size of microlens array [col,row]
imsize = pix_size.*flip(size(Ein)); % physical size of whole image [col,row]
p0 = (imsize-s)./2; % phisycal size of starting point (top left) of the microlens array in image [col,row]
p0_ind = flip(p0./pix_size); % starting point (top left) of the microlens array in image in pixels [row,col]
b_ind = flip(b./pix_size);% distance between centers (or sub image size) in pixels [col,row]
if any(s>imsize)
    error('MicroLens size must be smaller than the image size');
end
for i=1:n(1)
    for j=1:n(2)
        ind_row = round((p0_ind(1)+(i-1)*b_ind(1):p0_ind(1)+i*b_ind(1))); ind_col = round((p0_ind(2)+(j-1)*b_ind(2):p0_ind(2)+j*b_ind(2)));
        center(i,j,:) = [y(ceil(mean(ind_row)),1),x(1,ceil(mean(ind_col)))];
        lens = exp(-1i*k0/(2*f)*((x-center(i,j,2)).^2+(y-center(i,j,1)).^2));
        propagator(ind_row,ind_col) = lens(ind_row,ind_col);
    end 
end
Eout = propagator.*Ein;
centers = center;
sub_size = b_ind;
end