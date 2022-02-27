%% BPM symulation - Ariel Veler

clear all; clc; close all;
addpath(genpath(pwd)); % add all subfolders of current folder to path

% initials:
Nrow = 8192; Ncol = 8192;
lambda = 1030e-9; %[m] central wavelenght
n_row = 9; n_col = 9; n = [n_row, n_col]; % number of spots at the pinholes plane
theta = [0.97 0.97]; %[deg] seperation angle between spots (for the MS)
pixel_pitch = [1 1].*3.45e-6; % [m]
f1 = 100e-3; f2 = 100e-3; %[m] focal length of the lenses 
D = 120e-6; % cut in FWHM of the probe
b = 2e-3; % distance between pinholes
d = 20e-3;
fms = pi*b*D/4/lambda; %lens for collimating the spots comming from the MS
resolution = 10e-6; %resolution desired on object

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig = D./(2*(2*log(2))^0.5); % FWHM to sigma
sig_pix = sig./min(pixel_pitch); %for gauss probe
Im_size = [Nrow Ncol];
Center = Im_size./2;
M = f2./f1; % enlarging
FOV = n*b*d/f1; % FOV estimation on the object plane
% overlap = 1-(pi*D*b*d)/(lambda*sqrt(1+d*(lambda*D^2)/(pi*f1^2))); %overlap estimation between probes on the object plane

% definition of grid:   

k0 = 2*pi/lambda;
sensor = Im_size.*pixel_pitch; % sensor physical size

ax = sensor(2)/M; ay = sensor(1)/M; % pinhol plane physical size
dx = ax/Ncol; dy = ay/Nrow;
[x,y] = meshgrid(-ax/2:dx:ax/2-dx, -ay/2:dy:ay/2-dy);

Wk0x = 2*pi/dx; Wk0y = 2*pi/dy; % [Hz] largest element of frequencies for fourier transform
dkx = Wk0x/Ncol; dky = Wk0y/Nrow; % [Hz] smallest element of frequency
[kx,ky] = meshgrid(-Wk0x/2:dkx:Wk0x/2-dkx,-Wk0y/2:dky:Wk0y/2-dky);

%% settings
% initial wave at z=0

p = GaussProbe(sig, pixel_pitch, lambda, Im_size, Center);
[cen_col,cen_row] = meshgrid(Ncol/n_col/2:Ncol/n_col:Ncol-Ncol/n_col/2, Nrow/n_row/2:Nrow/n_row:Nrow-Nrow/n_row/2);
% centers = cat(3,cen_row,cen_col); %+rand(size(cen))*rf/(2*sqrt(2))

[centers_x,centers_y] = meshgrid(round(linspace(Center(2)-b*(n(2)-1)/2/pixel_pitch(1),Center(2)+b*(n(2)-1)/2/pixel_pitch(1),n(2))),round(linspace(Center(1)-b*(n(1)-1)/2/pixel_pitch(2),Center(1)+b*(n(1)-1)/2/pixel_pitch(2),n(1))));
centers = cat(3,centers_y,centers_x);
[im0,cen,sub_size] = p.pinhols(n,centers);

%% object

% load object:
Create_numbers_image();
amp = imresize(a7,256/size(a7,1));
phase = imresize(a8,256/size(a8,1));

[amp,phase] = cut_min_ind(amp,phase);
% amp = normalize_image(amp); phase = normalize_image(phase);
amp = normalize_image(abs(1-amp)); phase = normalize_image(abs(1-phase));

obj = amp.*exp(1i.*phase);
obj = imresize(obj,10);
obj = padarray(obj,(Im_size-size(obj))/2,1,'both');

%%

z1 = f1;
im1 = FreePropagation(im0,z1,lambda,kx,ky);
im2 = LensPropagation(im1,f1,lambda,x,y);
z2a = f1-d;
im3 = FreePropagation(im2,z2a,lambda,kx,ky);
im4 = im3.*obj;
z2b = f2+d;
im5 = FreePropagation(im4,z2b,lambda,kx,ky);
im6 = LensPropagation(im5,f2,lambda,x,y);
z3 = f2;
im7 = FreePropagation(im6,z3,lambda,kx,ky);


im8 = abs(im7).^2; % taking only intensity
im9 = rot90(im8,2); % flip to arrang the centers according to probes


% cutting probes from the image. result in intensity units(field magnitude square)
for i=1:n(1)
    for j=1:n(2)
        probes(:,:,i,j) = im9(cen(i,j,1)-round((sub_size(1)-1)/2):cen(i,j,1)+round((sub_size(1)-1)/2),cen(i,j,2)-round((sub_size(2)-1)/2):cen(i,j,2)+round((sub_size(2)-1)/2));
    end
end
s = size(probes,[1 2]);
% [probes,s] = cut_probes(im9,cen);
cen_obj = (cen-repmat(reshape(Im_size./2,[1 1 2]),n)).*(d./f1); % centes in object plane in size of the probe image
cen_obj(:,:,1) = cen_obj(:,:,1)*(M*b*dy)/(lambda*f2)+repmat(s(1)/2,n); % centes in object plane, when image is in size of cutted image on camera plane
cen_obj(:,:,2) = cen_obj(:,:,2)*(M*b*dx)/(lambda*f2)+repmat(s(2)/2,n); % centes in object plane, when image is in size of cutted image on camera plane

%% calibrate initial guess of probe
[pr0,~,~] = p.pinhols([1 1],reshape(Center,[1 1 2]));
pr1 = FreePropagation(pr0,z1,lambda,kx,ky);
pr2 = LensPropagation(pr1,f1,lambda,x,y);
pr3 = FreePropagation(pr2,z2a,lambda,kx,ky);
pr4 = imresize(pr3, (b*dx)/(lambda*f1));
pr5 = padarray(pr4,(s-size(pr4))/2,0,'both');

%% Reconstruction

c = s/2;
params.iterations = iterations;
params.eps = 0.05;
params.pixel_pitch = pixel_pitch;
params.f2 = f2;
params.lambda = lambda;

[object_rec, probe_rec, err] = ePIE_reconstruction(probes, cen_obj, pr5, c, params, obj);

figure('position',get(0,'ScreenSize')) 
subplot(2,3,1);
imshow(abs(object_rec), []), colorbar, title('Object Magnitude','FontSize',16); impixelinfo; % show reconstructed amplitude-image
subplot(2,3,2);
imshow(abs(angle(object_rec)), []), colorbar, title('Object Phase (absolut value)','FontSize',16); impixelinfo; % show reconstructed phase-image
subplot(2,3,4);
imshow(abs(probe_rec), []), colorbar, title('Probe Magnitude','FontSize',16); impixelinfo; % show reconstructed amplitude-probe
subplot(2,3,5);
imshow(angle(probe_rec), []), colorbar, title('Probe Phase','FontSize',16); impixelinfo; % show reconstructed phase-probe
subplot(2,3,3);
plot(err), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error'); % show Error VS iterations
subplot(2,3,6);
scatter(reshape(cen_obj(:,:,2),[1,n_row*n_col]),reshape(cen_obj(:,:,1),[1,n_row*n_col]),'x'); axis([1 s(2) 1 s(1)]); %show centers position
sgtitle('Reconstruction of object and Probe, amplitude and phase:','FontSize',24);
impixelinfo;

%% Help functions 

function probes = cut_probes1(image, n, margins, support)
% function for cutting the probes on the camera image into set of probes
% for SSP.
% input:
% image - camera intensity pattern ([NrowxNcol] complex)
% n - number of spots/cut ([1x2] double)
% output:
% probes - a 4D array with probes intensity on camera plane content by its 
%          order ([Nrow/n(1), Ncol/n(2), n(1), n(2)] complex)

if nargin == 2 
    margins = 0;
end

if margins >1          %if in percentage
    margins = margins/100;
end

image = rot90(image,2); % flip to arrang the centers according to probes
image = abs(image).^2; % taking only intensity

s = size(image)./n;
for i=1:n(1)
    for j=1:n(2)
        sub_image = image(floor(1+(i-1)*s(1):i*s(1)),floor(1+(j-1)*s(2):j*s(2)));
        probes(:,:,n(1)+1-i,n(2)+1-j) = sub_image(floor(1+s(1)*margins/2:s(1)*(1-margins/2)),floor(1+s(2)*margins/2:s(2)*(1-margins/2)));
    end
end

end

function l = draw_image_and_line(im,dist,l,xcor)
    subplot(1,2,1);axis tight; hold on;
    imagesc(abs(im));impixelinfo;
    title(['Field at z = ',num2str(dist),'[mm]']);
    subplot(1,2,2);  delete(l); l = line([1 1]*xcor,[1 257],'Color','green');
    drawnow;

end

function [imout1, imout2] = cut_min_ind(imin1,imin2)
    s1 = size(imin1); s2 = size(imin2);
    s = [min(s1(1),s2(1)), min(s1(2),s2(2))];
    imout1 = imin1(1:s(1),1:s(2)); imout2 = imin2(1:s(1),1:s(2));

end

function imout = normalize_image(imin)
    if max(abs(imin(:)))>min(abs(imin(:)))
        imout = (imin-min(imin(:)))/(max(imin(:))-min(imin(:)));
    else
        error('problems with normalizing image dynamic range');
    end
end

