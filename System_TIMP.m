%% BPM symulation - Ariel Veler

clear all; clc; close all;
addpath(genpath(pwd)); % add all subfolders of current folder to path

% initials:
Nrow = 4096; Ncol = 4096;
lambda = 520e-9; %[m] central wavelenght
n_row = 9; n_col = 9; n = [n_row, n_col]; % number of spots at the pinholes plane
theta = [0.97 0.97]; %[deg] seperation angle between spots (for the MS)
pixel_pitch = [1 1].*3.45e-6; % [m]
f1 = 100e-3; f2 = 100e-3; %[m] focal length of the lenses 
D = 120e-6; % cut in FWHM of the probe
b = 1.4e-3; % distance between pinholes
d = 16e-3;
fmla = pi*b*D/4/lambda;
resolution = 10e-6; %resolution desired on object

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sig = 10.*min(pixel_pitch); %for gauss probe
sig = D./(2*(2*log(2))^0.5); % FWHM to sigma
sig_pix = sig./min(pixel_pitch); %for gauss probe
Im_size = [Nrow Ncol];
Center = Im_size./2;
M = f2./f1; % enlarging
FOV = n*b*d/f1; % FOV estimation on the object plane
overlap = 1-(pi*D*b*d)/(lambda*sqrt(1+d*(lambda*D^2)/(pi*f1^2))); %overlap estimation between probes on the object plane

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

[centers_x,centers_y] = meshgrid(round(linspace(Center(2)-b*(n(2)-1)/2/pixel_pitch(1),Center(2)+b*(n(2)-1)/2/pixel_pitch(1),n(2))),round(linspace(Center(1)-b*(n(1)-1)/2/pixel_pitch(2),Center(1)+b*(n(1)-1)/2/pixel_pitch(2),n(1))));
centers = cat(3,centers_y,centers_x);

% parameters of the probes
order1.l = 1; order1.m = 1; order1.polarization = 'x';
order2.l = 1; order2.m = 1; order2.polarization = 'y';
fiber_params.CoreRad = 100e-6; % 100 microns core radius 
fiber_params.CladRad = 250e-6; % 250 microns clad radius
fiber_params.lambda = lambda; % 0 nm central wavelength
fiber_params.n1 = refrIndex('SILICA', fiber_params.lambda); % Silica is core material
fiber_params.n2 = refrIndex('SILICACLAD', fiber_params.lambda);

probe1 = LPmodeProbe(order1,fiber_params,pixel_pitch,Im_size,Center);
probe2 = LPmodeProbe(order2,fiber_params,pixel_pitch,Im_size,Center);


[im1,cen,sub_size] = probe1.pinhols(n,centers);
disp('prob 1 pinholes is created');
[im2,~,~] = probe2.pinhols(n,centers);
disp('prob 2 pinholes is created');
%% objects

% load objects:
Create_numbers_image();

% object 1:
amp1 = imresize(a7,256/size(a7,1));
phase1 = imresize(a8,256/size(a8,1));
[amp1,phase1] = cut_min_ind(amp1,phase1);
amp1 = normalize_image(amp1); phase1 = normalize_image(phase1);

obj1 = amp1.*exp(1i.*phase1);
obj1 = imresize(obj1,4);
obj1 = padarray(obj1,(Im_size-size(obj1))/2,1,'both');

% object 2:
amp2 = imresize(a4,256/size(a4,1));
phase2 = imresize(a5,256/size(a5,1));
[amp2,phase2] = cut_min_ind(amp2,phase2);
amp2 = normalize_image(amp2); phase2 = normalize_image(phase2);

obj2 = amp2.*exp(1i.*phase2);
obj2 = imresize(obj2,4);
obj2 = padarray(obj2,(Im_size-size(obj2))/2,1,'both');

obj(:,:,1) = obj1; obj(:,:,2) = obj2;

sig_new = (lambda*f1)/(pi*sig);

%% Propagate objects
% propagate object 1:
z1 = f1;
im1 = FreePropagation(im1,z1,lambda,kx,ky);
im1 = LensPropagation(im1,f1,lambda,x,y);
z2a = f1-d;
im1 = FreePropagation(im1,z2a,lambda,kx,ky);
im1 = im1.*obj1;
z2b = f2+d;
im1 = FreePropagation(im1,z2b,lambda,kx,ky);
im1 = LensPropagation(im1,f2,lambda,x,y);
z3 = f2;
im1 = FreePropagation(im1,z3,lambda,kx,ky);
im1 = abs(im1).^2; % taking only intensity
disp('prob 1 is propagated');

% propagate object 2:
z1 = f1;
im2 = FreePropagation(im2,z1,lambda,kx,ky);
im2 = LensPropagation(im2,f1,lambda,x,y);
z2a = f1-d;
im2 = FreePropagation(im2,z2a,lambda,kx,ky);
im2 = im2.*obj2;
z2b = f2+d;
im2 = FreePropagation(im2,z2b,lambda,kx,ky);
im2 = LensPropagation(im2,f2,lambda,x,y);
z3 = f2;
im2 = FreePropagation(im2,z3,lambda,kx,ky);
im2 = abs(im2).^2; % taking only intensity
disp('prob 2 is propagated');

% collect both transmitted probes images in camera
im = im1 + im2; % camera integrate over all probes*objects.
im = rot90(im,2); % flip to arrang the centers according to probes


% cutting probes from the image. result in intensity units(field magnitude square)
% probes = cut_probes(im7,n,0); % image, no. of probes, percent of margins
for i=1:n(1)
    for j=1:n(2)
        probes(:,:,i,j) = im(cen(i,j,1)-round((sub_size(1)-1)/2):cen(i,j,1)+round((sub_size(1)-1)/2),cen(i,j,2)-round((sub_size(2)-1)/2):cen(i,j,2)+round((sub_size(2)-1)/2));
    end
end

s = size(probes,[1 2]);
cen_obj = (cen-repmat(reshape(Im_size./2,[1 1 2]),n)).*(d./f1)+repmat(reshape(s./2,[1 1 2]),n); % centes in object plane in size of the probe image
cen_obj(:,:,1) = (cen_obj(:,:,1)-repmat(s(1)/2,n))*(M*b*dy)/(lambda*f2)+repmat(s(1)/2,n); % centes in object plane, when image is in size of cutted image on camera plane
cen_obj(:,:,2) = (cen_obj(:,:,2)-repmat(s(2)/2,n))*(M*b*dx)/(lambda*f2)+repmat(s(2)/2,n); % centes in object plane, when image is in size of cutted image on camera plane

%% calibrate initial guess of probe
[pr1,~,~] = probe1.pinhols([1 1],reshape(Center,[1 1 2]));
pr1 = FreePropagation(pr1,z1,lambda,kx,ky);
pr1 = LensPropagation(pr1,f1,lambda,x,y);
pr1 = FreePropagation(pr1,z2a,lambda,kx,ky);
pr1 = imresize(pr1, (b*dx)/(lambda*f1));
pr1 = padarray(pr1,(s-size(pr1))/2,0,'both');
disp('prob 1 initial guess is propagated');

[pr2,~,~] = probe2.pinhols([1 1],reshape(Center,[1 1 2]));
pr2 = FreePropagation(pr2,z1,lambda,kx,ky);
pr2 = LensPropagation(pr2,f1,lambda,x,y);
pr2 = FreePropagation(pr2,z2a,lambda,kx,ky);
pr2 = imresize(pr2, (b*dx)/(lambda*f1));
pr2 = padarray(pr2,(s-size(pr2))/2,0,'both');
disp('prob 2 initial guess is propagated');

pr(:,:,1) = pr1; pr(:,:,2) = pr2;
%% Reconstruction

c = s/2;
iterations = 50;

%%%%%%%%%%%%%%%%%% END EDITING %%%%%%%%%%%%%%%%%%%

% params.sig = sig_new;
% % params.sig = sig.*(1+rand(1)/3); % put an error in sigma
% params.pixel_pitch = pixel_pitch;
% params.lambda = lambda;
% params.im_size = s;
% params.center = c;
% 
% params.ProbeType = 'GaussProbe';
% params.LoadProbe = 'GaussProbe(params.sig,params.pixel_pitch,params.lambda,params.im_size,params.center);';
params.iterations = iterations;

[object_rec, probe_rec, err] = ePIE_reconstruction_TIMP(probes, cen_obj, pr, c, obj, 0.0001,params);

for ii=1:size(pr,3)
    figure('position',[100 100 1200 600]); 
    subplot(2,3,1);
    imshow(abs(object_rec(:,:,ii)), []), colorbar, title(['Object-',num2str(ii),' Magnitude'],'FontSize',16); impixelinfo; % show reconstructed amplitude-image
    subplot(2,3,2);
    imshow(abs(angle(object_rec(:,:,ii))), []), colorbar, title(['Object-',num2str(ii),' Phase (absolut value)'],'FontSize',16); impixelinfo; % show reconstructed phase-image
    subplot(2,3,4);
    imshow(abs(probe_rec(:,:,ii)), []), colorbar, title(['Probe-',num2str(ii),' Magnitude'],'FontSize',16); impixelinfo; % show reconstructed amplitude-probe
    subplot(2,3,5);
    imshow(angle(probe_rec(:,:,ii)), []), colorbar, title(['Probe-',num2str(ii),' Phase'],'FontSize',16); impixelinfo; % show reconstructed phase-probe
    subplot(2,3,3);
    plot(err(:,ii)), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error'); % show Error VS iterations
    subplot(2,3,6);
    scatter(reshape(cen_obj(:,:,2),[1,n_row*n_col]),reshape(cen_obj(:,:,1),[1,n_row*n_col]),'x'); axis([1 s(2) 1 s(1)]); %show centers position
    sgtitle(['Reconstruction of object-',num2str(ii),' and Probe-',num2str(ii),', amplitude and phase:'],'FontSize',24);
    impixelinfo;
end

%%
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