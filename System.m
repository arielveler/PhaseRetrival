%% BPM symulation - Ariel Veler

clear all; clc; close all;
addpath(genpath(pwd)); % add all subfolders of current folder to path

% initials:
N = 8192;
lambda = 1030e-9; %[m] central wavelenght
n_row = 5; n_col = 5; n = [n_row, n_col]; % number of spots at the pinholes plane
theta = [0.97 0.97]; %[deg] seperation angle between spots (for the MS)
pixel_pitch = 3.45e-6; % [m]
f1 = 200e-3; f2 = 100e-3; %[m] focal length of the lenses 
% D = 120e-6; % cut in FWHM of the probe
W0 = 0.001;
b = 3.4e-3; % distance between pinholes
d = 20e-3;
fms = 120e-3;
resolution = 10e-6; %resolution desired on object

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sig = D./(2*(2*log(2))^0.5); % FWHM to sigma
sig = lambda*fms/(pi*W0)/sqrt(2); % FWHM to sigma
sig_pix = sig./pixel_pitch; %for gauss probe
Im_size = [N N];
Center = Im_size./2;
M = f2./f1; % enlarging
FOV = n*b*d/f1; % FOV estimation on the object plane
% overlap = 1-(pi*D*b*d)/(lambda*sqrt(1+d*(lambda*D^2)/(pi*f1^2))); %overlap estimation between probes on the object plane

% definition of real grid:   

k0 = 2*pi/lambda;
sensor = Im_size.*pixel_pitch; % sensor physical size

dx = pixel_pitch; dy = pixel_pitch;
x_vec = -sensor(2)/2:dx:sensor(2)/2-dx; y_vec = -sensor(1)/2:dy:sensor(1)/2-dy;
[x,y] = meshgrid(x_vec, y_vec);

% kmax_x = 2*pi/dx; kmax_y = 2*pi/dy; % [Hz] largest element of frequencies for fourier transform
% dkx = kmax_x/Ncol; dky = kmax_y/Nrow; % [Hz] smallest element of frequency
% kx_vec = -kmax_x/2:dkx:kmax_x/2-dkx; ky_vec = -kmax_y/2:dky:kmax_y/2-dky;
% [kx,ky] = meshgrid(kx_vec, ky_vec);

% definition of reciprocoral grid:

dkx = 2*pi/(N*dx); dky = 2*pi/(N*dy); % [Hz] smallest element of frequency
kmax_x = dkx*N; kmax_y = dky*N; % [Hz] largest element of frequencies for fourier transform
kx_vec = -kmax_x/2:dkx:kmax_x/2-dkx; ky_vec = -kmax_y/2:dky:kmax_y/2-dky;
[kx,ky] = meshgrid(kx_vec, ky_vec);
%% settings
% initial wave at z=0

p = GaussProbe(sig, pixel_pitch*[1 1], lambda, Im_size, Center);

[centers_x,centers_y] = meshgrid(round(linspace(Center(2)-b*(n(2)-1)/2/pixel_pitch,Center(2)+b*(n(2)-1)/2/pixel_pitch,n(2))),round(linspace(Center(1)-b*(n(1)-1)/2/pixel_pitch,Center(1)+b*(n(1)-1)/2/pixel_pitch,n(1))));
centers = cat(3,centers_y,centers_x);
[im0,cen,sub_size] = p.pinhols(n,centers);

%% object

% load object:
% amp = double(imread('trees.tif')); % double(imread('liftingbody.png'));
% phase = double(imread('cameraman.tif')); % double(imread('westconcordorthophoto.png'));
Create_numbers_image();
amp = imresize(a7,256/size(a7,1));
phase = imresize(a8,256/size(a8,1));

[amp,phase] = cut_min_ind(amp,phase);
amp = normalize_image(amp); phase = normalize_image(phase);

obj = amp.*exp(1i.*phase);
obj = imresize(obj,4);
obj = padarray(obj,(Im_size-size(obj))/2,1,'both');


%% Propagate probes
im1 = FreePropagation(im0,f1,lambda,kx,ky);
im2 = LensPropagation(im1,f1,lambda,x,y);
im3 = FreePropagation(im2,f1-d,lambda,kx,ky);
im4 = im3.*obj;
im5 = FreePropagation(im4,f2+d,lambda,kx,ky);
im6 = LensPropagation(im5,f2,lambda,x,y);
im7 = FreePropagation(im6,f2,lambda,kx,ky);
im8 = abs(im7).^2; % taking only intensity
im9 = rot90(im8,2); % flip to arrang the centers according to probes

mid_cen = mean(cen(floor((n(1)+1)/2):ceil((n(1)+1)/2),floor((n(2)+1)/2):ceil((n(2)+1)/2),:),[1 2]); % center of image, point [0 0] 
cen = (cen-repmat(mid_cen,n))*M+repmat(mid_cen,n);

% cutting probes from the image. result in intensity units(field magnitude square)
[probes,s] = cut_probes(im9,cen);
probes = probes./max(probes,[],[1 2]); % normalize the images

factor = (dx*b*M)/(lambda*f2);
cen_obj = (cen-repmat(mid_cen,n)).*M.*(d./f2).*factor+repmat(reshape(s./2,[1 1 2]),n); % centes in object plane in size of the probe image
% cen_obj = (cen-repmat(reshape(Im_size./2,[1 1 2]),n)).*(d./f1)+repmat(reshape(s./2,[1 1 2]),n); % centes in object plane in size of the probe image
% cen_obj(:,:,1) = (cen_obj(:,:,1)-repmat(s(1)/2,n))*(M*b*dy)/(lambda*f2)+repmat(s(1)/2,n); % centes in object plane, when image is in size of cutted image on camera plane
% cen_obj(:,:,2) = (cen_obj(:,:,2)-repmat(s(2)/2,n))*(M*b*dx)/(lambda*f2)+repmat(s(2)/2,n); % centes in object plane, when image is in size of cutted image on camera plane

%% calibrate initial guess of probe
[pr0,~,~] = p.pinhols([1 1],reshape(Center,[1 1 2]));
pr1 = FreePropagation(pr0,f1,lambda,kx,ky);
pr2 = LensPropagation(pr1,f1,lambda,x,y);
pr3 = FreePropagation(pr2,f1-d,lambda,kx,ky);
pr4 = FreePropagation(pr3,d+f2,lambda,kx,ky);
pr5 = LensPropagation(pr4,f2,lambda,x,y);
pr6 = FreePropagation(pr5,f2,lambda,kx,ky);
pr7 = rot90(abs(pr6).^2,2);
pr8 = imresize(pr7, factor);
% if any(s>size(pr8))
%     pr9 = padarray(pr8,(s-size(pr8))/2,0,'both');
% else
%     pr9 = cut_probes(pr8,reshape(size(pr8)/2,[1 1 2]),s);
% end
pr_tmp = imresize(pr7, 3);
pr9 = cut_probes(pr_tmp,reshape(size(pr_tmp)/2,[1 1 2]),s);
%% Reconstruction

c = s/2;
iterations = 50;
params.pixel_pitch = pixel_pitch*[1 1];
params.f2 = f2;
params.lambda = lambda;

%%%%%%%%%%%%%%%%%% END EDITING %%%%%%%%%%%%%%%%%%%

% params.sig = sig_new;
% % params.sig = sig.*(1+rand(1)/3); % put an error in sigma
% params.center = c;

% params.ProbeType = 'GaussProbe';
% params.LoadProbe = 'GaussProbe(params.sig,params.pixel_pitch,params.lambda,params.im_size,params.center);';
params.iterations = iterations;
params.eps = 0.0001;

[object_rec, probe_rec, err] = ePIE_reconstruction(probes, cen_obj, pr9, c, params, obj);

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

%% visual propagation
% [Eout,imnumout] = VisualPropagation(type,params,zin,steps,imnumin)
to_save = 0; %save figure to create video or just render images to the screen?
if to_save
    figure('position',get(0,'ScreenSize'),'visible','off');axis tight; hold on;
else
    figure('position',get(0,'ScreenSize'));axis tight; hold on;
end
    
imnum = to_save;
z1 = f1;
params.Ein = im0; params.z = z1; params.lambda = lambda; params.k_x = kx; params.k_y = ky;
[im1,imnum] = VisualPropagation('Free',params,0,10,imnum);
params.Ein = im1; params.f = f1; params.lambda = lambda; params.x = x; params.y = y;
[im2,imnum] = VisualPropagation('Lens',params,z1,0,imnum);
z2a = f1-d;
params.Ein = im2; params.z = z2a; params.lambda = lambda; params.k_x = kx; params.k_y = ky;
[im3,imnum] = VisualPropagation('Free',params,z1,10,imnum);
im4 = im3.*obj;
z2b = f2+d;
params.Ein = im4; params.z = z2b; params.lambda = lambda; params.k_x = kx; params.k_y = ky;
[im5,imnum] = VisualPropagation('Free',params,z1+z2a,10,imnum);
params.Ein = im5; params.f = f2; params.lambda = lambda; params.x = x; params.y = y;
[im6,imnum] = VisualPropagation('Lens',params,z1+z2a+z2b,0,imnum);
z3 = f2;
params.Ein = im6; params.z = z3; params.lambda = lambda; params.k_x = kx; params.k_y = ky;
[im7,imnum] = VisualPropagation('Free',params,z1+z2a+z2b,10,imnum);

if to_save
    cd('.\movie\');
    MakeVideo(pwd,5)
end

%%
system = imread('system4f.PNG');

figure('position', [200 200 1500 600]);
subplot(1,2,2); 
imshow(system,[]);
title('system along "z" axis');
xcor = round([88:12.2:210,213:11.95:452,455:12.2:577]);
j=1;

subplot(1,2,1);axis tight; hold on;

% 1. free propagation befor lens (z=0[mm]-100[mm])
z1 = 100e-3;
steps = z1*1e3/10;
subplot(1,2,1);axis tight; hold on;
imagesc(abs(im0));impixelinfo;
subplot(1,2,2); l = line([1 1]*xcor(j),[1 257],'Color','green'); j=j+1;
% saveas(gcf,['.\movie\im' repelem(char(floor((j-2)/5)+65),mod((j-2),5)+1) ,'.png']);
for i=1:steps
    im1 = FreePropagation(im0,z1*(i/steps),lambda,kx,ky);
    l = draw_image_and_line(im1,z1*(i/steps)*1e3,l,xcor(j)); j=j+1;
%     saveas(gcf,['.\movie\im' repelem(char(floor((j-2)/5)+65),mod((j-2),5)+1) ,'.png']);
end

% 2. propagation throuth lens #1 (z=100[mm])
im2 = LensPropagation(im1,f1,lambda,x,y);
l = draw_image_and_line(im2,z1*1e3,l,xcor(j)); j=j+1;
% saveas(gcf,['.\movie\im' repelem(char(floor((j-2)/5)+65),mod((j-2),5)+1) ,'.png']);

% 3. free propagation after lens 1 (z=100[mm]-300[mm])
% z2 = 200e-3;
% steps = z2*1e3/10;
% for i=1:steps
%     im3 = FreePropagation(im2,z2*(i/steps),lambda,kx,ky);
%     l = draw_image_and_line(im3,(z1+z2*(i/steps))*1e3,l,xcor(j)); j=j+1;
% %     saveas(gcf,['.\movie\im' repelem(char(floor((j-2)/5)+65),mod((j-2),5)+1) ,'.png']);
% end

% propagation throuth object
z2a = (100-d)*1e-3;
z2b = d*1e-3;
steps = 9;
for i=1:steps
    im3a = FreePropagation(im2,z2a*(i/steps),lambda,kx,ky);
    l = draw_image_and_line(im3a,(z1+z2a*(i/steps))*1e3,l,xcor(j)); j=j+1;
%     saveas(gcf,['.\movie\im' repelem(char(floor((j-2)/5)+65),mod((j-2),5)+1) ,'.png']);
end

%multiply with object:
im3b = im3a.*obj;

im3b = FreePropagation(im3b,z2b,lambda,x,y);
l = draw_image_and_line(im4,(z1+z2a+z2b)*1e3,l,xcor(j)); j=j+1;

% 4. propagation throuth lens #2 (z=300[mm])
im4 = LensPropagation(im3b,f2,lambda,x,y);
l = draw_image_and_line(im4,(z1+z2a+z2b)*1e3,l,xcor(j)); j=j+1;
% saveas(gcf,['.\movie\im' repelem(char(floor((j-2)/5)+65),mod((j-2),5)+1) ,'.png']);

% 5. free propagation after lens #2 (z=300[mm]-400[mm])
z3 = 100e-3;
steps = z3*1e3/10;
for i=1:steps
    im5 = FreePropagation(im4,z3*(i/steps),lambda,kx,ky);
    l = draw_image_and_line(im5,(z1+z2a+z2b+z3*(i/steps))*1e3,l,xcor(j)); j=j+1;
%     saveas(gcf,['.\movie\im' repelem(char(floor((j-2)/5)+65),mod((j-2),5)+1) ,'.png']);
end

% cd('C:\Users\arielv\Desktop\תואר שני\חניכה\phase retrival algorithms\movie')
% MakeVideo(pwd,5)

%% propagation througth a multispot and lens

p.center = Center;
im0 = p.probe;
steps = 5;
imagesc(abs(im0));impixelinfo;
z1 = 50e-3; %[m] propagation distance
for i=1:steps
    im1 = FreePropagation(im0,(z1)*i/steps,lambda,kx,ky);
    imagesc(abs(im1));impixelinfo;
    title(['Field at z = ',num2str(z1*(i/steps)*1e3),'[mm]']);
    drawnow;
end

im2 = MultispotPropagation(im1,theta,n,lambda,x,y);
imagesc(abs(im2));impixelinfo;
title(['Field at z = ',num2str(z1*1e3), '[mm] rigth after the multispot']);
drawnow;

z2 = 100e-3;
steps = 10;
for i=1:steps
    im3 = FreePropagation(im2,(z2)*i/steps,lambda,kx,ky);
    imagesc(abs(im3));impixelinfo;
    title(['Field at z = ',num2str((z2)*i/steps*1e3),'[mm] after the Multispot']);
    drawnow;
end

im4 = LensPropagation(im3,f1,lambda,x,y);
imagesc(abs(im4));impixelinfo;
title(['Field at z = ',num2str(z*1e3), '[mm] rigth after the lens']);
drawnow;

z3 = 100e-3;
for i=1:steps
    im5 = FreePropagation(im4,(z3)*i/steps,lambda,kx,ky);
    imagesc(abs(im5));impixelinfo;
    title(['Field at z = ',num2str((z3)*i/steps*1e3),'[mm] after the lens']);
    drawnow;
end


%% Plan an experiment - restrictions:
lambda = 1030e-6; %wavelength - must!
dx_f = 10e-6;
FOV = 100*dx_f;
dx_t = 3.45e-6;
alpha = 3;

%conditions:
f2 = alpha*FOV*dx_t/lambda;

f1_b = dx_f./lambda;




D = 25e-6; % pinhole width
f1 = 50e-3; f2 = 50e-3;
b = 10e-3;
d = 18e-3;
bt = b*d/f1;
n = [7 7];

% calculate d
W = (lambda*f1/(pi*D))*sqrt(1+d*pi*D^2/(lambda*f1^2));
overlap = (W/10-bt)/W/10

%% Help functions 

% function probes = cut_probes(image, n, margins, support)
% % function for cutting the probes on the camera image into set of probes
% % for SSP.
% % input:
% % image - camera intensity pattern ([NrowxNcol] complex)
% % n - number of spots/cut ([1x2] double)
% % output:
% % probes - a 4D array with probes intensity on camera plane content by its 
% %          order ([Nrow/n(1), Ncol/n(2), n(1), n(2)] complex)
% 
% if nargin == 2 
%     margins = 0;
% end
% 
% if margins >1          %if in percentage
%     margins = margins/100;
% end
% 
% image = rot90(image,2); % flip to arrang the centers according to probes
% image = abs(image).^2; % taking only intensity
% 
% s = size(image)./n;
% for i=1:n(1)
%     for j=1:n(2)
%         sub_image = image(floor(1+(i-1)*s(1):i*s(1)),floor(1+(j-1)*s(2):j*s(2)));
%         probes(:,:,n(1)+1-i,n(2)+1-j) = sub_image(floor(1+s(1)*margins/2:s(1)*(1-margins/2)),floor(1+s(2)*margins/2:s(2)*(1-margins/2)));
%     end
% end
% 
% end

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

