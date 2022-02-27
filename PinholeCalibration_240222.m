clear all; clc; close all;
addpath(genpath(pwd)); % add all subfolders of current folder to path
fullfile = "C:\Users\arielv\Desktop\Master\images\pinhole_240222";
im = load(strcat(fullfile,"\pinhole_0mm.mat")); im = im.frame;
[cen,subimage,theta,bt] = FindRefCenters(im,[5 5],10,30);

% initials:
Nim = size(im);
lambda = 1030e-9; %[m] central wavelenghtimagesc(im)

n_row = 5; n_col = 5; n = [n_row, n_col]; % number of spots at the pinholes plane
theta = [0.97 0.97]; %[deg] seperation angle between spots (for the MS)
pixel_pitch = 3.45e-6; % [m]
f1 = 100e-3; f2 = 100e-3; %[m] focal length of the lenses 
d = 1.5e-3; % pinhole diameter

%% definitions of grid:   

% K-grid - camera plane:

k0 = 2*pi/lambda;
sensor = Nim.*pixel_pitch; % sensor physical size

N = 10*round(mean(subimage)/10); % round subimage up to 10 pixels
dk = pixel_pitch/(lambda*f2); % frequency in Fourier plane (=sensor plane) is nu=x/(lam*f) im [1/m]
kx = dk*(-N/2:N/2-1); % grid of frequencies according to image # of pixels. also (-1) becase of harmonic functions [1/m]
ky = dk*(-N/2:N/2-1); % [1/m]
[KX, KY] = meshgrid(kx,ky);

% X-grid - object plane:

dx = 1/(N*dk); % Smallest element on object plane corresponds to the largest element in image plane
x = dx*(-N/2:N/2-1); % same amount of pixels in FFT2
y = dx*(-N/2:N/2-1);
[X ,Y] = meshgrid(x,y);

% ax = sensor(2); ay = sensor(1); % pinhol plane physical size
% dx = ax/N(1); dy = ay/N(2);
% [x,y] = meshgrid(-ax/2:dx:ax/2-dx, -ay/2:dy:ay/2-dy);
% 
% Wk0x = 2*pi/dx; Wk0y = 2*pi/dy; % [Hz] largest element of frequencies for fourier transform
% dkx = Wk0x/N(1); dky = Wk0y/N(2); % [Hz] smallest element of frequency
% [kx,ky] = meshgrid(-Wk0x/2:dkx:Wk0x/2-dkx,-Wk0y/2:dky:Wk0y/2-dky);
[probes,~] = cut_probes(im,cen,[N N]);

%% Compare each probe detected to a pinhole FFT

im_sim = double(sqrt(X.^2+Y.^2)<(d/2));
imagesc(x,y,im_sim);impixelinfo;
im_sim_k = abs(fftshift(fft2(im_sim))).^2;
im_sim_k = im_sim_k./max(im_sim_k(:)); % normalize
rmse = zeros(n);
for i=1:n(1)
    for j=1:n(2)
        im_tmp = double(probes(:,:,i,j));
%         figure(1); imagesc(kx,ky,im_tmp);impixelinfo; title(['recorded probe no. ',num2str(i),'-',num2str(j)]);
%         axis square
        noise = mean(mean([im_tmp(1:10,1:10);im_tmp(1:10,end-9:end);im_tmp(end-9:end,1:10);im_tmp(end-9:end,end-9:end)]))*2;
        im_tmp = im_tmp-noise; % removing BG noise
        im_tmp = im_tmp./max(im_tmp(:)); % normalizing probe image
%         figure(2); imshowpair(im_sim_k,im_tmp,'montage');impixelinfo;
        figure(3); subplot(n(1),n(2),i+n(2)*(j-1)); plot(kx,im_tmp(251,:),kx,im_sim_k(251,:)); %title(['comparison between simulated and recorded probe no. ',num2str(i),'-',num2str(j)])
        title(['no. ',num2str(i),'-',num2str(j)]);
        sgtitle('Comparison between simulated and recorded probe');
        rmse(i,j) = sqrt(mean((im_tmp(251,:)-im_sim_k(251,:)).^2));
    end
end

% don't forget to remove background and normalize!!!!!!!!
%% calculate the distance from fourier plane:
% p = RoundProbe(d,Nim,Nim/2);
% b = p.find_next_center(0.7);
% find_b = load(strcat(fullfile,"\find_b.mat")); find_b = find_b.frame;
% [cen,~,~,bt] = FindRefCenters(find_b,[5 5],10,15);
% bt = mean(bt(:));
% d_from_forier = f1*b/(bt*pixel_pitch);

%% Reconstruction of 32 mm (~70% overlap)
imd = load(strcat(fullfile,"\pinhole_32mm.mat")); imd = imd.frame;
imd = rot90(im,2); % flip image to be in the rigth orientation with the object 
Dpix = 30;
Thr = 10;
[cen,subimage,theta,bt] = FindRefCenters(imd,n,Thr,Dpix);

N = 10*round(mean(subimage)/10); % round subimage
bt = mean(bt.*pixel_pitch);
[probes,~] = cut_probes(imd,cen,[N N]); probes = double(probes);
noise = squeeze(mean([probes(1:10,1:10,:,:);probes(1:10,end-9:end,:,:);probes(end-9:end,1:10,:,:);probes(end-9:end,end-9:end,:,:)],[1 2])*2);
probes = probes-repmat(reshape(noise,[1 1 n(1) n(2)]),[N,N,1,1]); % remove BG noise
probes = probes./max(probes,[],[1 2]); % normalize the image

mid_cen = mean(cen(floor((n(1)+1)/2):ceil((n(1)+1)/2),floor((n(2)+1)/2):ceil((n(2)+1)/2),:),[1 2]); % center of image, point [0 0] 
factor = (dx*bt)/(lambda*f2);
cen_obj = (cen-repmat(mid_cen,n)).*pixel_pitch.*(d./f2).*factor+repmat(reshape([N N]./2,[1 1 2]),n); % centes in object plane in size of the probe image
% cen_obj = (cen-repmat(mid_cen,n)).*(d./f2)+repmat(reshape(subimage./2,[1 1 2]),n); % centes in object plane in size of the probe image
% cen_obj(:,:,1) = (cen_obj(:,:,1)-repmat(subimage(1)/2,n))*(bt)/(lambda*f2)+repmat(subimage(1)/2,n); % centes in object plane, when image is in size of cutted image on camera plane

% c = round(reshape(cen(ceil(n(1)/2),ceil(n(2)/2),:),[1 2]));
probe_guess = im_sim;
% probe_guess = im(c(1)-floor(subimage(1)/2):c(1)+floor(subimage(1)/2),c(2)-floor(subimage(2)/2):c(2)+floor(subimage(2)/2));
% probe_guess = imresize(probe_guess, (bt*pixel_pitch(1))/(lambda*f1));
% probe_guess = padarray(probe_guess,(subimage-size(probe_guess))/2,0,'both');

%%
%%%%%%%%%%%%%%%%%% END EDITING %%%%%%%%%%%%%%%%%%%

params.iterations = 50;
params.eps = 0.05;
params.pixel_pitch = [1 1]*pixel_pitch;
params.f2 = f2;
params.lambda = lambda;

[object_rec, probe_rec, err] = ePIE_reconstruction(probes, cen_obj, probe_guess, ceil((subimage-1)/2), params);

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