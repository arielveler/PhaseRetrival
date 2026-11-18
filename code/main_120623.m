clear all; clc; close all;
addpath(genpath('code')); % add all subfolders of Master folder
addpath(genpath(pwd));                              % add all subfolders of current folder to path

% initials:
N = 2^13;
Nim = 440;                                            % Block size
lambda = 1030e-9;                                     %[m] central wavelenght
n = [5, 5];                                           % number of spots at the pinholes plane
theta = 0.97;                                         %[deg] seperation angle between spots (for the MS)
f1 = 0.1; f2 = 0.06; fms = 0.15;                      %[m] focal length of the lenses 
pixel_pitch = 3.45e-6;                                % [m]
W0_init = 0.00245;                                    % 1/e of radius of beam at entrance of the system
b = fms*tan(2*pi*theta/360);                          % distance between adjacent "pinholes"/probe-centers in pihole plane
sig_init = W0_init/sqrt(2);                           % sigma of beam at entrance of the system
Im_size = [N N];
Center = ceil((Im_size+1)./2);
M = f2/f1;                                            % enlarging
d = 0.015;

% definition of real grid:   
dx = pixel_pitch; dy = pixel_pitch;
x_vec = dx*(-N/2:N/2-1); y_vec = dy*(-N/2:N/2-1);
[x,y] = meshgrid(x_vec, y_vec);

% definition of reciprocoral grid:
dkx = 1/(N*dx); dky = 1/(N*dy); % [Hz] smallest element of frequency
kx_vec = (-N/2:N/2-1); ky_vec = (-N/2:N/2-1);
[kx,ky] = meshgrid(kx_vec, ky_vec);
kx = fftshift(kx)/(N*dx); ky = fftshift(ky)/(N*dx);

% real grid after scaling by the f2 lens, with block number of pixels:
dk = pixel_pitch/(lambda*f2);
dx_im = 1/(Nim*dk);
x_vec_im = dx_im*(-Nim/2:Nim/2-1); y_vec_im = dx_im*(-Nim/2:Nim/2-1);
[x_im,y_im] = meshgrid(x_vec_im, y_vec_im);

% filePath = dir(matlab.desktop.editor.getActiveFilename);
% fullfile = filePath.folder;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fullfile = "code\Figures";
%% set a simulated probe
load(strcat(fullfile,"\sim_ms-5x5_fms-150_f1-100_f2-60.mat"));
[centers_x,centers_y] = meshgrid(linspace(Center(2)-b*(n(2)-1)/2/pixel_pitch,Center(2)+b*(n(2)-1)/2/pixel_pitch,n(2)),linspace(Center(1)-b*(n(1)-1)/2/pixel_pitch,Center(1)+b*(n(1)-1)/2/pixel_pitch,n(1)));
centers = cat(3,centers_y,centers_x);
subimage_rec = Nim*[1 1]; % take only square blocks, for comfort

%% load first guess for probes:
probe_rec_1 = load(strcat(fullfile,"\reconstructions_010922\probe_no1_reconstruction.mat")); probe_rec_1 = probe_rec_1.probe_rec_1;
probe_rec_2 = load(strcat(fullfile,"\reconstructions_010922\probe_no2_reconstruction.mat")); probe_rec_2 = probe_rec_2.probe_rec_2;
probe_rec_3 = load(strcat(fullfile,"\reconstructions_010922\probe_no3_reconstruction.mat")); probe_rec_3 = probe_rec_3.probe_rec_3;
probe_guess(:,:,1) = probe_rec_1; probe_guess(:,:,2) = probe_rec_2; probe_guess(:,:,3) = probe_rec_3; % this specific case....

%%
%%%%%%%%%%%%%%%%%%%%%% SSP Reconstruction of reference frames %%%%%%%%%%%%%%%%%%%%%% 

%% 1. sp_holes: reconstruction of SSP of individual fibers after shouting the pumps:
imcen = load(strcat(fullfile,"\Fourier_plane\centers.mat")); imcen = double(imcen.frame);
imdark = load(strcat(fullfile,"\Fourier_plane\dark.mat")); imdark = double(imdark.frame);
imcen_n = noise_removal(imcen,imdark); % remove noise
% imcen_n = padarray(imcen_n,[10 10],0,"both");

[sim_cam_comp,imcen_n]=CutMinInd(abs(sim_ms),imcen_n);
[cen_rec,~,~,bt] = FindRefCenters(imcen_n,n,10,30);
[probes_ref,~] = cut_probes(imcen_n,cen_rec(end:-1:1,end:-1:1,:),subimage_rec); % cut the image without object for registration and energy normalization 
[cen_sim,~,~,bt_sim] = FindRefCenters(sim_cam_comp,n,10,30);
[probes_sim,~] = cut_probes(sim_cam_comp,cen_sim(end:-1:1,end:-1:1,:),subimage_rec); % cut the image without object for registration and energy normalization 

[xx,yy] = meshgrid(-Nim/2:Nim/2-1,-Nim/2:Nim/2-1);

distances = [0.010 0.011 0.012];
thresh = [0 0.01 0.05 0.1 0.15 0.2 0.3 0.5];

%%
% load reference and dark image:
% fib = [2 3 1];
count = 0;
fib = 1;
for ii=[0.5 1 2 3 4 5]
for dist = 1:1
    for thr = 1:1

        imrec = load(strcat(fullfile,"\1+2+3\att3\probe3_after_gain7.mat")); imrec = imgaussfilt(double(imrec.frame),ii,"FilterSize",[9 9]);
        imdark =load(strcat(fullfile,"\1+2+3\att3\dark7.mat")); imdark = imgaussfilt(double(imdark.frame),ii,"FilterSize",[9 9]);
        imrec_n = noise_removal(imrec,imdark,'thresholding',thresh(thr)); % remove noise
        d=distances(dist);
        probe_init = probe_guess(:,:,fib);

        % registration of the measured image:
        [probes_real,~] = cut_probes(imrec_n,cen_rec(end:-1:1,end:-1:1,:),subimage_rec); % cut the probes of measured image according to the reference image. invert centers order for reconstruction - top left block is suitable for bottom right probe

        % get parameters for measured image registration, normalize and register each block independetly:
        output_ref = zeros(n(1),n(2),4);
        for i=1:n(1)
            for j=1:n(2)
                [output_ref(i,j,:), probe_tmp] = dftregistration(fft2(probes_sim(:,:,i,j)),fft2(probes_ref(:,:,i,j)),10); %produce block by block registration, due to small errors in Z-axis of the camera (so the centers are not exactly at Fourier plane)
                probes_ref_reg(:,:,i,j) = abs(ifft2(probe_tmp));
                probes_real(:,:,i,j) = normalize_energy(probes_real(:,:,i,j),probes_ref(:,:,i,j)); % normalize energy
                probes_real_reg(:,:,i,j) = abs(ShiftProbefft(probes_real(:,:,i,j),[output_ref(i,j,3) output_ref(i,j,4)])); %  block-block registration
            end
        end

        % compute centers in object plane:
        s_cam = size(probes_real_reg,[1 2]); % size of the blocks
        cen_obj = (cen_rec-repmat(reshape(size(imcen_n)./2,[1 1 2]),n)).*(d/f2)*(dy^2*subimage_rec(1))/(lambda*f2)+repmat(ceil(s_cam(1)/2),n)+1; % centes in object plane in size of the probe image

        %---------------------  END OF MEASURED IMAGE PART -----------------------------%

        % Reconstruction
        clear params;
        params.iterations = 201;
        params.pixel_pitch = [1 1]*pixel_pitch;
        params.f1 = f1;
        params.f2 = f2;
        params.b = M*b;
        params.lambda = lambda;
        params.probe_init_pos = floor(s_cam/2)+1;
        params.show_internal = 1;
        params.change_iter = 30;
        params.is_support = 1;
        params.support_thr = 0.02;

        [object_rec, probe_rec, err] = ePIE_reconstruction(probes_real_reg, cen_obj, probe_init/1000, params); % reconstruct simulation data

        count = count+1;
        object_rec_mat(:,:,count) = object_rec;
        probe_rec_mat(:,:,count) = probe_rec;

    end
end
end
%%
for i=1:count
    figure; phplot(object_rec_mat(180:260,180:260,i),1,1); xlabel('mm'); ylabel('mm'); colorbar; colormap('jet'); title(['d= ' num2str(distances(ceil(i/4))*1000) 'mm']);
end
%% for drawing
figure;
a = numel(distances); b = numel(thresh);
for i=1:a
    for j=1:b
        subplot(a,b,(i-1)*b+j)
        imagesc(x_vec_im(180:260),x_vec_im(180:260),abs(object_rec_mat(180:260,180:260,(i-1)*b+j))); title(strcat("d=",num2str(distances(i)),"-thr=",num2str(thresh(j)))); axis square; colorbar; 
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%% TIMP Reconstruction of 3 frames %%%%%%%%%%%%%%%%%%%%%% 

%% 2. sp_holes: TIMP reconstruction of holes
imcen = load(strcat(fullfile,"\1+2+3\centers.mat")); imcen = double(imcen.frame);
imdark = load(strcat(fullfile,"\1+2+3\dark.mat")); imdark = double(imdark.frame);
imcen_n = noise_removal(imcen,imdark); % remove noise

[sim_cam_comp,imcen_n]=CutMinInd(abs(sim_ms),imcen_n);
[cen_rec,~,~,bt] = FindRefCenters(imcen_n,n,10,30);
[probes_ref,~] = cut_probes(imcen_n,cen_rec(end:-1:1,end:-1:1,:),subimage_rec); % cut the image without object for registration and energy normalization 
[cen_sim,~,~,bt_sim] = FindRefCenters(sim_cam_comp,n,10,30);
[probes_sim,~] = cut_probes(sim_cam_comp,cen_sim(end:-1:1,end:-1:1,:),subimage_rec); % cut the image without object for registration and energy normalization 

distances = [0.009 0.010 0.011];
thresh = [0.01 0.05 0.15 0.18];
% obj guess for the 3 objects:
[xx,yy] = meshgrid(-Nim/2:Nim/2-1,-Nim/2:Nim/2-1);
obj_guess = cat(3,2*xx/Nim,zeros(Nim)+1,2*yy/Nim)*10*eps;
ind = 180:260;
%% first reconstruction using higher threshold level - thr=0.18
% save_file = strcat(fullfile,'\1+2+3\att3\');
count=0;
fib = [2 3 1];
for ii=1:10
for dist = 2:2
    for thr = 2:2

        % load reference and dark image:
        imrec = load(strcat(fullfile,"\1+2+3\att3\probe1+2+3_during_gain24.mat")); sat_mask = imrec.frame >= 4095; imrec = imgaussfilt(double(imrec.frame),0.6,"FilterSize",[9 9]);
        imdark24 =load(strcat(fullfile,"\1+2+3\att3\dark24.mat")); imdark24 = imgaussfilt(double(imdark24.frame),0.6,"FilterSize",[9 9]);
        imrec_n = noise_removal(imrec,imdark,'thresholding',0.18); % remove noise
        probe_init = probe_guess;
        d=distances(dist);

        % registration of the measured image:
        [probes_real,~] = cut_probes(imrec_n,cen_rec(end:-1:1,end:-1:1,:),subimage_rec); % cut the probes of measured image according to the reference image. invert centers order for reconstruction - top left block is suitable for bottom right probe
        [sat_mask_real,~] = cut_probes(sat_mask,cen_rec(end:-1:1,end:-1:1,:),subimage_rec); % cut the probes of measured image according to the reference image. invert centers order for reconstruction - top left block is suitable for bottom right probe

        % get parameters for measured image registration, normalize and register each block independetly:
        output_ref = zeros(n(1),n(2),4);
        for i=1:n(1)
            for j=1:n(2)
                [output_ref(i,j,:), probe_tmp] = dftregistration(fft2(probes_sim(:,:,i,j)),fft2(probes_ref(:,:,i,j)),10); %produce block by block registration, due to small errors in Z-axis of the camera (so the centers are not exactly at Fourier plane)
                probes_ref_reg(:,:,i,j) = abs(ifft2(probe_tmp));
                probes_real(:,:,i,j) = normalize_energy(probes_real(:,:,i,j),probes_ref(:,:,i,j)); % normalize energy
                probes_real_reg(:,:,i,j) = abs(ShiftProbefft(probes_real(:,:,i,j),[output_ref(i,j,3) output_ref(i,j,4)])); %  block-block registration
                sat_mask_real_reg(:,:,i,j) = abs(ShiftProbefft(sat_mask_real(:,:,i,j),[output_ref(i,j,3) output_ref(i,j,4)]))>0.5; %  block-block registration of saturation mask
            end
        end

        % compute centers in object plane:
        s_cam = size(probes_real_reg,[1 2]); % size of the blocks
        cen_obj = (cen_rec-repmat(reshape(size(imcen_n)./2,[1 1 2]),n)).*(d/f2)*(dy^2*subimage_rec(1))/(lambda*f2)+repmat(ceil(s_cam(1)/2),n)+1; % centes in object plane in size of the probe image
        %---------------------  END OF MEASURED IMAGE PART -----------------------------%
        % %%
        % draw_blocks(probes_real_reg)
        % %%
        % Reconstruction
        clear params;
        params.iterations = 101;
        params.pixel_pitch = [1 1]*pixel_pitch;
        params.f1 = f1;
        params.f2 = f2;
        params.b = M*b;
        params.lambda = lambda;
        params.probe_init_pos = floor(s_cam/2)+1;
        params.show_internal = 1;
        params.change_iter = 41;
        params.is_support = 1;
        params.support_thr = 0.02;
        params.obj_guess = obj_guess(:,:,fib);
        params.sat_mask = sat_mask_real_reg;

        [object_rec, probe_rec, err] = ePIE_reconstruction_TIMP(probes_real_reg, cen_obj, probe_init(:,:,fib)/50, params); % reconstruct simulation data
        count = count+1;
        object_rec_mat(:,:,:,count) = object_rec;
    end
end
end
% save(strcat(save_file,"object_rec_d=10mm_thresh=[0_0.01_0.03_0.04_0.05_0.06_0.07].mat"),"object_rec_mat");
% save(strcat("C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\new_reconstructions_no_saturations_051023\","object_rec_d=10mm_thresh=0.18_10attempts.mat"),"object_rec_mat");
%% second reconstruction using one of previous higher threshold reconstruction as first guess, and low reconstruction to emphasis small features (thr=0.0, 0.01, 0.02) - several attempts each.
% load('C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\object_rec_d=10mm_thresh=0.18_filt=0.6_9x9_GT.mat')
load(strcat(fullfile,"\1+2+3\att3\new_reconstructions_no_saturations_051023_object_rec_d=10mm_thresh=0.18_10attempts.mat"))
params.obj_guess = object_rec_mat(:,:,:,5);
params.iterations = 201;
clear object_rec_mat probe_rec_mat

count = 0;
thr = [0 0 0 0 0 0.01 0.01 0.01 0.01 0.01 0.02 0.02 0.02 0.02 0.02];
for ii=1:numel(thr)
imrec = load(strcat(fullfile,"\1+2+3\att3\probe1+2+3_during_gain24.mat")); imrec = imgaussfilt(double(imrec.frame),0.6,"FilterSize",[9 9]);
imdark24 =load(strcat(fullfile,"\1+2+3\att3\dark24.mat")); imdark24 = imgaussfilt(double(imdark24.frame),0.6,"FilterSize",[9 9]);
imrec_n = noise_removal(imrec,imdark24,'thresholding',thr(ii)); % remove noise

% registration of the measured image:
[probes_real,~] = cut_probes(imrec_n,cen_rec(end:-1:1,end:-1:1,:),subimage_rec); % cut the probes of measured image according to the reference image. invert centers order for reconstruction - top left block is suitable for bottom right probe

% get parameters for measured image registration, normalize and register each block independetly:
output_ref = zeros(n(1),n(2),4);
for i=1:n(1)
    for j=1:n(2)
        [output_ref(i,j,:), probe_tmp] = dftregistration(fft2(probes_sim(:,:,i,j)),fft2(probes_ref(:,:,i,j)),10); %produce block by block registration, due to small errors in Z-axis of the camera (so the centers are not exactly at Fourier plane)
        probes_ref_reg(:,:,i,j) = abs(ifft2(probe_tmp));
        probes_real(:,:,i,j) = normalize_energy(probes_real(:,:,i,j),probes_ref(:,:,i,j)); % normalize energy
        probes_real_reg(:,:,i,j) = abs(ShiftProbefft(probes_real(:,:,i,j),[output_ref(i,j,3) output_ref(i,j,4)])); %  block-block registration
    end
end

[object_rec, probe_rec, err] = ePIE_reconstruction_TIMP(probes_real_reg, cen_obj, probe_init(:,:,fib)/50, params); % reconstruct simulation data
count = count+1;
object_rec_mat(:,:,:,count) = object_rec;
probe_rec_mat(:,:,:,count) = probe_rec;
end

%% Help functions 

function imout = normalize_energy(imin,imref)
% function for normalize a block of measurement of an object according to
% its initial energy (same block of measurement without the object).

    if ~exist('imref', 'var') % normalize itself - without object
        imref = imin;
    end
    if any(~isreal(imin(:))) || any(~isreal(imref(:)))
        error('image to normalize or reference has to have real values');
    end
    % normalize according to norm2:
    imout = imin./(sqrt(sum(imref(:).^2)));
end


