%% Show mixing of image and phase
I1 = imread('cameraman.tif');
I2 = imread('trees.tif');
s = size(I1);
I2 = I2(1:s(1),1:s(2));

fft_I1 = fftshift(fft2(I1));
fft_I2 = fftshift(fft2(I2));
amp_fft_I1 = abs(fft_I1);
amp_fft_I2 = abs(fft_I2);
phase_fft_I1 = angle(fft_I1);
phase_fft_I2 = angle(fft_I2);

amp_I1_phase_I2 = amp_fft_I1.*exp(1i*phase_fft_I2);
amp_I2_phase_I1 = amp_fft_I2.*exp(1i*phase_fft_I1);

rec_I1_eiI2 = ifft2(ifftshift(amp_I1_phase_I2));
rec_I2_eiI1 = ifft2(ifftshift(amp_I2_phase_I1));
figure(1); 
%imshow(I1,[]);
%imshow(log(1+abs(fft_I1)),[]);
imshow(rec_I1_eiI2,[]); impixelinfo; title('Magnitude of cameraman with trees phase');
% 
figure(2); 
% %imshow(I2,[]);
% %imshow(log(1+abs(fft_I2)),[]);
imshow(rec_I2_eiI1,[]); impixelinfo; title('Magnitude of trees with cameraman phase');

%% GS reconstruction
clear; clc; close all;

amp = double(imread('cameraman.tif'));
phase = double(imread('trees.tif'));
s = size(amp);
phase = phase(1:s(1),1:s(2));
mask = zeros(s); mask(65:end-64,65:end-64)=1;
amp = (amp-min(amp(:)))/(max(amp(:))-min(amp(:))).*mask;
phase = (phase-min(phase(:)))/(max(phase(:))-min(phase(:))).*mask;

I = amp.*exp(1i.*phase);
I_dp = abs(fftshift(fft2(I))).^2;
[phase_rec, err] = GS_reconstruction(amp, I_dp, 1e-24);

figure
imshow(phase_rec, []), colorbar, title('sample phase');

figure
plot(log10(err)), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error - log scale');

%% FienupHIO reconstruction
clear; clc; close all;
cd('C:\Users\arielv\Desktop\תואר שני\חניכה\phase retrival algorithms');

amp = double(imread('trees.tif'));
phase = double(imread('cameraman.tif'));
s = size(phase);
div = 0.3; % how much of image there is in every axis
div_par = -0.5*div+0.5;
amp = amp(1:s(1),1:s(2));
mask = ones(s).*zeros(s); mask(s(1)*div_par+1:end-s(1)*div_par,s(2)*div_par+1:end-s(2)*div_par)=1;
amp = (amp-min(amp(:)))/(max(amp(:))-min(amp(:))).*mask;
if max(abs(phase(:)))>min(abs(phase(:)))
    phase = (phase-min(phase(:)))/(max(phase(:))-min(phase(:))).*mask;
end

I = amp.*exp(1i.*phase);
I_dp = abs(fftshift(fft2(I))).^2;
s = size(amp);
[phase_rec, err] = FienupHIO_reconstruction(I_dp,0.7, div, 1e-12);

figure
imshow(phase_rec, []), colorbar, title('sample phase');

figure
plot(log10(err)), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error - log scale');

%% ePIE - Round probe
clear; clc; close all;
addpath(genpath(pwd)); % add all subfolders of current folder to path

%%%%%%%%%%%%%%%%%% BEGIN EDITING %%%%%%%%%%%%%%%%%%%

% creating the data
amp = double(imread('liftingbody.png')); % double(imread('trees.tif'));
phase = double(imread('westconcordorthophoto.png')); % double(imread('cameraman.tif'));

[amp,phase] = cut_min_ind(amp,phase);
%mask = ones(s).*zeros(s); mask(s(1)*div_par+1:end-s(1)*div_par,s(2)*div_par+1:end-s(2)*div_par)=1;
amp = (amp-min(amp(:)))/(max(amp(:))-min(amp(:)));
if max(abs(phase(:)))>min(abs(phase(:)))
    phase = (phase-min(phase(:)))/(max(phase(:))-min(phase(:)));
end
s = size(amp);
% parameters of the probe
r = 25;
n_row = 9; n_col = 9;
percent = 70;
err = percent/7;
c = s/2;
iterations = 50;

%%%%%%%%%%%%%%%%%% END EDITING %%%%%%%%%%%%%%%%%%%

% params.radius = r;
params.radius = r.*(1+rand(1)/3); % put an error in radius
params.im_size = s;
params.center = c;

params.ProbeType = 'RoundProbe';
params.LoadProbe = 'RoundProbe(params.radius,params.im_size,params.center);';
params.iterations = iterations;

probe = RoundProbe(r,s,c);
[probes,cen] = probe.make_probes([n_row,n_col],percent,err);
% cen_plot = reshape(cen,[size(cen,1)*size(cen,2),size(cen,3)]);
% figure();
% scatter(cen_plot(:,1),cen_plot(:,2),'+')
% p=zeros(s);
% for i=1:n_row
%     for j=1:n_col
%         p = p+probes(:,:,i,j);
%     end
% end
% imshow(p>0,[]);impixelinfo;

obj = amp.*exp(i*phase);
probes = repmat(obj,[1,1,size(probes, 3), size(probes, 4)]).*probes;
for i=1:size(probes, 3)
    for j=1:size(probes,4)
        %I_dp(:,:,i,j) = abs(fftshift(fft2(probes(:,:,i,j)))/s(1)).^2;
        I_dp(:,:,i,j) = abs(fftshift(fft2(probes(:,:,i,j)))).^2;
    end
end
[object_rec, probe_rec, err] = ePIE_reconstruction(I_dp, cen, obj, 0.0001,params);

figure('position',get(0,'ScreenSize')); 
subplot(2,2,1);
imshow(abs(object_rec), []), colorbar, title('Object Magnitude','FontSize',16); %show current phase-image
subplot(2,2,2);
imshow(angle(object_rec), []), colorbar, title('Object Phase','FontSize',16); %show current phase-image
subplot(2,2,3);
imshow(abs(probe_rec), []), colorbar, title('Probe Magnitude','FontSize',16); %show current phase-image
subplot(2,2,4);
imshow(angle(probe_rec), []), colorbar, title('Probe Phase','FontSize',16); %show current phase-image
impixelinfo;

figure();
plot(err), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error');
%% ePIE - Gaussian probe
clear; clc; close all;
addpath(genpath(pwd)); % add all subfolders of current folder to path

%%%%%%%%%%%%%%%%%% BEGIN EDITING %%%%%%%%%%%%%%%%%%%

% creating the data
amp = double(imread('trees.tif')); % double(imread('liftingbody.png'));
phase = double(imread('cameraman.tif')); % double(imread('westconcordorthophoto.png'));

[amp,phase] = cut_min_ind(amp,phase);
amp = (amp-min(amp(:)))/(max(amp(:))-min(amp(:)));
if max(abs(phase(:)))>min(abs(phase(:)))
    phase = (phase-min(phase(:)))/(max(phase(:))-min(phase(:)));
end
s = size(amp);

amp = imresize(amp,819/s(1));phase = imresize(phase,819/s(1));
s = size(amp);

% parameters of the probe
sig = 20*3.5; sig = sig.*3.45e-6; % sig pixels
pixel_pitch = [3.45e-6 3.45e-6];
lambda = 850e-9;
n_row = 5; n_col = 5;
percent = 70;
err = 0;%percent/10;
c = s/2;
iterations = 50;

%%%%%%%%%%%%%%%%%% END EDITING %%%%%%%%%%%%%%%%%%%

% params.sig = sig;
params.sig = sig;%.*(1+rand(1)/3); % put an error in sigma
params.pixel_pitch = pixel_pitch;
params.lambda = lambda;
params.im_size = s;
params.center = c;

params.ProbeType = 'GaussProbe';
params.LoadProbe = 'GaussProbe(params.sig,params.pixel_pitch,params.lambda,params.im_size,params.center);';
params.iterations = iterations;

probe = GaussProbe(sig,pixel_pitch,lambda,s,c);
[probes,cen] = probe.make_probes([n_row,n_col],percent,err);

obj = amp.*exp(i*phase);
probes = repmat(obj,[1,1,size(probes, 3), size(probes, 4)]).*probes;
for i=1:size(probes, 3)
    for j=1:size(probes,4)
        I_dp(:,:,i,j) = abs(fftshift(fft2(probes(:,:,i,j)))).^2;
    end
end
[object_rec, probe_rec, err] = ePIE_reconstruction(I_dp, cen, probe.probe, probe.center, obj, 0.0001,params);

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
scatter(reshape(cen(:,:,2),[1,n_row*n_col]),reshape(cen(:,:,1),[1,n_row*n_col]),'x'); axis([1 s(2) 1 s(1)]); %show centers position
sgtitle('Reconstruction of object and Probe, amplitude and phase:','FontSize',24);
impixelinfo;

% figure();
% plot(err), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error'); % show Error VS iterations

%% ePIE - LP probe
clear; clc; close all;
addpath(genpath(pwd)); % add all subfolders of current folder to path

%%%%%%%%%%%%%%%%%% BEGIN EDITING %%%%%%%%%%%%%%%%%%%

% creating the data
amp = double(imread('trees.tif')); % double(imread('liftingbody.png'));
phase = double(imread('cameraman.tif')); % double(imread('westconcordorthophoto.png'));
[amp,phase] = cut_min_ind(amp,phase);
%mask = ones(s).*zeros(s); mask(s(1)*div_par+1:end-s(1)*div_par,s(2)*div_par+1:end-s(2)*div_par)=1;
amp = (amp-min(amp(:)))/(max(amp(:))-min(amp(:)));
if max(abs(phase(:)))>min(abs(phase(:)))
    phase = (phase-min(phase(:)))/(max(phase(:))-min(phase(:)));
end
s = size(amp);
% parameters of the probe
order.l = 1; order.m = 1; order.polarization = 'x';
fiber_params.CoreRad = 100e-6; % 100 microns core radius 
fiber_params.CladRad = 250e-6; % 250 microns clad radius
fiber_params.lambda = 850e-9; % 850 nm central wavelength
fiber_params.n1 = refrIndex('SILICA', fiber_params.lambda); % Silica is core material
fiber_params.n2 = refrIndex('SILICACLAD', fiber_params.lambda);
im_size = s;
center = s/2;
pixel_pitch = [3.45e-6 3.45e-6];

n_row = 8; n_col = 8;
percent = 80;
err = percent/10;
iterations = 100;

%%%%%%%%%%%%%%%%%% END EDITING %%%%%%%%%%%%%%%%%%%

params.order = order;
params.fiber_params = fiber_params;
params.fiber_params.CoreRad = params.fiber_params.CoreRad.*(1+rand(1)/5); %ambiguity in probe size
params.im_size = im_size;
params.center = center;
params.pixel_pitch = pixel_pitch;

params.ProbeType = 'LPmodeProbe';
params.LoadProbe = 'LPmodeProbe(params.order,params.fiber_params,params.pixel_pitch,params.im_size,params.center);';
params.iterations = iterations;

probe = LPmodeProbe(order,fiber_params,pixel_pitch,im_size,center);
[probes,cen] = probe.make_probes([n_row,n_col],percent,err);

obj = amp.*exp(i*phase);
probes = repmat(obj,[1,1,size(probes, 3), size(probes, 4)]).*probes;
for i=1:size(probes, 3)
    for j=1:size(probes,4)
        I_dp(:,:,i,j) = abs(fftshift(fft2(probes(:,:,i,j)))).^2;
    end
end
[object_rec, probe_rec, err] = ePIE_reconstruction(I_dp, cen, probe.probe, probe.center, obj, 0.0001,params);

figure('position',[100 100 1200 600]); 
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
scatter(reshape(cen(:,:,2),[1,n_row*n_col]),reshape(cen(:,:,1),[1,n_row*n_col]),'x'); axis([1 s(2) 1 s(1)]); %show centers position
sgtitle('Reconstruction of object and Probe, amplitude and phase:','FontSize',24);
impixelinfo;

% figure();
% plot(err), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error'); % show Error VS iterations

%% TIMP reconstruction - 2 objects

clear; clc; close all;
addpath(genpath(pwd)); % add all subfolders of current folder to path

%%%%%%%%%%%%%%%%%% BEGIN EDITING %%%%%%%%%%%%%%%%%%%
Create_numbers_image();

% creating the data
amp1 = imresize(a7,256/size(a7,1)); % double(imread('trees.tif')); % double(imread('liftingbody.png'));
phase1 = imresize(a8,256/size(a8,1)); %double(imread('cameraman.tif')); % double(imread('westconcordorthophoto.png'));
[amp1,phase1] = cut_min_ind(amp1,phase1);
amp1 = normalize_image(amp1); phase1 = normalize_image(phase1);
obj1 = amp1.*exp(i*phase1);

amp2 = imresize(a4,256/size(a4,1)); % double(imread('liftingbody.png'));
phase2 = imresize(a5,256/size(a5,1)); %double(imread('westconcordorthophoto.png'));
[amp2,phase2] = cut_min_ind(amp2,phase2);
amp2 = normalize_image(amp2); phase2 = normalize_image(phase2);
obj2 = amp2.*exp(i*phase2);

[obj1,obj2] = cut_min_ind(obj1,obj2);
s = size(obj1);

% parameters of the probes
order1.l = 1; order1.m = 1; order1.polarization = 'x';
order2.l = 1; order2.m = 1; order2.polarization = 'y';
fiber_params.CoreRad = 100e-6; % 100 microns core radius 
fiber_params.CladRad = 250e-6; % 250 microns clad radius
fiber_params.lambda = 850e-9; % 850 nm central wavelength
fiber_params.n1 = refrIndex('SILICA', fiber_params.lambda); % Silica is core material
fiber_params.n2 = refrIndex('SILICACLAD', fiber_params.lambda);
im_size = s;
center = s/2;
pixel_pitch = [3.45e-6 3.45e-6];

n_row = 8; n_col = 8;
percent = 80;
err = percent/10;
iterations = 50;

%%%%%%%%%%%%%%%%%% END EDITING %%%%%%%%%%%%%%%%%%%

% params.order = order1;
% params.fiber_params = fiber_params;
% params.fiber_params.CoreRad = params.fiber_params.CoreRad.*(1+rand(1)/5); %ambiguity in probe size
% params.im_size = im_size;
% params.center = center;
% params.pixel_pitch = pixel_pitch;
% 
% params.ProbeType = 'LPmodeProbe';
% params.LoadProbe = 'LPmodeProbe(params.order,params.fiber_params,params.pixel_pitch,params.im_size,params.center);';
params.iterations = iterations;

% make probes
probe1 = LPmodeProbe(order1,fiber_params,pixel_pitch,im_size,center);
[probes1,cen] = probe1.make_probes([n_row,n_col],percent,err);
probes1 = repmat(obj1,[1,1,size(probes1, 3), size(probes1, 4)]).*probes1;

probe2 = LPmodeProbe(order2,fiber_params,pixel_pitch,im_size,center);
[probes2,~] = probe2.make_probes([n_row,n_col],percent,err,cen);
probes2 = repmat(obj2,[1,1,size(probes2, 3), size(probes2, 4)]).*probes2;

for i=1:size(probes1, 3)
    for j=1:size(probes1,4)
        I_dp(:,:,i,j) = abs(fftshift(fft2(probes1(:,:,i,j)))).^2 + abs(fftshift(fft2(probes2(:,:,i,j)))).^2;
    end
end
probe_guess(:,:,1) = probe1.probe; probe_guess(:,:,2) = probe2.probe;
obj_comp(:,:,1) = obj1; obj_comp(:,:,2) = obj2;
[object_rec, probe_rec, err] = ePIE_reconstruction_TIMP(I_dp, cen, probe_guess, probe1.center, obj_comp, 0.0001,params);

for ii=1:size(probe_guess,3)
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
    scatter(reshape(cen(:,:,2),[1,n_row*n_col]),reshape(cen(:,:,1),[1,n_row*n_col]),'x'); axis([1 s(2) 1 s(1)]); %show centers position
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
