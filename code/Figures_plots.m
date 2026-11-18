%% Figures
scriptFolder = fileparts(mfilename('fullpath'));

N = 2^13;
Nim = 440;
lambda = 1030e-9; %[m] central wavelenght
n = [5, 5]; % number of spots at the pinholes plane
f1 = 0.1; f2 = 0.06; fms = 0.15; %[m] focal length of the lenses 
pixel_pitch = 3.45e-6; % [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% definition of real grid:   
dx = pixel_pitch; dy = pixel_pitch;
x_vec = dx*(-N/2:N/2-1); y_vec = dy*(-N/2:N/2-1);
[x,y] = meshgrid(x_vec, y_vec);

% definition of reciprocoral grid:
dkx = 1/(N*dx); dky = 1/(N*dy); % [Hz] smallest element of frequency
kx_vec = (-N/2:N/2-1); ky_vec = (-N/2:N/2-1);
[kx,ky] = meshgrid(kx_vec, ky_vec);
kx = fftshift(kx)/(N*dx); ky = fftshift(ky)/(N*dx);

% definition of the reciprocoral and real grid after scaling by the f2 lens:
dk = pixel_pitch/(lambda*f2);
kx_vec_sim = dk*(-N/2:N/2-1); ky_vec_sim = dk*(-N/2:N/2-1);
[kx_sim,ky_sim] = meshgrid(kx_vec_sim, ky_vec_sim);

% real grid after scaling by the f2 lens, with block number of pixels:
dx_im = 1/(Nim*dk);
x_vec_im = dx_im*(-Nim/2:Nim/2-1); y_vec_im = dx_im*(-Nim/2:Nim/2-1);
[x_im,y_im] = meshgrid(x_vec_im, y_vec_im);

[xx,yy] = meshgrid(-Nim/2:Nim/2-1,-Nim/2:Nim/2-1);

Factor = dx^2*Nim/(lambda*f2); % dx_cam = dx_obj * Factor
ind = 180:259;
load('Figures\sim_ms-5x5_fms-150_f1-100_f2-60.mat');

%% load recorded intensities from object plane and HDR reconstructions
probes_file = 'Figures\reconstructions_010922';
% reconstructions of probes from HDR data:
probe_rec_1 = load(strcat(probes_file,"\probe_no1_reconstruction.mat")); probe_rec_1 = probe_rec_1.probe_rec_1;
probe_rec_2 = load(strcat(probes_file,"\probe_no2_reconstruction.mat")); probe_rec_2 = probe_rec_2.probe_rec_2;
probe_rec_3 = load(strcat(probes_file,"\probe_no3_reconstruction.mat")); probe_rec_3 = probe_rec_3.probe_rec_3;
% recordings of probes - average of 10 images:
probe_recorded_1 = load(strcat(probes_file,"\probe_no1_recorded.mat")); probe_recorded_1 = probe_recorded_1.pr5;
probe_recorded_2 = load(strcat(probes_file,"\probe_no2_recorded.mat")); probe_recorded_2 = probe_recorded_2.pr5;
probe_recorded_3 = load(strcat(probes_file,"\probe_no3_recorded.mat")); probe_recorded_3 = probe_recorded_3.pr5;
% reconstructions of objects from HDR data using recorded probes as initial guess:
object_rec_1 = load(strcat(probes_file,"\object_no1_reconstruction.mat")); object_rec_1 = object_rec_1.obj_rec_1;
object_rec_2 = load(strcat(probes_file,"\object_no2_reconstruction.mat")); object_rec_2 = object_rec_2.obj_rec_2;
object_rec_3 = load(strcat(probes_file,"\object_no3_reconstruction.mat")); object_rec_3 = object_rec_3.obj_rec_3;

f = figure; imagesc(x_vec_im(ind), y_vec_im(ind), probe_recorded_1(ind,ind)); title("Probe 1 recorded"); axis square tight; f.OuterPosition = [601.0000  334.6000  454.4000  506.4000];
f = figure; imagesc(x_vec_im(ind), y_vec_im(ind), probe_recorded_2(ind,ind)); title("Probe 2 recorded"); axis square tight; f.OuterPosition = [601.0000  334.6000  454.4000  506.4000];
f = figure; imagesc(x_vec_im(ind), y_vec_im(ind), probe_recorded_3(ind,ind)); title("Probe 3 recorded"); axis square tight; f.OuterPosition = [601.0000  334.6000  454.4000  506.4000];

%% load probes reconstructed from 1+2+3 single shot with known object
file = "Figures\031122_reconstructions_3_fibers_thresh_0.15_dist_0.014_object_7_count_33_probes1+2+3_exp1.mat";
probes = load(file);
probes = gather(probes.object_rec_mat3);

clear params
params.ind = 180:259;
params.layout = [1 3];
params.Padding = "Compact";
params.plotType = "complex";
params.InnerPosition = [365.8,405.8,871.2,334.4];
params.TileSpacing = "compact";
params.XData = x_vec_im*10^3; params.YData = y_vec_im*10^3;
params.axisType = "axis square tight";
params.xlab = "X [mm]"; params.ylab = "Y [mm]";
params.FontSize = 14;


%% object for #33  with linear phase correction
% plot object #7 (using 3rd probe, full HDR, correction linear phaes:)
[xx, yy] = meshgrid(-220:219,-220:219);
load("Figures\031122_reconstructions_object_rec_mat_mask5_ang0_1-3_fiber1_4-6_fiber2_7-9_fiber3_1-no-ldr_2-ldr-p5_3-ldr-p25_031122.mat");
obj_plot = object_rec_mat1(:,:,7).*exp(-1i*18*xx(:,:)*2*pi/440-1i*2*pi/3);
obj_plot = abs(obj_plot).^2.*exp(1i*angle(obj_plot));
scale = 20;
figure; 
h = phplot(imresize(obj_plot(ind,ind),scale,"nearest"),1,1,"BottomLeft"); 
image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj_plot,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj_plot,1)*scale)*10^3,AddScaleBar(h.CData,dx_im/scale,0.4)); 
axis off square tight;

%%
% draw reconstructed probes using the #7 object above as known object,
% linear phase corrected:
params.layout = [1 3];
reconst_probes = cat(3,probes(:,:,1).*exp(1i*2*pi*(10*xx(:,:)-7*yy(:,:))/440),probes(:,:,2).*exp(1i*2*pi*(10*xx(:,:)-7*yy(:,:)+330)/440),probes(:,:,3).*exp(1i*2*pi*(12*xx(:,:)-12*yy(:,:)+180)/440));
draw_tiles(reconst_probes,params)


%% For IEEE Photonic Journal Article:
ind = 190:249;
probe_reconstructed_1 = probes(ind,ind,1).*exp(1i*2*pi*(10*xx(ind,ind)-7*yy(ind,ind))/440);
probe_reconstructed_2 = probes(ind,ind,2).*exp(1i*2*pi*(10*xx(ind,ind)-7*yy(ind,ind)+330)/440);
probe_reconstructed_3 = probes(ind,ind,3).*exp(1i*2*pi*(12*xx(ind,ind)-12*yy(ind,ind)+180)/440);
probe_reconstructed_1 = abs(probe_reconstructed_1).^2.*exp(1i*angle(probe_reconstructed_1)); % intensity with phase
probe_reconstructed_2 = abs(probe_reconstructed_2).^2.*exp(1i*angle(probe_reconstructed_2)); % intensity with phase
probe_reconstructed_3 = abs(probe_reconstructed_3).^2.*exp(1i*angle(probe_reconstructed_3)); % intensity with phase
XData = x_vec_im*10^3; YData = y_vec_im*10^3; sfont = 30;

% % 1a. reconstructions of probes, option 1 - with axes
% figure; tiledlayout(1,3,'TileSpacing',"tight",'Padding',"compact");nexttile;
% ftmp = phplot(imresize(probe_reconstructed_1,scale,"nearest"),1,1,"BottomLeft"); 
% image(linspace(XData(ind(1)),XData(ind(end)),size(probe_reconstructed_1,2)*scale),linspace(YData(ind(1)),YData(ind(end)),size(probe_reconstructed_1,1)*scale),ftmp.CData); axis square tight off; ylabel("Y [mm]"); set(gca,'xtick',[]); set(gca,'ytick',[-1 0 1]); set(gca,"FontSize",sfont);nexttile; % xlabel("X [mm]"); 
% ftmp = phplot(probe_reconstructed_3);
% image(XData(ind),YData(ind),ftmp.CData); axis square tight off;  set(gca,"FontSize",sfont);nexttile; %xlabel("X [mm]"); ylabel("Y [mm]");
% ftmp = phplot(probe_reconstructed_2);
% image(XData(ind),YData(ind),ftmp.CData); axis square tight off;  set(gca,"FontSize",sfont); % xlabel("X [mm]"); ylabel("Y [mm]");

% 1a. reconstructions of probes, option 1 - with axes
figure; tiledlayout(1,3,'TileSpacing',"tight",'Padding',"compact");nexttile;
ftmp = phplot(imresize(probe_reconstructed_1,scale,"nearest"),1,1,"BottomLeft"); image(ftmp.CData); axis square tight off;nexttile;
ftmp = phplot(probe_reconstructed_3); image(ftmp.CData); axis square tight off; nexttile;
ftmp = phplot(probe_reconstructed_2); image(ftmp.CData); axis square tight off;

% 1b. reconstructions of probes, option 2 - without axes
figure; tiledlayout(1,3,'TileSpacing',"tight",'Padding',"compact");nexttile;
ftmp = phplot(imresize(probe_reconstructed_1,scale,"nearest"),1,1,"BottomLeft"); image(AddScaleBar(ftmp.CData,dx_im/scale,0.4)); axis square tight off;nexttile;
ftmp = phplot(probe_reconstructed_3); image(ftmp.CData); axis square tight off; nexttile;
ftmp = phplot(probe_reconstructed_2); image(ftmp.CData); axis square tight off;

% draw recorded probes intensities from object plane
probes_file = 'C:\Users\arielv\Desktop\Master\images\Fibers\probe_measurements_130723';
% recordings of probes - average of 10 images:
probe1 = load_mean_image(strcat(probes_file,"\probe1_0p5ms")); imdark1 = load(strcat(probes_file,"\imdark_0p5ms.mat")); imdark1 = double(imdark1.frame); probe1 = noise_removal(probe1,imdark1);
probe2 = load_mean_image(strcat(probes_file,"\probe2_0p9ms")); imdark2 = load(strcat(probes_file,"\imdark_0p9ms.mat")); imdark2 = double(imdark2.frame); probe2 = noise_removal(probe2,imdark1);
probe3 = load_mean_image(strcat(probes_file,"\probe3_0p8ms")); imdark3 = load(strcat(probes_file,"\imdark_0p8ms.mat")); imdark3 = double(imdark3.frame); probe3 = noise_removal(probe3,imdark3);
% find center of the probes and move them to the center:
[cen_probe,~,~,~] = FindRefCenters(probe1,[1 1],10,30);
probe1 = abs(ShiftProbefft(probe1,size(probe1)/2-squeeze(cen_probe))); probe2 = abs(ShiftProbefft(probe2,size(probe2)/2-squeeze(cen_probe))); probe3 = abs(ShiftProbefft(probe3,size(probe3)/2-squeeze(cen_probe)));
% create vector of camera and find indexes for drawing: 
x_cam_vec = (-1500:1499)*dx; [~,ind1] = min(abs(x_cam_vec-x_vec_im(ind(1)))); [~,ind2] = min(abs(x_cam_vec-x_vec_im(ind(end)))); indx = ind1:ind2;

% draw recorded probes (without axes):
figure; tiledlayout(1,3,'TileSpacing',"tight",'Padding',"compact"); nexttile;
imagesc(AddScaleBar(probe1(indx,indx).^1,dx,0.4)); axis square tight off; colormap(magma); nexttile;
imagesc(probe2(indx,indx).^1); axis square tight off; colormap(magma); nexttile;
imagesc(probe3(indx,indx).^1); axis square tight off; colormap(magma);

%% For Optics Letters Article:
% 2a. recorded intensities of probes, option 1 - with axes
figure; tiledlayout(1,3,'TileSpacing',"tight",'Padding',"compact"); nexttile;
imagesc(XData(ind),YData(ind), probe_recorded_2(ind,ind).^1); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma); nexttile;
imagesc(XData(ind),YData(ind), probe_recorded_3(ind,ind).^1); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma); nexttile; %
imagesc(XData(ind),YData(ind), probe_recorded_1(ind,ind).^1); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma); % 

% 2b. reconstructions of probes, option 2 - without axes
figure; tiledlayout(1,3,'TileSpacing',"tight",'Padding',"compact"); nexttile;
imagesc(AddScaleBar(imresize(probe_recorded_2(ind,ind).^1,scale,"nearest"),dx_im/scale,0.3)); axis square tight off; colormap(magma); nexttile;
imagesc(probe_recorded_3(ind,ind).^1); axis square tight off; colormap(magma); nexttile;
imagesc(probe_recorded_1(ind,ind).^1); axis square tight off; colormap(magma);

% draw recorded probes intensities from object plane
probes_file = 'C:\Users\arielv\Desktop\Master\images\Fibers\probe_measurements_130723';
% recordings of probes - average of 10 images:
probe1 = load_mean_image(strcat(probes_file,"\probe1_0p5ms")); imdark1 = load(strcat(probes_file,"\imdark_0p5ms.mat")); imdark1 = double(imdark1.frame); probe1 = noise_removal(probe1,imdark1);
probe2 = load_mean_image(strcat(probes_file,"\probe2_0p9ms")); imdark2 = load(strcat(probes_file,"\imdark_0p9ms.mat")); imdark2 = double(imdark2.frame); probe2 = noise_removal(probe2,imdark1);
probe3 = load_mean_image(strcat(probes_file,"\probe3_0p8ms")); imdark3 = load(strcat(probes_file,"\imdark_0p8ms.mat")); imdark3 = double(imdark3.frame); probe3 = noise_removal(probe3,imdark3);
% find center of the probes and move them to the center:
[cen_probe,~,~,~] = FindRefCenters(probe1,[1 1],10,30);
probe1 = abs(ShiftProbefft(probe1,size(probe1)/2-squeeze(cen_probe))); probe2 = abs(ShiftProbefft(probe2,size(probe2)/2-squeeze(cen_probe))); probe3 = abs(ShiftProbefft(probe3,size(probe3)/2-squeeze(cen_probe)));
% create vector of camera and find indexes for drawing: 
x_cam_vec = (-1500:1499)*dx; [~,ind1] = min(abs(x_cam_vec-x_vec_im(ind(1)))); [~,ind2] = min(abs(x_cam_vec-x_vec_im(ind(end)))); indx = ind1:ind2;
% draw recorded probes (option a - with axes):
figure; tiledlayout(1,3,'TileSpacing',"tight",'Padding',"compact"); nexttile;
imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe1(indx,indx).^1); axis square tight off; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,'ytick',[-1 0 1]); set(gca,"FontSize",sfont); colormap(magma); nexttile;
imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe2(indx,indx).^1); axis square tight off; xlabel("X [mm]"); set(gca,'ytick',[]); set(gca,"FontSize",sfont); colormap(magma); nexttile; % ylabel("Y [mm]");
imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe3(indx,indx).^1); axis square tight off; xlabel("X [mm]"); set(gca,'ytick',[]); set(gca,"FontSize",sfont); colormap(magma); % ylabel("Y [mm]");
% draw recorded probes (option b - without axes):
figure; tiledlayout(1,3,'TileSpacing',"tight",'Padding',"compact"); nexttile;
imagesc(AddScaleBar(probe1(indx,indx).^1,dx,0.3)); axis square tight off; colormap(magma); nexttile;
imagesc(probe2(indx,indx).^1); axis square tight off; colormap(magma); nexttile;
imagesc(probe3(indx,indx).^1); axis square tight off; colormap(magma);

% draw recorded probes (option c - without all):
figure; tiledlayout(1,3,'TileSpacing',"tight",'Padding',"compact"); nexttile;
imagesc(probe1(indx,indx).^1); axis square tight off; colormap(magma); nexttile;
imagesc(probe2(indx,indx).^1); axis square tight off; colormap(magma); nexttile;
imagesc(probe3(indx,indx).^1); axis square tight off; colormap(magma);

%%
params.FontSize = 14;
params.layout = [1 3];
params.plotType = "real";
params.colormapType = 'magma';
% params.ScaleBar = 0.5;
draw_tiles(cat(3,probe_recorded_2,probe_recorded_1,probe_recorded_3),params)

%% draw the recorded diffraction pattern:
imcen = load("C:\Users\arielv\Desktop\Master\images\Fibers\031122\L3M_probe_centers_0p8ms_031122.mat"); imcen = double(imcen.frame);
imdark_cen = load("C:\Users\arielv\Desktop\Master\images\Fibers\031122\image_dark_0p8ms_031122.mat"); imdark_cen = double(imdark_cen.frame);
imcen_n = noise_removal(imcen,imdark_cen); % remove noise
imrec = load("C:\Users\arielv\Desktop\Master\images\Fibers\031122\fiber_3M_no1+no2+no3_031122\fibers_1-2-3_measurements_mask5_target_sp_ldr_p5_0p035ms_031122.mat"); imrec = double(imrec.frame);
imdark = load("C:\Users\arielv\Desktop\Master\images\Fibers\031122\fiber_3M_no1+no2+no3_031122\fibers_1-2-3_image_dark_mask5_sp_ldr_p5_0p035ms_031122.mat"); imdark = double(imdark.frame);
imrec_n = noise_removal(imrec,imdark,'thresholding',0.15); % remove noise

imrec_n = remove_burned_pixels(imrec_n);

[sim_cam_comp,imcen_n]=CutMinInd(abs(sim_ms),imcen_n);
[cen_rec,~,ang_tilt,bt] = FindRefCenters(imcen_n,n,10,30);
subimage_rec = [1 1].*440; % take only square blocks, for comfort
[probes_ref,~] = cut_probes(imcen_n,cen_rec,subimage_rec); % cut the image without object for registration and energy normalization
[cen_sim,~,~,~] = FindRefCenters(sim_cam_comp,n,10,30);
[probes_sim,~] = cut_probes(sim_cam_comp,cen_sim,subimage_rec); % cut the image without object for registration and energy normalization
[probes_real,~] = cut_probes(imrec_n,cen_rec(end:-1:1,end:-1:1,:),subimage_rec); % cut the probes of measured image according to the reference image. invert centers order for reconstruction - top left block is suitable for bottom right probe

% get parameters for measured image registration, normalize and register each block independetly:
output_ref = zeros(n(1),n(2),4);
% imdraw = zeros(n.*subimage_rec);
for i=1:n(1)
    for j=1:n(2)
        [output_ref(i,j,:), probe_tmp] = dftregistration(fft2(probes_sim(:,:,i,j)),fft2(probes_ref(:,:,i,j)),10); %produce block by block registration, due to small errors in Z-axis of the camera (so the centers are not exactly at Fourier plane)
        probes_ref_reg(:,:,i,j) = abs(ifft2(probe_tmp));
%         probes_real(:,:,i,j) = normalize_energy(probes_real(:,:,i,j),probes_ref(:,:,i,j)); % normalize energy
        probes_real_reg(:,:,i,j) = abs(ShiftProbefft(probes_real(:,:,i,j),[output_ref(i,j,3) output_ref(i,j,4)])); %  block-block registration
%         imdraw((i-1)*subimage_rec(1)+1:i*subimage_rec(1),(j-1)*subimage_rec(2)+1:j*subimage_rec(2)) = probes_real_reg(:,:,i,j);
    end
end

%%
figure; imagesc(AddScaleBar(imrec_n.^0.25,dx,0.3,"BottomLeft"),[0 6]); axis off tight square; colormap(magma)
figure; imagesc(AddScaleBar(probes_real(:,:,4,3).^0.25,dx,0.3,"BottomLeft"),[0 6]); axis off tight square; colormap(magma)
%%
clear params
params.layout = [5 5];
% params.ind = 191:250;
params.TileSpacing = "none";
params.Padding = "Compact";
params.plotType = "real";
params.InnerPosition = [74.6000   85.8000  703.2000  687.2000];
draw_tiles(probes_real_reg.^0.5,params)
%%
% % draw of the recorded diffraction patterns of 3 holes TIMP recostructions:
% load('C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\object_rec_d=10mm_thresh=0.05_filt=0.6_9x9.mat');
% ind = 190:250;
% scale = 20;
% sfont = 18;
% obj1 = abs(object_rec(ind,ind,1)).^(1).*exp(1i*angle(object_rec(ind,ind,1))).*exp(1i*2*pi*(-5*xx(ind,ind)-2*yy(ind,ind))/440).*exp(1i*2*pi*0.15);
% obj2 = abs(object_rec(ind,ind,2)).^(1).*exp(1i*angle(object_rec(ind,ind,2))).*exp(1i*2*pi*(-18*xx(ind,ind)+4*yy(ind,ind))/440).*exp(1i*2*pi*1);
% obj3 = abs(object_rec(ind,ind,3)).^(1).*exp(1i*angle(object_rec(ind,ind,3))).*exp(1i*2*pi*(-18*xx(ind,ind)+2*yy(ind,ind))/440).*exp(1i*2*pi*0.4);
% figure; h1 = phplot(imresize(obj1,scale,"nearest"),1,1,'BottomLeft');
% image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj1,1)*scale)*10^3,h1.CData);
% axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
% figure; h2 = phplot(imresize(obj2,scale,"nearest"));
% image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj2,1)*scale)*10^3,h2.CData);
% axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
% figure; h3 = phplot(imresize(obj3,scale,"nearest"));
% image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj3,1)*scale)*10^3,h3.CData);
% axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
% %%
% %comparison to ssp reconstructed image:
% load('C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\probe1-after_gain7_object_rec_d=10mm_thresh=0.1.mat')
% obj_comp = object_rec(ind,ind).^(2).*exp(1i*angle(object_rec(ind,ind,1))).*exp(1i*2*pi*(-7*xx(ind,ind)-4*yy(ind,ind))/440).*exp(1i*2*pi*0.6);
% figure; h4 = phplot(imresize(obj_comp,scale,"nearest"));
% image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj3,1)*scale)*10^3,h4.CData);
% axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
%% For Article "Single-Shot Ptychographic Imaging of Non-Repetitive Ultrafast Events"
% draw of the recorded diffraction patterns of 3 holes TIMP recostructions:
% load('C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\object_rec_d=10mm_thresh=[0_0_0_0_0_0.01_0.01_0.01_0.01_0.01_0.02_0.02_0.02_0.02_0.02]_filt=0.6_9x9_b.mat');
% % obj_num = [2 4 5 7 12 14]; %options
% obj_num = 2;
load('C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\new_reconstructions_no_saturations_051023\object_rec_d=10mm_thresh=[0_0_0_0_0_0.01_0.01_0.01_0.01_0.01_0.02_0.02_0.02_0.02_0.02]_filt=0.9_9x9_d_using_previous_reconst.mat');
% obj_num = [4 6 10 14 15]; %options
obj_num = 4;
object_rec = object_rec_mat(:,:,:,obj_num);
ind = 190:250;
scale = 20;
sfont = 14;
obj1 = abs(object_rec(ind,ind,1)).^(2).*exp(1i*angle(object_rec(ind,ind,1))).*exp(1i*2*pi*(-8*xx(ind,ind)-4*yy(ind,ind))/440).*exp(1i*2*pi*0.18);
obj2 = abs(object_rec(ind,ind,2)).^(2).*exp(1i*angle(object_rec(ind,ind,2))).*exp(1i*2*pi*(-20*xx(ind,ind)+4*yy(ind,ind))/440).*exp(1i*2*pi*0.75);
obj3 = abs(object_rec(ind,ind,3)).^(2).*exp(1i*angle(object_rec(ind,ind,3))).*exp(1i*2*pi*(-20*xx(ind,ind)+2*yy(ind,ind))/440).*exp(1i*2*pi*0.6);
figure; 
subplot(1,3,1); h1 = phplot(imresize(obj1,scale,"nearest"),1,1,'BottomLeft');
image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj1,1)*scale)*10^3,h1.CData);
axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
subplot(1,3,2); h2 = phplot(imresize(obj2,scale,"nearest"));
image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj2,1)*scale)*10^3,h2.CData);
axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
subplot(1,3,3); h3 = phplot(imresize(obj3,scale,"nearest"));
image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj3,1)*scale)*10^3,h3.CData);
axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
%%
%comparison to ssp reconstructed image:
load('C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\probe3-afer_gain7_object_rec_d=10mm_thresh=0_filt=[0.5_1_2_3_4_5]_9x9.mat');
obj_num = 6;
object_rec = object_rec_mat(:,:,obj_num);
obj_comp = abs(object_rec(ind,ind)).^(2).*exp(1i*angle(object_rec(ind,ind,1))).*exp(1i*2*pi*(-15*xx(ind,ind)+5*yy(ind,ind))/440).*exp(1i*2*pi*0.7);
figure; h4 = phplot(imresize(obj_comp,scale,"nearest"));
image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj3,1)*scale)*10^3,h4.CData);
axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
%% probe images
% draw of the recorded diffraction patterns of 3 probes TIMP recostructions:
% load('C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\probe_rec_d=10mm_thresh=[0_0_0_0_0_0.01_0.01_0.01_0.01_0.01_0.02_0.02_0.02_0.02_0.02]_filt=0.6_9x9_b.mat');
load('C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\new_reconstructions_no_saturations_051023\probe_rec_d=10mm_thresh=[0_0_0_0_0_0.01_0.01_0.01_0.01_0.01_0.02_0.02_0.02_0.02_0.02]_filt=0.9_9x9_d_using_previous_reconst.mat');
figure;
for ii=1:numel(obj_num)
probe_rec = probe_rec_mat(:,:,:,4);
ind = 190:250;
scale = 20;
sfont = 18;
obj1 = abs(probe_rec(ind,ind,1)).^(2).*exp(1i*angle(probe_rec(ind,ind,1))).*exp(1i*2*pi*(-0*xx(ind,ind)-0*yy(ind,ind))/440).*exp(1i*2*pi*0.0);
obj2 = abs(probe_rec(ind,ind,2)).^(2).*exp(1i*angle(probe_rec(ind,ind,2))).*exp(1i*2*pi*(+13*xx(ind,ind)+0*yy(ind,ind))/440).*exp(1i*2*pi*0.65);
obj3 = abs(probe_rec(ind,ind,3)).^(2).*exp(1i*angle(probe_rec(ind,ind,3))).*exp(1i*2*pi*(+3*xx(ind,ind)+2*yy(ind,ind))/440).*exp(1i*2*pi*0.6);
subplot(1,3,1); h1 = phplot(imresize(obj1,scale,"nearest"),1,1,'BottomLeft');
image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj1,1)*scale)*10^3,h1.CData);
axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
subplot(1,3,2); h2 = phplot(imresize(obj2,scale,"nearest"));
image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj2,1)*scale)*10^3,h2.CData);
axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
subplot(1,3,3); h3 = phplot(imresize(obj3,scale,"nearest"));
image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(obj1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj3,1)*scale)*10^3,h3.CData);
axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);

% option b:
figure;
subplot(1,3,1); h1 = phplot(imresize(obj1,scale,"nearest"),1,1,'BottomLeft');
image(AddScaleBar(h1.CData,dx_im/scale,0.5)); axis square tight off;
subplot(1,3,2); h2 = phplot(imresize(obj2,scale,"nearest"));
image(h2.CData); axis square tight off;
subplot(1,3,3); h3 = phplot(imresize(obj3,scale,"nearest"));
image(h3.CData); axis square tight off;
end
%%
% draw recorded probes intensities from object plane
probes_file = 'C:\Users\arielv\Desktop\Master\images\Fibers\probe_measurements_130723';
% recordings of probes - average of 10 images:
probe1 = load_mean_image(strcat(probes_file,"\probe1_0p5ms")); imdark1 = load(strcat(probes_file,"\imdark_0p5ms.mat")); imdark1 = double(imdark1.frame); probe1 = noise_removal(probe1,imdark1);
probe2 = load_mean_image(strcat(probes_file,"\probe2_0p9ms")); imdark2 = load(strcat(probes_file,"\imdark_0p9ms.mat")); imdark2 = double(imdark2.frame); probe2 = noise_removal(probe2,imdark1);
probe3 = load_mean_image(strcat(probes_file,"\probe3_0p8ms")); imdark3 = load(strcat(probes_file,"\imdark_0p8ms.mat")); imdark3 = double(imdark3.frame); probe3 = noise_removal(probe3,imdark3);
% find center of the probes and move them to the center:
[cen_probe,~,~,~] = FindRefCenters(probe1,[1 1],10,30);
probe1 = abs(ShiftProbefft(probe1,size(probe1)/2-squeeze(cen_probe))); probe2 = abs(ShiftProbefft(probe2,size(probe2)/2-squeeze(cen_probe))); probe3 = abs(ShiftProbefft(probe3,size(probe3)/2-squeeze(cen_probe)));
% create vector of camera and find indexes for drawing: 
x_cam_vec = (-1500:1499)*dx; [~,ind1] = min(abs(x_cam_vec-x_vec_im(ind(1)))); [~,ind2] = min(abs(x_cam_vec-x_vec_im(ind(end)))); indx = ind1:ind2;
% draw recorded probes:
figure;
subplot(1,3,1); imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe1(indx,indx).^1); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma);
subplot(1,3,2); imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe2(indx,indx).^1); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma);
subplot(1,3,3); imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe3(indx,indx).^1); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma);
% draw recorded probes (option b):
figure;
subplot(1,3,1); imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, AddScaleBar(probe1(indx,indx).^1,dx,0.5)); axis square tight off; colormap(magma);
subplot(1,3,2); imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe2(indx,indx).^1); axis square tight off; colormap(magma);
subplot(1,3,3); imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe3(indx,indx).^1); axis square tight off; colormap(magma);
%% option b in one figure:
f= figure; tiledlayout(2,3,'TileSpacing','tight','Padding','tight'); nexttile; 
h1 = phplot(imresize(obj1,scale,"nearest"),1,1,'BottomLeft');
image(AddScaleBar(h1.CData,dx_im/scale,0.5)); axis square tight off; nexttile;
h2 = phplot(imresize(obj2,scale,"nearest"));
image(h2.CData); axis square tight off; nexttile;
h3 = phplot(imresize(obj3,scale,"nearest"));
image(h3.CData); axis square tight off; nexttile;
imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, AddScaleBar(probe1(indx,indx).^1,dx,0.5)); axis square tight off; colormap(magma); nexttile;
imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe2(indx,indx).^1); axis square tight off; colormap(magma); nexttile;
imagesc(x_cam_vec(indx)*10^3,x_cam_vec(indx)*10^3, probe3(indx,indx).^1); axis square tight off; colormap(magma);

%% 
% draw the image on the camera before reconstruction
imrec1 = load("C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\probe1+2+3_during_gain24.mat"); imrec1 = double(imrec1.frame);
imrec2 = load("C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\probe1+2+3_during_gain24.mat"); imrec2 = imgaussfilt(double(imrec2.frame),0.6,"FilterSize",[9 9]);
imdark24 =load("C:\Users\arielv\Desktop\Master\images\carbon_target\120623\1+2+3\att3\dark24.mat"); imdark24 = imgaussfilt(double(imdark24.frame),0.6,"FilterSize",[9 9]);
imrec_n1 = noise_removal(imrec2,imdark24,'thresholding',0.0); % remove noise
imrec_n2 = noise_removal(imrec2,imdark24,'thresholding',0.18); % remove noise

power = 0.25;
scale = 0.5;
%option a - with axis scale
figure;
subplot(1,4,1); imagesc(x_vec(end/2-size(imrec1,2)/2:end/2+size(imrec1,2)/2-1)*10^3,y_vec(end/2-size(imrec1,1)/2:end/2+size(imrec1,1)/2-1)*10^3,imrec1.^power); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma); ax = gca; ax.XAxisLocation = 'top'; a = gca; clims = a.CLim;
subplot(1,4,2); imagesc(x_vec(end/2-size(imrec2,2)/2:end/2+size(imrec2,2)/2-1)*10^3,y_vec(end/2-size(imrec2,1)/2:end/2+size(imrec2,1)/2-1)*10^3,imrec2.^power); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma); ax = gca; ax.XAxisLocation = 'top';
subplot(1,4,3); imagesc(x_vec(end/2-size(imrec_n1,2)/2:end/2+size(imrec_n1,2)/2-1)*10^3,y_vec(end/2-size(imrec_n1,1)/2:end/2+size(imrec_n1,1)/2-1)*10^3,imrec_n1.^power); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma); ax = gca; ax.XAxisLocation = 'top';
subplot(1,4,4); imagesc(x_vec(end/2-size(imrec_n2,2)/2:end/2+size(imrec_n2,2)/2-1)*10^3,y_vec(end/2-size(imrec_n2,1)/2:end/2+size(imrec_n2,1)/2-1)*10^3,imrec_n2.^power); axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont); colormap(magma); ax = gca; ax.XAxisLocation = 'top';
%option b - with scale bar
figure;
subplot(1,4,1); imagesc(x_vec(end/2-size(imrec1,2)/2:end/2+size(imrec1,2)/2-1)*10^3,y_vec(end/2-size(imrec1,1)/2:end/2+size(imrec1,1)/2-1)*10^3,AddScaleBar(imresize(imrec1.^power,scale,"nearest"),dx/scale,0.5,"TopLeft"),clims); axis square tight off; colormap(magma);
subplot(1,4,2); imagesc(x_vec(end/2-size(imrec2,2)/2:end/2+size(imrec2,2)/2-1)*10^3,y_vec(end/2-size(imrec2,1)/2:end/2+size(imrec2,1)/2-1)*10^3,imrec2.^power); axis square tight off; colormap(magma); 
subplot(1,4,3); imagesc(x_vec(end/2-size(imrec_n1,2)/2:end/2+size(imrec_n1,2)/2-1)*10^3,y_vec(end/2-size(imrec_n1,1)/2:end/2+size(imrec_n1,1)/2-1)*10^3,imrec_n1.^power); axis square tight off; colormap(magma); 
subplot(1,4,4); imagesc(x_vec(end/2-size(imrec_n2,2)/2:end/2+size(imrec_n2,2)/2-1)*10^3,y_vec(end/2-size(imrec_n2,1)/2:end/2+size(imrec_n2,1)/2-1)*10^3,imrec_n2.^power); axis square tight off; colormap(magma); 
% option c - 2x2 with tiles:
f= figure; tiledlayout(1,4,'TileSpacing','tight','Padding','tight'); nexttile; 
imagesc(x_vec(end/2-size(imrec1,2)/2:end/2+size(imrec1,2)/2-1)*10^3,y_vec(end/2-size(imrec1,1)/2:end/2+size(imrec1,1)/2-1)*10^3,AddScaleBar(imresize(imrec1.^power,scale,"nearest"),dx/scale,0.5,"TopRight"),clims); axis square tight off; colormap(magma); nexttile; 
imagesc(x_vec(end/2-size(imrec2,2)/2:end/2+size(imrec2,2)/2-1)*10^3,y_vec(end/2-size(imrec2,1)/2:end/2+size(imrec2,1)/2-1)*10^3,imrec2.^power); axis square tight off; colormap(magma); nexttile; 
imagesc(x_vec(end/2-size(imrec_n1,2)/2:end/2+size(imrec_n1,2)/2-1)*10^3,y_vec(end/2-size(imrec_n1,1)/2:end/2+size(imrec_n1,1)/2-1)*10^3,imrec_n1.^power); axis square tight off; colormap(magma); nexttile; 
imagesc(x_vec(end/2-size(imrec_n2,2)/2:end/2+size(imrec_n2,2)/2-1)*10^3,y_vec(end/2-size(imrec_n2,1)/2:end/2+size(imrec_n2,1)/2-1)*10^3,imrec_n2.^power); axis square tight off; colormap(magma); 

indx = 1850:2150; indy = 1450:1750;
figure; imagesc(imrec1(indx,indy).^power); axis square tight off; colormap(magma);
figure; imagesc(imrec2(indx,indy).^power); axis square tight off; colormap(magma);
figure; imagesc(imrec_n1(indx,indy).^power); axis square tight off; colormap(magma);
figure; imagesc(imrec_n2(indx,indy).^power); axis square tight off; colormap(magma);

%option b - with scale bar
scale = 5;
figure; imagesc(AddScaleBar(imresize(imrec1(indx,indy).^power,scale,"nearest"),dx/scale,0.5,"BottomLeft"),clims); axis square tight off; colormap(magma);
% figure; imagesc(aaa(indx,indy).^power); axis square tight off; colormap(magma);
%%
% inline images for SSP system
theta = 0.97;
p = GaussProbe(0.0003, [1 1]*pixel_pitch, lambda, [1 1]*N, [1 1]*N/2+1);
impinhole = FreePropagation(LensPropagation(FreePropagation(MultispotPropagation(p.probe,theta,n,lambda,x,y),fms,lambda,kx,ky),fms,lambda,x,y),fms,lambda,kx,ky);
figure; imagesc(abs(CutMinInd(impinhole,ones(2400))));axis tight square off; colormap(jet);

rec = load("C:\Users\arielv\Desktop\Master\images\Fibers\280822\fiber_3M_no2_280822\measurements_mask5_target_hdr_280822.mat");
imrec_hdr = double(rec.frameHDR.frameArray); 
exposure = double(rec.frameHDR.exposure);
imdark_hdr = load("C:\Users\arielv\Desktop\Master\images\Fibers\280822\fiber_3M_no2_280822\image_dark_mask5_hdr_280822.mat"); 
imdark_hdr = double(imdark_hdr.frameArrayDark);
imrec_hdr_n = make_hdr(imrec_hdr,exposure,imdark_hdr); % remove noise and make hdr image
figure; imagesc(abs(CutMinInd(ShiftProbefft(imrotate(imrec_hdr_n,0.6),[30 -100]),ones(1500))).^0.5);axis tight square off; colormap(jet) % amplitude of field on camera

%%
% % draw probes
% probes_file = 'C:\Users\arielv\Desktop\Master\images\Fibers\010922\reconstructions';
% probe3 = load(strcat(probes_file,"\probe_no1_reconstruction.mat")); probe3 = probe3.probe_rec_1(ind,ind);
% probe1 = load(strcat(probes_file,"\probe_no2_reconstruction.mat")); probe1 = probe1.probe_rec_2(ind,ind);
% probe2 = load(strcat(probes_file,"\probe_no3_reconstruction.mat")); probe2 = probe2.probe_rec_3(ind,ind).*exp(1i*2*pi*(+15*xx(ind,ind)-1*yy(ind,ind))/440);
% 
% figure; h1 = phplot(imresize(probe1,scale,"nearest"));
% image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(probe1,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj1,1)*scale)*10^3,h1.CData);
% axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
% figure; h2 = phplot(imresize(probe2,scale,"nearest"));
% image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(probe2,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj2,1)*scale)*10^3,h2.CData);
% axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);
% figure; h3 = phplot(imresize(probe3,scale,"nearest"));
% image(linspace(x_vec_im(ind(1)),x_vec_im(ind(end)),size(probe3,2)*scale)*10^3,linspace(y_vec_im(ind(1)),y_vec_im(ind(end)),size(obj3,1)*scale)*10^3,h3.CData);
% axis square tight; xlabel("X [mm]"); ylabel("Y [mm]"); set(gca,"FontSize",sfont);

%% help functions
function draw_tiles(blocks,params)
    if ~exist("params","var") params = []; end
    fields = fieldnames(params);
    for i=1:numel(fields)
        eval([fields{i},' = params.',fields{i},';']);
    end
    if ~exist("layout","var") layout = size(blocks,[3 4]); end
    if ~exist("ind","var") ind = 1:size(blocks,1); end
    if ~exist("TileSpacing","var") TileSpacing = "tight"; end
    if ~exist("Padding","var") Padding = "tight"; end
    if ~exist("plotType","var") plotType = "real"; end
    if ~exist("axisType","var") axisType = "axis square tight off"; end

    blocks = reshape(blocks,[size(blocks,1),size(blocks,2),size(blocks,3)*size(blocks,4)]);

    f= figure; tiledlayout(layout(1),layout(2),'TileSpacing',TileSpacing,'Padding',Padding);
    for j=1:size(blocks,3)
        nexttile; 
        if plotType == "real" 
            if exist("XData","var")
                imagesc(XData,YData,blocks(ind,ind,j)); 
            else
                imagesc(blocks(ind,ind,j));
            end
        end
        if plotType == "complex"
            if exist("XData","var")
                ftmp = phplot(blocks(ind,ind,j));
                image(XData(ind),YData(ind),ftmp.CData);
            else
                phplot(blocks(ind,ind,j),1,1); 
            end
        end
        eval(axisType);
        if exist("titles","var") title(titles(j)); end
        if exist("FontSize","var") set(gca,"FontSize",FontSize); end
        if exist("xlab","var") xlabel(xlab); ylabel(ylab); end
        if exist("colormapType","var") colormap(colormapType); end
        if exist("InnerPosition","var") f.InnerPosition = InnerPosition; end
    end
end

function support = find_probe_support(probe)
    %support = (probe > max(probe(:))/exp(0.5)) + (probe < -max(probe(:))/exp(0.5)); % Sigma
    support = (probe > max(probe(:)).*0.05) + (probe < -max(probe(:)).*0.05);
    s = size(probe);
    [x,y] = meshgrid(1:s(2),1:s(1));
    dia1 = find(sum(support,1),1,"last")-find(sum(support,1),1); % diameter of active columns
    dia2 = find(sum(support,2),1,"last")-find(sum(support,2),1); % diameter of active rows
    support = (x-s(2)/2).^2./(dia1/2).^2+(y-s(1)/2).^2./(dia2/2).^2<=1;
end

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

function imout = load_mean_image(fullfile)
% function for load and averaging set of image in fullfile
% fullfile - path of *.mat files to load
% imout - averaged image [Nr X Nc double] 

    files = dir(fullfile);
    for i=3:length(files)
        file = load(strcat(fullfile,'\',files(i).name)); file = double(file.frame);
        if ~exist("imout","var")
            imout = zeros(size(file));
        end
        imout = imout+file;
    end
    imout = double(imout./(length(files)-2));
end

function im = remove_burned_pixels(im)
    im(1476:1479,1605:1609) = 0;
end