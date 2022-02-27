%% calculate mode patterns using different HG modes:

clear all; clc; close all;

% initials:
lambda = 500e-9; %[m] central wavelenght
W0 = 15e-6; %[m] waist 
z = 5e-3; %[m] propagation distance
n0 = 1;
dn = 0;
steps = 1000; %how many steps in a 
dz = z/steps;
zvec = dz:dz:z;

% definition of grid:
N = 256;
a = 10*W0; %width of image - about 0.75[mm] 
k0 = 2*pi*n0/lambda;
dx = a/N; dy = a/N;
[x,y]= meshgrid(-a/2:dx:a/2-dx,-a/2:dy:a/2-dy);
Wk0x = 2*pi/dx; Wk0y = 2*pi/dy; % [Hz] largest element of wavelenth for fourier transform
dkx = Wk0x/N; dky = Wk0y/N; % [Hz] smallest element of frequency
[kx,ky] = meshgrid(-Wk0x/2:dkx:Wk0x/2-dkx,-Wk0y/2:dky:Wk0y/2-dky);
r = sqrt(x.^2 + y.^2);
%% create different modes at z=0 positions
% n = 10; %how many modes to use
% psi0 = 1.*exp(-(r.^2)./W0.^2);
% iter = 0;
% % create first n^2 mpdes
% for ii=1:n
%     for jj=1:n
%         Psi(:,:,ii,jj) = psi0.*hermiteH(ii-1,sqrt(2).*x/W0).*hermiteH(jj-1,sqrt(2).*y/W0);
%         iter = iter+1;
%         fprintf('Iteration: %d\n', iter); 
%     end
% end
%% Read 100 HG patterns
Psi = load('C:\Users\arielv\Desktop\תואר שני\Speckle-fibers\pulseLaser\100_HG_modes.mat').Psi;
n = 10; %how many modes to use
%% Normalize the modes patterns with L2 norm
PsiNorm = Psi;
for ii=1:n
    for jj=1:n
        PsiNorm(:,:,ii,jj) = PsiNorm(:,:,ii,jj)/norm(PsiNorm(:,:,ii,jj),2);
        sums(ii,jj) = sum(abs(PsiNorm(:,:,ii,jj)).^2,'all');
    end
end


%% present the modes
fig = figure(1); 
for ii=1:n
     for jj=1:n
         subplot(n,n,(ii-1)*n+jj);
         imagesc(PsiNorm(:,:,ii,jj)); impixelinfo;
         %title(['HG_',num2str(ii-1),'_',num2str(jj-1),' Mode']);
         axis off
     end
end

%% try and present different modes combinations
Nmodes = 100; % number of modes
Tmodes = 5; % optional modes in each axis
shifts = 10;
modes = randi(Tmodes,[Nmodes,2]); %random Nmodes indexes
PsiTot = zeros(N); %initiate total modes
for nn=1:Nmodes
    PsiTot = PsiTot + rand(1)*imrotate(circshift(PsiNorm(:,:,modes(nn,1),modes(nn,2)).*exp(1i*ones(N)*(2*rand(1)-1)*pi),...
        [randi(shifts)-shifts randi(shifts)-shifts]),360*rand(1),'bilinear','crop'); % add all modes together, random phase, random shifts, random rotations
    PsiTot = PsiTot./norm(PsiTot,1);
end
figure(2);
imagesc(abs(PsiTot)/max(abs(PsiTot(:)))), colorbar, impixelinfo;
title(['Intensity pattern of combination of ',num2str(Nmodes),' modes (not nessecarily different) out of ',num2str(Tmodes^2), ' optional modes']);
