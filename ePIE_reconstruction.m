function [object, probe, err] =ePIE_reconstruction(I_fft, cen, probe_init, probe_init_pos, params, obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ePIE (extended PtychographIcal Engine) ITERATIVE PHASE RETRIEVAL ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is based on the ePIE itertaive phase retrieval 
% algorithm as described by Andrew M.Maiden and John M.Rodenburg in 
% "An improved ptychographical phase retrieval algorithm for diffractive
% imaging", Ultramicroscopy 109 (2009) 1256–1262.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Ariel Veler, July-2021
% inputs:
% I_fft - input diffraction patterns of probes and object (magnitude square)
% cen - centers of the probes at image plane ([n(1)x n(2)x 2] double)
% obj - object at image plane (amplitude and phase) for comparison
% Eps - allowed total erroe in reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iterations = params.iterations; % number of iterative Iterations
Eps = params.eps;
pixel_pitch = params.pixel_pitch;
lambda = params.lambda;
f2 = params.f2;
s = size(I_fft(:,:,1,1)); 
[nrow,ncol] = size(cen(:,:,1)); % number of probes [row,col]
%h1 = [round(s(1)/2-s(1)/5):round(s(1)/2+s(1)/5)]; h2 = [round(s(2)/2-s(2)/5):round(s(2)/2+s(2)/5)];
alpha = 1; beta = 0;

% create a grid for drawings:
dkx = pixel_pitch(1)/(lambda*f2); dky = pixel_pitch(2)/(lambda*f2);
dx = 1/(s(2)*dkx); dy = 1/(s(1)*dky);
x = -s(2)/2*dx:dx:s(2)/2*dx-dx; y = -s(1)/2*dy:dy:s(1)/2*dy-dy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make stack of probes
amp_fft = sqrt(I_fft);
cen = reshape(cen,[size(cen,1)*size(cen,2),size(cen,3)]);
amp_fft = reshape(amp_fft, [size(amp_fft,1),size(amp_fft,2),size(amp_fft,3)*size(amp_fft,4)]);
h1 = round(cen(1,1):cen(end,1)); h2 = round(cen(1,2):cen(end,2));

% sig = EstimateRadius(amp_fft,0.1)./sqrt(2*log(2)); % initial guess of sigma of the probe - FWHM/sqrt(2*log(2))
% p = eval(params.LoadProbe);
probe = probe_init; % Probe initial guess
probe_support = find_probe_support(probe);

% object initialization from the diffraction pattern:
object = ifft2(ifftshift(amp_fft(:,:,ceil(size(amp_fft,3)/2))))./(ShiftProbefft(probe,cen(ceil(size(cen,1)/2),:)-probe_init_pos)+0.01); %guess - center probe diffraction pattern

% initiate for reconstruction: 
psi_real = zeros(size(amp_fft));
psi_fft = zeros(size(amp_fft));
psi_real_tmp = zeros(size(amp_fft));

% support:
support = zeros(s);
for jj=1:size(cen,1)
    support = support+abs(ShiftProbefft(probe_support,cen(jj,:)-probe_init_pos)); % support - FWHM
end
support(support>0.5)=1;support(support<=0.5)=0;

error_step = 2*Eps;
%general step:
iter=1;
figure('position',get(0,'ScreenSize'))
while iter < Iterations % && error_step > Eps  
    fprintf('Iteration: %d\n', iter);
    
    % Shuffeling:
    pnum = randperm(size(cen,1));
    c = cen(pnum,:);
    a_fft = amp_fft(:,:,pnum);
    
    for ii=1:numel(pnum) 
    % apply 4-steps of algorithm
        probe = ShiftProbefft(probe,c(ii,:)-probe_init_pos);
        psi_real(:,:,ii) = probe.*object;
        psi_fft(:,:,ii) = fftshift(fft2(psi_real(:,:,ii))); % step 1 - fft
        psi_fft(:,:,ii) = a_fft(:,:,ii).*psi_fft(:,:,ii)./abs(psi_fft(:,:,ii)); % step 2 - constraint k-space
        psi_real_tmp(:,:,ii) = ifft2(ifftshift(psi_fft(:,:,ii))); % step 3 - ifft
        % object update
        object_new = object + alpha*(conj(probe)./max(max(abs(probe).^2))).*(psi_real_tmp(:,:,ii)-psi_real(:,:,ii));
        % probe update
        probe_new = probe + beta*(conj(object)./max(max(abs(object).^2))).*(psi_real_tmp(:,:,ii)-psi_real(:,:,ii));
        
        if any(isnan([abs(object_new(:)),angle(object_new(:)),abs(probe_new(:)),angle(probe_new(:))]))
            error('NaN value');
        end

        object = object_new.*support;
        probe = ShiftProbefft(probe_new,-c(ii,:)+probe_init_pos).*probe_support;
        
    
    end
    if iter == 10 
        beta = 1;
    end
    
    % Calculate error - Compare to original object
    if exist('obj','var')
        gamma = sum(obj(h1,h2).*conj(object(h1,h2)),'all')./sum(abs(object(h1,h2)).^2,'all');
    %     error_step = sum(abs(obj(h1,h2)-gamma.*object(h1,h2))^2,'all')./sum(abs(obj(h1,h2)).^2,'all');
        error_step = sum(abs(abs(obj(h1,h2))-abs(gamma).*abs(object(h1,h2))).^2,'all')./sum(abs(obj(h1,h2)).^2,'all');
        err(iter) = error_step; % count error every iteration
    else
        error_step = 0; err(iter) = 0;
    end

    object_phase = angle(object);
    object_magnitude = abs(object);
    probe_phase = angle(probe);
    probe_magnitude = abs(probe);
    
    if any(isnan([object_magnitude(:),object_phase(:),probe_magnitude(:),probe_phase(:)]))
        error('NaN value');
    end

    if iter<10 || mod(iter,10)==0
        subplot(2,3,1);
        imagesc(x,y,object_magnitude), colorbar, title('Object Magnitude','FontSize',16); axis tight; xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
        subplot(2,3,2);
        imagesc(x,y,object_phase), colorbar, title('Object Phase','FontSize',16); axis tight; xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
        subplot(2,3,4);
        imagesc(x,y,probe_magnitude), colorbar, title('Probe Magnitude','FontSize',16); axis tight; xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
        subplot(2,3,5);
        imagesc(x,y,probe_phase), colorbar, title('Probe Phase','FontSize',16); axis tight; xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
        subplot(2,3,3);
        plot(err), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error'); %show error vs iterations
        subplot(2,3,6);
        scatter(cen(:,2),cen(:,1),'x'); axis([1 s(2) 1 s(1)]); %show centers position
        sgtitle(['Sample progress - iter: ', num2str(iter), ' error: ' , num2str(error_step)],'FontSize',24);
        impixelinfo;
        drawnow;
    end

    iter = iter+1;
end
end
% return reconstruction and number of iteratiom
%Phase_image = angle(z_real);

function Pout = ShiftProbefft(Pin,c)
% Image shift throuth fft and linear phase in c pixels
% c - pixels to move ([1x2] double)
    s = size(Pin);    
    sr = 0:s(1)-1;
    sc = 0:s(2)-1;
    [sc,sr] = meshgrid(sc,sr);
    Pout = ifft2(ifftshift(fftshift(fft2(Pin)).*exp(-1i*2*pi*(c(1)*sr/s(1)+c(2)*sc/s(2)))));
       
end

function Pout = ShiftProbe(Pin,c) 
% shift center of probe to "c", when image size "s" and center of Pin is [1,1]
    Pout = circshift(Pin,flip(c));
end
 
function [row,col] = orig_ind(r,n)
    row = mod(n,r);
    col = ceil(n/r);
end

function r = EstimateRadius(amp_fft,thr)
    r = mean(sqrt(sum(ifft2(ifftshift(amp_fft))>thr,[1 2])/pi));
end

function support = find_probe_support(probe)
    %support = (probe > max(probe(:))/exp(0.5)) + (probe < -max(probe(:))/exp(0.5)); % Sigma
    support = (probe > max(probe(:)).*0.05) + (probe < -max(probe(:)).*0.05);
    s = size(probe);
    [x,y] = meshgrid(1:s(2),1:s(1));
    dia1 = max(sum(support,1))-min(sum(support,1)); % diameter of active columns
    dia2 = max(sum(support,2))-min(sum(support,2)); % diameter of active rows
    support = (x-s(2)/2).^2./(dia1/2).^2+(y-s(1)/2).^2./(dia2/2).^2<=1;
end