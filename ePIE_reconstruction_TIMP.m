function [object, probe, err] =ePIE_reconstruction_TIMP(I_fft, cen, probe_init, probe_init_pos, obj, Eps, params)
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
% cen - centers of the probes at object plane (Fourier plane of the 4f system) at calibrated distance to real space 
% probe_init - initial guess for all probes by there order [s(1)x s(2)x k]
% probe_init_pos - initial positions for all probes [1 x 2] (usually center of the image)
% obj - objects at image plane (amplitude and phase) for comparison
% Eps - allowed total erroe in reconstruction
% params - optional parametens to use (e.g. data on the probes, no. of iteration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('params', 'var')
    Iterations = params.iterations; % number of iterative Iterations
else
    Iterations = 50;
end
s = size(I_fft(:,:,1,1)); 
[nrow,ncol] = size(cen(:,:,1)); % number of probes [row,col]
alpha = 1; beta = 0; % controll of convergence rate of the objects\probes respectively
k_num = size(probe_init,3); % no. of initial guess for all probes. 
% Here we assume that no. of objects equal to no. of probes!!! meaning
% there are k_nim different probes matching to k_num different objects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make stack of amplitude images (measurement images) and centers
amp_fft = sqrt(I_fft);
amp_fft = reshape(amp_fft, [size(amp_fft,1),size(amp_fft,2),size(amp_fft,3)*size(amp_fft,4)]);
cen = reshape(cen,[size(cen,1)*size(cen,2),size(cen,3)]);
h1 = round(cen(1,1):cen(end,1)); h2 = round(cen(1,2):cen(end,2)); % choose the support for calculating the error to be the square boundary by the centers 

% sig = EstimateRadius(amp_fft,0.1)./sqrt(2*log(2)); % initial guess of sigma of the probe - FWHM/sqrt(2*log(2))
% p = eval(params.LoadProbe);
probe = probe_init; % Probe initial guess
probe_support = find_probe_support(probe);

% object initialization from the diffraction pattern:    ************ check ******************
for k=1:k_num
    object(:,:,k) = ifft2(ifftshift(amp_fft(:,:,ceil(size(amp_fft,3)/2))))./(ShiftProbefft(probe(:,:,k),cen(ceil(size(cen,1)/2),:)-probe_init_pos)+0.01)/sqrt(k_num); %guess - center probe diffraction pattern
end

% support:
support = zeros(s(1),s(2),k_num);
for k=1:k_num
    for jj=1:size(cen,1)
        support(:,:,k) = support(:,:,k)+abs(ShiftProbefft(probe_support(:,:,k),cen(jj,:)-probe_init_pos)); % support - FWHM
    end
end
support(support>0.5)=1;support(support<=0.5)=0;


%general step:
iter=1;
figure('position',get(0,'ScreenSize'))
while iter < Iterations %&& sqrt(sum((abs(z_fft)-amp_fft).^2,'all')/sum(amp_fft.^2,'all')) > Eps  
    fprintf('Iteration: %d\n', iter);
    % Shuffeling:
    pnum = randperm(size(cen,1));
    c = cen(pnum,:);
    a_fft = amp_fft(:,:,pnum);

    for ii=1:numel(pnum) 
        probe = ShiftProbefft(probe,c(ii,:)-probe_init_pos);
        psi_real = probe.*object;

        % apply 4-steps of algorithm   
        psi_fft = fftshift(fft2(psi_real)); % step 1 - fft
        psi_fft = a_fft(:,:,ii).*psi_fft./sqrt(sum(abs(psi_fft).^2,3)); % step 2 - constraint k-space
        psi_real_tmp = ifft2(ifftshift(psi_fft)); % step 3 - ifft
        % object update
        object_new = object + alpha*(conj(probe)./max(max(abs(probe).^2))).*(psi_real_tmp-psi_real);
        % probe update
        probe_new = probe + beta*(conj(object)./max(max(abs(object).^2))).*(psi_real_tmp-psi_real);
        
        if any(isnan([abs(object_new(:)),angle(object_new(:)),abs(probe_new(:)),angle(probe_new(:))]))
            error('NaN value');
        end

        object = object_new.*support;
        probe = ShiftProbefft(probe_new,-c(ii,:)+probe_init_pos).*probe_support;

    end
    if iter == 10 
        beta = 1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%% continue from here - error and drawing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate error - Compare to original object
    gamma = sum(obj(h1,h2,:).*conj(object(h1,h2,:)),[1,2])./sum(abs(object(h1,h2,:)).^2,[1,2]);
    error_step = sum(abs(obj(h1,h2,:))-abs(gamma).*abs(object(h1,h2,:)).^2,[1,2])./sum(abs(obj(h1,h2)).^2,[1,2]);
    err(iter,:) = error_step; % count error every iteration


    object_phase = angle(object);
    object_magnitude = abs(object);
    probe_phase = angle(probe);
    probe_magnitude = abs(probe);

    if any(isnan([object_magnitude(:),object_phase(:),probe_magnitude(:),probe_phase(:)]))
        error('NaN value');
    end

%     if iter<10 || mod(iter,10)==0
%         subplot(2,3,1);
%         imshow(object_magnitude, []), colorbar, title('Object Magnitude','FontSize',16); %show current phase-image
%         subplot(2,3,2);
%         imshow(object_phase, []), colorbar, title('Object Phase','FontSize',16); %show current phase-image
%         subplot(2,3,4);
%         imshow(probe_magnitude, []), colorbar, title('Probe Magnitude','FontSize',16); %show current phase-image
%         subplot(2,3,5);
%         imshow(probe_phase, []), colorbar, title('Probe Phase','FontSize',16); %show current phase-image
% %         subplot(2,3,3);
% %         plot(err), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error'); %show error vs iterations
%         subplot(2,3,6);
%         scatter(cen(:,2),cen(:,1),'x'); axis([1 s(2) 1 s(1)]); %show centers position
%         sgtitle(['Sample progress - iter: ', num2str(iter), ' error: ' , num2str(error_step)],'FontSize',24);
%         impixelinfo;
%         drawnow;
%     end
    iter = iter+1;
end
end
% return reconstruction and number of iteratiom
%Phase_image = angle(z_real);

function Pout = ShiftProbefft(Pin,c)
% Image shift throuth fft and linear phase in c pixels
% Pin - input image (may be 3D) 
% c - pixels to move ([1x2] double)
    for k=1:size(Pin,3)
        s = size(Pin(:,:,k));    
        sr = 0:s(1)-1;
        sc = 0:s(2)-1;
        [sc,sr] = meshgrid(sc,sr);
        Pout(:,:,k) = ifft2(ifftshift(fftshift(fft2(Pin(:,:,k))).*exp(-1i*2*pi*(c(1)*sr/s(1)+c(2)*sc/s(2)))));
    end     
end

function support = find_probe_support(probe)
% find support for each prob initial pattern
% input:
% probe - a [m x n x k] matrice, when k is the number of different
% probes\objects in the system. each [n x m] image is an initial guess for
% the probes used in the propagation.
    for k=1:size(probe,3)
        %support = (probe > max(probe(:))/exp(0.5)) + (probe < -max(probe(:))/exp(0.5)); % Sigma
        support_k = (probe(:,:,k) > max(max(probe(:,:,k))).*0.05) + (probe(:,:,k) < -max(max(probe(:,:,k))).*0.05); % 95% of max value used
        s = size(probe(:,:,k));
        [x,y] = meshgrid(1:s(2),1:s(1));
        dia1 = max(sum(support_k,1))-min(sum(support_k,1)); % diameter of active columns
        dia2 = max(sum(support_k,2))  -min(sum(support_k,2)); % diameter of active rows
        support(:,:,k) = (x-s(2)/2).^2./(dia1/2).^2+(y-s(1)/2).^2./(dia2/2).^2<=1; % an eliptical support for the k'th probe according to active area
    end
end