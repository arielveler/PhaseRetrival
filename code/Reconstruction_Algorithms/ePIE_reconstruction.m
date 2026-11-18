function [object, probe, err] =ePIE_reconstruction(I_fft, cen, probe_init, params, obj)
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

fields = fieldnames(params);
for i=1:numel(fields)
    eval([fields{i},' = params.',fields{i},';']);
end
% must have fields:
% iterations, Eps
% optional fields (for drawing):
% pixel_pitch, lambda, f1, f2, b 

s = size(I_fft(:,:,1,1)); 
if ~exist('alpha', 'var')  alpha = 1; end
if ~exist('beta', 'var')  beta = 0; end
if ~exist('change_iter', 'var')  change_iter = 0; end
if ~exist('support_thr', 'var')  support_thr = 0.05; end

if exist('f2','var')
    % definition of the reciprocoral and real grid after scaling by the f2 lens:
    dk = pixel_pitch(1)/(lambda*f2);
    
    % create a grid for drawings:
    % dx = pixel_pitch(1); dy = pixel_pitch(2);
    % factor = [(b*dx)/(lambda*f2), (b*dy)/(lambda*f2)];
    % x = dx*(-s(2)/2:s(2)/2-1)/factor(1); y = dy*(-s(1)/2:s(1)/2-1)/factor(2);
    dx_im = 1/(s(2)*dk); dy_im = 1/(s(1)*dk);
    x_vec_im = dx_im*(-s(2)/2:s(2)/2-1); y_vec_im = dy_im*(-s(1)/2:s(1)/2-1);
    [x_im,y_im] = meshgrid(x_vec_im, y_vec_im);
else
    x_vec_im = 1:s(2); y_vec_im = 1:s(1);
end
if ~exist('probe_init_pos','var')
    probe_init_pos = floor(size(probe_init)/2)+1;
end

if ~exist('is_support', 'var')  is_support = 0; end % if we want to activate support on region every iteration
if ~exist('show_internal', 'var')  show_internal = 0; end % if we want to see in zoom the object region during iterations
is_GPU = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make stack of amplitude images (measurement images) and centers
amp_fft = sqrt(I_fft);
amp_fft = reshape(amp_fft, [size(amp_fft,1),size(amp_fft,2),size(amp_fft,3)*size(amp_fft,4)]);
cen = reshape(cen,[size(cen,1)*size(cen,2),size(cen,3)]);
h1 = round(cen(1,1):cen(end,1)); h2 = round(cen(1,2):cen(end,2));

% Probe initial guess - input parameter:
probe = probe_init; 
probe_support = find_probe_support(probe,support_thr); 

if exist('obj_guess','var')
    object = obj_guess;
else
% object initial guess - initialization from the diffraction pattern:
% object = ifftshift(ifft2(ifftshift(amp_fft(:,:,ceil(size(amp_fft,3)/2)))))./(ShiftProbefft(probe,cen(ceil(size(cen,1)/2),:)-probe_init_pos)+0.01); %guess - center probe diffraction pattern
    object = zeros(size(amp_fft(:,:,1)))+0.000000000001; %object(ceil((size(object,1)+1)/2),ceil((size(object,2)+1)/2)) = 0.0000001;
%     object = zeros(size(amp_fft(:,:,1)))+0.01; %object(ceil((size(object,1)+1)/2),ceil((size(object,2)+1)/2)) = 0.0000001;
end


% initiate for reconstruction: 
psi_real = zeros(size(amp_fft));
psi_fft = zeros(size(amp_fft));
psi_real_tmp = zeros(size(amp_fft));

% support:
support = zeros(s);
for jj=1:size(cen,1)
    support_tmp = abs(ShiftProbefft(probe_support,cen(jj,:)-probe_init_pos));% support - FWHM
    support = support + support_tmp>0.5;
end
support(support>0.5)=1;support(support<=0.5)=0; 

if ~is_support
    probe_support = ones(size(probe_support)); %%% check if needed
    support = ones(size(support)); %%% check if needed
end

if is_GPU
    [amp_fft, object, probe, cen] = GPU_up(amp_fft, object, probe, cen);
end

%general step:
iter=1;
figure('position',get(0,'ScreenSize'))

while iter < iterations % && error_step > Eps  
    if mod(iter,100) == 0 fprintf('Iteration: %d\n', iter); end
    
    % Shuffeling:
    pnum = randperm(size(cen,1));
    c = cen(pnum,:);
    a_fft = amp_fft(:,:,pnum);
    
    for ii=1:numel(pnum) 
        probe = ShiftProbefft(probe,c(ii,:)-probe_init_pos);
        psi_real = probe.*object;
        % apply 4-steps of algorithm
        psi_fft = fft2(ifftshift(psi_real)); % step 1 - fft
        psi_fft = ifftshift(a_fft(:,:,ii)).*exp(1i*angle(psi_fft)); % step 2 - constraint k-space
        psi_real_tmp = fftshift(ifft2(psi_fft)); % step 3 - ifft
        % object update
        object_new = object + alpha*(conj(probe)./max(max(abs(probe).^2))).*(psi_real_tmp-psi_real);
        % probe update
        probe_new = probe + beta*(conj(object)./max(max(abs(object).^2))).*(psi_real_tmp-psi_real);
        
        if any(isnan([abs(object_new(:)),angle(object_new(:)),abs(probe_new(:)),angle(probe_new(:))]))
            error('NaN value');
        end

%         object = exp(1i*angle(object_new));
        object = object_new.*support;
        probe = ShiftProbefft(probe_new,-c(ii,:)+probe_init_pos).*probe_support;
        
    
    end
%     Correcting probes-centers (Gradient Decend):
%     option 1 - only 2 convolutions:
%     dOdx = diff(O)./dx;
%     dOdy = diff(O')'./dy;
%     circshift(dOdx,m,n)
%
%     option 2 - n*m convolutions:
%     temp = 0;
%     Dx = [-1,1;-1,1]; Dy = [-1,-1;1,1];
%     gamma = 1;
%     for m,n
%       probe = ShiftProbefft(probe,c(ii,:)-probe_init_pos);
%       psi_real = probe.*object;
%       psi_fft = fft2(ifftshift(psi_real)); % step 1 - fft
%       psi_fft = ifftshift(a_fft(:,:,ii)).*exp(1i*angle(psi_fft)); % step 2 - constraint k-space
%       psi_real_tmp = fftshift(ifft2(psi_fft)); % step 3 - ifft
%       theta = atan2(n,m); 
%       grad = Dx*cos(theta) + Dy*sin(theta); 
%       gradO = conv2(O,grad,'same')*sqrt(m^2+n^2);
%       temp = temp + tan(atand(b/f1))/2*sum(conj(object.*probe-psi_real_tmp).*probe.*gradO)  *{gamma/(max(abs(object))*max(abs(probe)))};
%      end
%      d_new = d + temp;
%      update c accordign to d_new, ang=atand(b/f1)
%
% option 3:
% inside the while loop:
%      temp = 0;
%      Dx = [-1,1;-1,1]; Dy = [-1,-1;1,1];
%      gamma = 1;
%      for ii=1:numel(pnum)
%          probe = ShiftProbefft(probe,c(ii,:)-probe_init_pos);
%          psi_real = probe.*object;
%          psi_fft = fft2(ifftshift(psi_real)); % step 1 - fft
%          psi_fft = ifftshift(a_fft(:,:,ii)).*exp(1i*angle(psi_fft)); % step 2 - constraint k-space
%          psi_real_tmp = fftshift(ifft2(psi_fft)); % step 3 - ifft
%          phi = atan(imag(psi_fft)./real(psi_fft));
%          d_phi(x/y) = 1i*k_x/y.*abs(fft2(ifftshift(probe_init.*object))).^2./(real(psi_fft)).^2;
%          c_new(ii,x/y) = c(ii,x/y) + sum(conj(object.*probe-psi_real_tmp).*(probe.*diff(object,1/2)+1i*fftshift(ifft2(a_fft(:,:,ii).*d_phi(x/y).*exp(1i*phi)))   *{gamma/(max(abs(object))*max(abs(probe)))};
%      end
    if iter == 30 && change_iter
        beta = 0.2;
    end
    
    % Calculate error - Compare to original object
    if exist('obj','var') % simulation - there is an object to compare
        gamma = sum(obj(h1,h2).*conj(object(h1,h2)),'all')./sum(abs(object(h1,h2)).^2,'all');
%         error_step = sum(abs(obj(h1,h2)-gamma.*object(h1,h2))^2,'all')./sum(abs(obj(h1,h2)).^2,'all'); % error from all - same ass paper
        error_step = sum(abs(abs(obj(h1,h2))-abs(gamma).*abs(object(h1,h2))).^2,'all')./sum(abs(obj(h1,h2)).^2,'all'); %error only from absolute value
        err(iter) = error_step; % count error every iteration
    else % real object - no object to compare - error using the original diffraction patterns as ground-trouth.
        for jj=1:size(amp_fft,3)
            psi_fft_comp(:,:,jj) = fftshift(fft2(object.*ShiftProbefft(probe,cen(jj,:)-probe_init_pos)));  %%%check fftshift-fft2 order%%
        end
%         error_step = sum(sum(abs(amp_fft(h1,h2,:)-abs(psi_fft_comp(h1,h2,:))).^2,[1 2]))./sum(sum(abs(amp_fft(h1,h2,:)).^2,[1 2])); 
        error_step = sum(sqrt(sum(abs(amp_fft-abs(psi_fft_comp)).^2,[1 2]))./sqrt(sum(abs(amp_fft).^2,[1 2])))./numel(pnum); 
        err(iter) = error_step;
    end

    object_phase = angle(object);
    object_magnitude = abs(object).^2; % intensity
    probe_phase = angle(probe);
    probe_magnitude = abs(probe).^2; % intensity, for comparison
    
    if any(isnan([object_magnitude(:),object_phase(:),probe_magnitude(:),probe_phase(:)]))
        error('NaN value');
    end

    if mod(iter,10)==0  %iter<10 || mod(iter,10)==0
        color_map = 'gray';
        if show_internal
            subplot(2,3,1);
            imagesc(x_vec_im(sum(support,1)>0),y_vec_im(sum(support,2)>0),object_magnitude(sum(support,1)>0,sum(support,2)>0)), colorbar, title('Object Magnitude','FontSize',16); axis tight; colormap(color_map); xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
            subplot(2,3,2);
            imagesc(x_vec_im(sum(support,1)>0),y_vec_im(sum(support,2)>0),object_phase(sum(support,1)>0,sum(support,2)>0)), colorbar, title('Object Phase','FontSize',16); axis tight; colormap(color_map); xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
            subplot(2,3,4);
            imagesc(x_vec_im(sum(support,1)>0),y_vec_im(sum(support,2)>0),probe_magnitude((sum(support,1)>0),(sum(support,2)>0))), colorbar, title('Probe Magnitude','FontSize',16); axis tight; colormap(color_map); xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
            subplot(2,3,5);
            imagesc(x_vec_im(sum(support,1)>0),y_vec_im(sum(support,2)>0),probe_phase((sum(support,1)>0),(sum(support,2)>0))), colorbar, title('Probe Phase','FontSize',16); axis tight; colormap(color_map); xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
        else
            subplot(2,3,1);
            imagesc(x_vec_im,y_vec_im,object_magnitude.*support), colorbar, title('Object Magnitude','FontSize',16); axis tight; colormap(color_map); xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
            subplot(2,3,2);
            imagesc(x_vec_im,y_vec_im,object_phase), colorbar, title('Object Phase','FontSize',16); axis tight; colormap(color_map); xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
            subplot(2,3,4);
            imagesc(x_vec_im,y_vec_im,probe_magnitude), colorbar, title('Probe Magnitude','FontSize',16); axis tight; colormap(color_map); xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
            subplot(2,3,5);
            imagesc(x_vec_im,y_vec_im,probe_phase), colorbar, title('Probe Phase','FontSize',16); axis tight; colormap(color_map); xlabel('Size [m]'); ylabel('Size [m]'); %show current phase-image
        end
        subplot(2,3,3);
        plot(err(11:end)), title('error Vs iterations'), xlabel('Iterations'), ylabel('Error'); %show error vs iterations
        subplot(2,3,6);
        scatter(cen(:,2),cen(:,1),'x'); axis([1 s(2) 1 s(1)]); %show centers position
        sgtitle(['iter: ', num2str(iter), ' error: ' , num2str(error_step)],'FontSize',24);
        impixelinfo;
        drawnow;
    end
    iter = iter+1;
end
end

% return reconstruction and number of iteratiom
function Pout = ShiftProbefft(Pin,c)
% Image shift throuth fft and linear phase in c pixels
% c - pixels to move ([1x2] double)
    s = size(Pin);  
    sr = ifftshift(-fix(s(1)/2):ceil(s(1)/2)-1);
    sc = ifftshift(-fix(s(2)/2):ceil(s(2)/2)-1);
    [sc,sr] = meshgrid(sc,sr);

    Pout = ifft2(fft2(Pin).*exp(-1i*2*pi*(c(1)*sr/s(1)+c(2)*sc/s(2))));
       
end

function varargout = GPU_up(varargin)
     for i=1:nargin
         varargout{i} = gpuArray(varargin{i});
     end
end

% function support = find_probe_support(probe) % old version
%     %support = (probe > max(probe(:))/exp(0.5)) + (probe < -max(probe(:))/exp(0.5)); % Sigma
%     support = (probe > max(probe(:)).*0.05) + (probe < -max(probe(:)).*0.05);
%     s = size(probe);
%     [x,y] = meshgrid(1:s(2),1:s(1));
%     dia1 = max(sum(support,1))-min(sum(support,1)); % diameter of active columns
%     dia2 = max(sum(support,2))-min(sum(support,2)); % diameter of active rows
%     support = (x-s(2)/2).^2./(dia1/2).^2+(y-s(1)/2).^2./(dia2/2).^2<=1;
% end

function support = find_probe_support(probe,thr)
% find support for each prob initial pattern
% input:
% probe - a [m x n x k] matrice, when k is the number of different
% probes\objects in the system. each [n x m] image is an initial guess for
% the probes used in the propagation.
    if ~exist('thr','var') thr = 0.05; end % default of 5%
    for k=1:size(probe,3)
        %support = (probe > max(probe(:))/exp(0.5)) + (probe < -max(probe(:))/exp(0.5)); % Sigma
        support = (probe > max(max(abs(probe))).*thr) + (probe < -max(max(abs(probe))).*thr); % thresholding of max value
        s = size(probe);
        [x,y] = meshgrid(1:s(2),1:s(1));
        dia1 = find(sum(support,1),1,"last")-find(sum(support,1),1); % diameter of active columns
        dia2 = find(sum(support,2),1,"last")-find(sum(support,2),1); % diameter of active rows
        support = (x-s(2)/2).^2./(dia1/2).^2+(y-s(1)/2).^2./(dia2/2).^2<=1; % an eliptical support for the k'th probe according to active area
    end
end