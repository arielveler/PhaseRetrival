function [Phase_image, err] = GS_reconstruction(amp_real, I_fft, Eps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GERCHBERG-SAXTON ITERATIVE PHASE RETRIEVAL ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is based on the Gerchberg-Saxton itertaive phase retrieval 
% algorithm as described by R. W. Gerchberg and W. O. Saxton in 
% "A practical algorithm for the determination of phase from image and 
% diffraction plane pictures" Optik 35(2), pages 237 - 246 (1972) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Ariel Veler, June-2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iterations = 1000; % number of iterative Iterations, typically 200 iterations are enough
p = 0.0001;           % time to pause, otherwise images are not shown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = size(amp_real);
if s~=size(I_fft)
    error('not the same size');
end
%initialization:
amp_fft = sqrt(I_fft);
z_real = amp_real.*exp(1i*(2*rand(s)-1)*pi); % random phase between -pi:pi
z_fft = fftshift(fft2(z_real)); % step 1 - fft

%general step:
iter=1;
figure('position',[100 100 2400 1200]); 
while sum(abs(abs(z_fft)-amp_fft).^2,'all')/sum(amp_fft.^2,'all') > Eps && iter < Iterations
    
    fprintf('Iteration: %d\n', iter);  
    error_step = sqrt(sum(abs(abs(z_fft)-amp_fft).^2,'all')/sum(amp_fft.^2,'all'));
    err(iter) = error_step; % count error every iteration

    z_fft = amp_fft.*z_fft./abs(z_fft); % step 2 - constraint k-space
    z_real = ifft2(ifftshift(z_fft)); % step 3 - ifft
    
    z_phase = angle(z_real);

    if mod(iter,10)==0
        subplot(1,2,1);
        imshow(abs(z_real), []), colorbar, title('Magnitude','FontSize',16); %show current phase-image
        subplot(1,2,2);
        imshow(z_phase, []), colorbar, title('Phase','FontSize',16); %show current phase-image
        sgtitle(['Sample progress - iter: ', num2str(iter), ' error: ' , num2str(error_step)],'FontSize',24);
        impixelinfo;
    end
    
    z_real = amp_real.*z_real./abs(z_real); % step 4 - constraint real-space
    iter = iter+1;
    pause(p);
    z_fft = fftshift(fft2(z_real)); % step 1 - fft

end
% return reconstruction and number of iteratiom
Phase_image = angle(z_real);