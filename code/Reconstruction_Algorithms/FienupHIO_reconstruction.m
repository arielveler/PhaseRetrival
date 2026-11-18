function [Phase_image, err] = FienupHIO_reconstruction(I_fft, beta, supp, Eps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIENUP-HIO (HYBRID INPUT OUTPUT) ITERATIVE PHASE RETRIEVAL ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is based on the Fienup-HIO itertaive phase retrieval 
% algorithm as described by J. Fienup in "Phase retrieval algorithms: 
% A comparison", Appl. Opt., vol. 21, no. 15, pp. 2758–2769 (1982).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Ariel Veler, June-2021
% inputs:
% I_fft - input diffraction pattern (magnitude square)
% beta - the parameter for the real-space constraint step, usually 0.9.
% supp -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iterations = 1000; % number of iterative Iterations, typically 200 iterations are enough
p = 0.00001;           % time to pause, otherwise images are not shown
s = size(I_fft); 
div_par = -0.5*supp+0.5; % how much of image there is in every axis
% support = zeros(s); support(s(1)*div_par+1:end-s(1)*div_par,s(2)*div_par+1:end-s(2)*div_par)=1;
support = ones(s).*zeros(s); support(ceil(s(1)*div_par)+1:end-ceil(s(1)*div_par),ceil(s(2)*div_par)+1:end-ceil(s(2)*div_par))=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialization:
amp_fft = sqrt(I_fft);
z_fft = amp_fft.*exp(1i*(2*rand(s)-1)*pi); % random phase between -pi:pi
z_real_tmp = ifft2(ifftshift(z_fft)); % real space image
z_real = z_real_tmp.*support; % step 4 - constraint real-space - support HIO
z_fft = fftshift(fft2(z_real)); % step 1 - fft

%general step:
iter=1;
figure('position',[100 100 2400 1200]); 
while iter < Iterations && sqrt(sum((abs(z_fft)-amp_fft).^2,'all')/sum(amp_fft.^2,'all')) > Eps  
    
    fprintf('Iteration: %d\n', iter);  
    
    
    % apply 4-steps of algorithm
    z_fft = amp_fft.*z_fft./abs(z_fft); % step 2 - constraint k-space
    z_real_tmp = ifft2(ifftshift(z_fft)); % step 3 - ifft
    z_real = z_real_tmp.*support + (z_real-beta*z_real_tmp).*(1-support);  % step 4 - constraint real-space - support HIO
    z_fft = fftshift(fft2(z_real)); % step 1 - fft
    
    error_step = sqrt(sum(abs(abs(z_fft)-amp_fft).^2,'all')/sum(amp_fft.^2,'all'));
    err(iter) = error_step; % count error every iteration
        
        z_phase = angle(z_real);
        z_abs = abs(z_real);

        if mod(iter,50)==0
            subplot(1,2,1);
            imshow(z_abs, []), colorbar, title('Magnitude','FontSize',16); %show current phase-image
            subplot(1,2,2);
            imshow(z_phase, []), colorbar, title('Phase','FontSize',16); %show current phase-image
            sgtitle(['Sample progress - iter: ', num2str(iter), ' error: ' , num2str(error_step)],'FontSize',24);
            impixelinfo;
        end
    
    iter = iter+1;
    pause(p);

end
% return reconstruction and number of iteratiom
Phase_image = angle(z_real);