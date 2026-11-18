function [cen,subimage,angle_tilt,b] = FindRefCenters(im,n,Thr,Dpix)
% function to find centers of probes, indexes of cutting to subimages and
% angle of image rotation.
% inputs:
% im - probe calibration image (taken without object)
% n - number of probes ([1x2] int)
% Dpix - extimation of probe diameter FWHM in pixels (double) 
% Thr - Threshold division from maximum image value (double)
% outputs:
% cen - centers of image in subpixel measurement ([n(1) x n(2) x2] double)
% subimage - size of the subimages to cut ([1x2] int)
% angle_tilt - tilted angle (in degrees) of the data relative to the vertical axis
% b - calculated distance between adjacent probes vertically and
% horizontally seperatelly in subpixel manner ([1x2] double)

% created by Ariel Veler, Jan-2022

    % default values
    if ~exist('n','var'); n=[5 5]; end 
    if ~exist('Thr','var'); Thr=10; end 
    if ~exist('Dpix','var'); Dpix=0; end 

    Grid = n(1)*n(2);
    GW = kron(gausswin(Dpix), gausswin(Dpix)'); % filter with a gaussian window
    ref = conv2(im, GW, 'same')/sum(GW(:));
    thresh = max(ref(:))/Thr;
    flag = 1;

    while flag % make sure circN is higher or equal to Grid   
        ref_tmp = ref;
        ref_tmp(ref_tmp<thresh) = 0;
        BW = imbinarize(gather(ref_tmp), gather(thresh));
        s = regionprops(BW,ref,'WeightedCentroid');
        circN = size(s,1);
    
        if circN <Grid
%            error('Threshold is too high. try lower threshold');
            thresh = thresh*0.95;
            continue
        end
        flag = 0;
    end

    
    while circN>=Grid 
        if circN==Grid
            tmp_s = s;
            break
        end
        thresh = thresh*1.05;
        ref_tmp = ref;
        ref_tmp(ref_tmp<thresh) = 0;
        BW = imbinarize(ref_tmp, thresh);
        s = regionprops(BW,ref,'WeightedCentroid');
        circN = size(s,1);
%         fprintf('circN = %d thresh = %d\n', circN, thresh);
    
    end
    tmp_s = struct2cell(tmp_s);
    s = fliplr(reshape(cell2mat(tmp_s), 2, Grid)');
    for ii=1:n(2)
        tmp_col = s(1+(ii-1)*n(1):ii*n(1),:);
        [~,I] = sort(tmp_col(:,1),"descend");
        s(1+(ii-1)*n(1):ii*n(1),:) = tmp_col(I,:);
    end
    cen = cat(3,flipud(reshape(s(:,1),n)),flipud(reshape(s(:,2),n))); % 1st dimention - row coordinates, 2nd dimention - col coordinates.
    
    % find mean distances between probes - "b", and the tilt angle of the image.
    if any(size(cen(:,:,1))>[1 1])
        if size(cen(:,:,1),1)>1
            cen_diff_rows = cen(2:end,:,:)-cen(1:end-1,:,:);
        else
            cen_diff_rows = zeros(1,1,2); cen_diff_rows(1,1,2) = eps;% in case only 1 column
        end
        if size(cen(:,:,1),2)>1
            cen_diff_cols = cen(:,2:end,:)-cen(:,1:end-1,:);
        else
            cen_diff_cols = zeros(1,1,2); cen_diff_cols(1,1,2) = eps; % in case only 1 row
        end
        b_rows = mean(sqrt(sum(cen_diff_rows.^2,3)),'all');
        b_cols = mean(sqrt(sum(cen_diff_cols.^2,3)),'all');
        angle_rows = atand(cen_diff_rows(:,:,2)./cen_diff_rows(:,:,1)); %frox vertical axis
        angle_cols = atand(cen_diff_cols(:,:,1)./cen_diff_cols(:,:,2)); %from horizontal axis
        angle_tilt = sign(mean(angle_rows(:)))*mean([abs(mean(angle_rows(:))),abs(mean(angle_cols(:)))]); %from vertical axis
%         if abs(abs(angle)-abs(mean(angle_rows(:))))/mean([abs(angle),abs(mean(angle_rows(:)))])>0.01 || abs(abs(angle)-abs(mean(angle_cols(:))))/mean([abs(angle),abs(mean(angle_cols(:)))])>0.01
%             error('your data is sqewed, please reacquire the data'); % more than 1% of sqew 
%         end
        b = [b_rows b_cols];
    end
    
    % find subimages size
    if exist('angle_tilt','var')
        subimage = [round(mean(cen_diff_rows(:,:,1),'all')),round(mean(cen_diff_cols(:,:,2),'all'))];
        subimage = subimage-~mod(subimage,2); % make subimage odd numbers
    else % if only one probe - default values
        subimage = size(im); 
        angle_tilt = 0;
        b = [0 0];
    end
end