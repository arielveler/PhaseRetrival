function [cen,subimage,angle,b] = FindRefCenters(im,n,Thr,Dpix)
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
% angle - tilted angle (in degrees) of the data relative to the vertical axis
% b - calculated distance between adjacent probes vertically and
% horizontally seperatelly in subpixel manner ([1x2] double)

% created by Ariel Veler, Jan-2022

    % default values
    switch nargin
        case 1
            n = [5 5];
            Thr = 10;
            Dpix = 30;
        case 2
            Thr = 10;
            Dpix = 30;
        case 3
            Dpix = 30;
    end

   
    area_est = pi*(Dpix/2)^2;
    Grid = n(1)*n(2);
    GW = kron(gausswin(Dpix), gausswin(Dpix)'); % filter with a gaussian window
    ref = conv2(im, GW, 'same')/sum(GW(:));
    
    thresh = max(ref(:))/Thr;
    
    ref(ref<thresh) = 0;
    BW = imbinarize(ref, thresh);
    s = regionprops(BW,ref,'WeightedCentroid','Area');
    circN = size(s,1);
    if circN <Grid || circN==1
       error('Threshold is too high. try lower threshold');
    end
    flag =1 ;
    
    while circN>=Grid 
        if circN==Grid
            if flag
                tmp_s = s;
                flag = 0;
                continue
            else
            area_tmp = squeeze(cat(tmp_s(1).Area, tmp_s(:).Area));
            area_s = squeeze(cat(s(1).Area, s(:).Area));
            ind = abs(area_est-area_s)<abs(area_est-area_tmp);
            tmp_s(ind) = s(ind);
            end
        end
        thresh = thresh*1.05;
        ref(ref<thresh) = 0;
        BW = imbinarize(ref, thresh);
        s = regionprops(BW,ref,'WeightedCentroid','Area');
        circN = size(s,1);
%         fprintf('circN = %d thresh = %d\n', circN, thresh);
    
    end
    tmp_s = struct2cell(tmp_s);
    s = fliplr(reshape(cell2mat(tmp_s(2,:)), 2, Grid)');
    for ii=1:n(2)
        tmp_col = s(1+(ii-1)*n(1):ii*n(1),:);
        [~,I] = sort(tmp_col(:,1),"descend");
        s(1+(ii-1)*n(1):ii*n(1),:) = tmp_col(I,:);
    end
    cen = cat(3,flipud(reshape(s(:,1),n)),flipud(reshape(s(:,2),n))); % 1st dimention - row coordinates, 2nd dimention - col coordinates.
    
    % find mean distances between probes - "b", and the tilt angle of the image.
    if any(size(cen(:,:,1))>[1 1])
        cen_diff_rows = cen(2:end,:,:)-cen(1:end-1,:,:);
        cen_diff_cols = cen(:,2:end,:)-cen(:,1:end-1,:);
        b_rows = mean(sqrt(sum(cen_diff_rows.^2,3)),'all');
        b_cols = mean(sqrt(sum(cen_diff_cols.^2,3)),'all');
        angle_rows = atand(cen_diff_rows(:,:,2)./cen_diff_rows(:,:,1)); %frox vertical axis
        angle_cols = atand(cen_diff_cols(:,:,1)./cen_diff_cols(:,:,2)); %from horizontal axis
        angle = sign(mean(angle_rows(:)))*mean([abs(mean(angle_rows(:))),abs(mean(angle_cols(:)))]); %from vertical axis
%         if abs(abs(angle)-abs(mean(angle_rows(:))))/mean([abs(angle),abs(mean(angle_rows(:)))])>0.01 || abs(abs(angle)-abs(mean(angle_cols(:))))/mean([abs(angle),abs(mean(angle_cols(:)))])>0.01
%             error('your data is sqewed, please reacquire the data'); % more than 1% of sqew 
%         end
        b = [b_rows b_cols];
    end
    
    % find subimages size
    if exist('angle','var')
        subimage = [round(mean(cen_diff_rows(:,:,1),'all')),round(mean(cen_diff_cols(:,:,2),'all'))];
        subimage = subimage-~mod(subimage,2); % make subimage odd numbers
    else % if only one probe - default values
        subimage = size(im); 
        angle = 0;
        b = [0 0];
    end
end