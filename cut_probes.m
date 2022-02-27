function [probes, subimage] = cut_probes(image, cen, subim)
% function for cutting the probes on the camera image into set of probes
% for SSP.
% input:
% image - camera intensity pattern ([NrowxNcol] complex)
% cen - centers of probes ([n(1)x n(2)x 2] double)
% subim (optional) - desired subimages of the probes
% output:
% probes - a 4D array with probes intensity on camera plane content by its 
%          order ([subimage(1), subimage(2), n(1), n(2)] complex)
% subimage - size of subimage in pixels ([1x2] int)

    if nargin < 2 
        error('Not enought input arguments');
    end
    n = size(cen,[1 2]); % number of spots to cut in every axis ([1x2] double)

    if nargin == 3
        subimage = subim;
    else
        % find subimages size
        if n == [1,1]  % if only one probe - default values
            subimage = size(im);
        else % if more than one probe
            cen_diff_rows = cen(2:end,:,:)-cen(1:end-1,:,:);
            cen_diff_cols = cen(:,2:end,:)-cen(:,1:end-1,:);
            subimage = [round(mean(cen_diff_rows(:,:,1),'all')),round(mean(cen_diff_cols(:,:,2),'all'))];
            subimage = subimage-~mod(subimage,2); % make subimage odd numbers
        end
    end
    cen = round(cen);
    
    for i=1:n(1)
        for j=1:n(2)
            cur_cen = reshape(cen(i,j,:),[1 2]);
            sub_image = image(cur_cen(1)-ceil(subimage(1)/2):cur_cen(1)+floor(subimage(1)/2)-1,cur_cen(2)-ceil(subimage(1)/2):cur_cen(2)+floor(subimage(2)/2)-1);
            probes(:,:,i,j) = sub_image;
        end
    end


end