classdef GaussProbe < RealProbe
    properties
        sigma
    end
     properties (Dependent)
      probe
    end
    methods
        function obj = GaussProbe(sig,pixel_pitch,lambda,s,c) % Constructor
            switch nargin
                case 0
                    sig = 20.*3.45.*1e-6; % [m] default - 20 pixels
                    pixel_pitch = [3.45,3.45].*1e-6; % [m] default - acA2440-35um - Basler ace
                    lambda = 850e-9; % [m] default
                    s = [256,256]; % default
                    c  = s/2; % default
                case 1
                    pixel_pitch = [3.45,3.45].*1e-6; % [m] default - acA2440-35um - Basler ace
                    lambda = 850e-9; % [m] default
                    s = [256,256]; % default
                    c  = s/2; % default
                case 2
                    lambda = 850e-9; % [m] default
                    s = [256,256]; % default
                    c  = s/2; % default
                case 3
                    s = [256,256]; % default
                    c  = s/2; % default
                case 4
                    c = s/2; % default
                otherwise
                    if nargin > 5
                        error('Too many input parameters!');
                    end
            end
            if sig<0
                error('Standart deviation must be nonegative');
            end
            obj@RealProbe(pixel_pitch,lambda,s,c);
            obj.sigma = sig;
        end
        
        % get methods

%        function ind = get.probe(obj)
%           cam_size = obj.im_size.*fliplr(obj.pix); % [y,x]
%           dx = obj.pix(2); dy = obj.pix(1);
%           [x,y] = meshgrid(0:dx:cam_size(2)-dx,0:dy:cam_size(1)-dy);
%           ind = 1*exp(-((x-cam_size(2)/2).^2+(y-cam_size(1)/2).^2)./(2.*obj.sigma.^2));
%           ind = obj.shift_center_fft(ind,obj.center);
%           ind = ind./max(ind(:));
%           
%           probe_support = ind>max(ind(:)).*0.5; % x_cut = FWHM
%           ind = ind.*probe_support;
%           ind = ind./max(ind(:));
%         end 

        function ind = get.probe(obj)
          cam_size = obj.im_size.*fliplr(obj.pix); % [y,x]
          center = obj.center.*fliplr(obj.pix); % [y,x]
          dx = obj.pix(2); dy = obj.pix(1);
          [x,y] = meshgrid(0:dx:cam_size(2)-dx,0:dy:cam_size(1)-dy);
          ind = 1*exp(-((x-center(2)).^2+(y-center(1)).^2)./(2.*obj.sigma.^2));
          ind = ind./max(ind(:));
          
%           probe_support = ind>max(ind(:)).*0.5; % x_cut = FWHM
%           ind = ind.*probe_support;
%           ind = ind./max(ind(:));
        end 
        
        % set methods
        function obj = set.sigma(obj,sig)
            if sig < 0
              error('Standart deviation must be nonegative')
            end
            obj.sigma = sig;
        end

        % end of get and set methods
       
        % find methods
%         function cen = find_next_center(obj, overlap)
%         % find distance between center of adjacent probs, based on the overlaping
%         % area (percentage or part) and prob radius
%         % input:
%         % area - desired percentage of overlapping area
%         % radius - the radius of the probs
%         % output:
%         % cen - dintance of next center
%             if overlap >1          %if in percentage
%                 overlap = overlap/100;
%             end
% 
%             fun = @(hsize,sig,c) sum(min(abs(obj.shift_center_fft(fspecial('gaussian',hsize,sig),obj.center)),abs(obj.shift_center_fft(fspecial('gaussian',hsize,sig),obj.center+[0 c]))),'all')/sum(fspecial('gaussian',hsize,sig),'all');
%             cen_fun = @(c) fun(obj.im_size,obj.sigma./min(obj.pix),c)-overlap;
%             cen = abs(fzero(cen_fun,0.00001));
% 
%         end

        function dist = find_next_center(obj, overlap)
        % find distance between center of adjacent probs, based on the overlaping
        % area (percentage or part) and prob sigma
        % input:
        % overlap - desired percentage of overlapping area
        % output:
        % dist - dintance of next center
        
        % calculations: 
        % (sigma here ~ 1/e of signal)
        % overlap(dist,sig) = (1/sqrt(pi*sig^2))*integral(-inf,inf,min(exp(-((x-dist/2)/sig)^2),exp(-((x+dist/2)/sig)^2))dx)
        % overlap(dist,sig) = (2/sqrt(pi*sig^2))*integral(0,inf,exp(-((x+dist/2)/sig)^2)dx) = erfc(dist/(2*sig))
        % dist = erfcinv(overlap)*2*sig;
        
            if overlap >1          %if in percentage
                overlap = overlap/100;
            end
            
            dist = erfcinv(overlap)*2*(obj.sigma);

        end
        
        function per = find_prob_overlap(obj1,obj2)
        % find overlap index between binary images ind1 & ind2, return
        % ratio of overlapping pixels (between 0 to 1) of the probes
            per = sum(min(obj1.probe(:),obj2.probe(:)))/sum(obj1.probe(:));
        end
        
        function center = find_first_center(obj, overlap, n)
        %find the first probe center, according to the desired overlap area and number of probes at each axis 
            scan_area = round([obj.find_next_center(overlap)*(n(1)-1), obj.find_next_center(overlap)*(n(2)-1)]);
            if any(scan_area > obj.im_size)
                error('Scan area is bigger than the image');
            end
            center = round((obj.im_size-scan_area)/2);
        end
        
        % make methods
        function [probes,cen] = make_probes(obj,n,p,rf,cen_in)
        % function to produce "n" probs  of image size "s" with prob radius "r" 
        % and "p" percent overlaps 
        % input:
        % s - size of the image ([1x2] int)
        % r - radius of the probs (double)
        % n - number of probes, row and col ([1x2] int)
        % p - percent of overlap between probes (double)
        % rf - randomize factor in probe overlaping plus minus in percent (double)
        % cen_in (optional) - centers given to make probes accordingly ([n(1) x n(2) x 2] double)
        % output:
        % probes - matrix of probes ([s(1) x s(2) x n(1) x n(2)] double)
        % cen - centers of probes ([n(1) x n(2) x 2] double)
        
            probes = zeros(obj.im_size(1),obj.im_size(2),n(1),n(2));
            cen = zeros(n(1),n(2),2);
            obj.center = obj.find_first_center(p,n);
            cen_dist = obj.find_next_center(p);
            [x,y] = meshgrid(obj.center(2):cen_dist:obj.center(2)+(n(2)-1)*cen_dist,obj.center(1):cen_dist:obj.center(1)+(n(1)-1)*cen_dist);
            cen = cat(3,y,x)+rand(size(cen))*rf/(2*sqrt(2));
            
            if exist('cen_in', 'var')
                cen = cen_in;
            end
            
            for ii=1:n(1)
                for jj=1:n(2)
                  obj.center = reshape(cen(ii,jj,:),[1,2]);
                  probes(:,:,ii,jj) = obj.probe;
                end
            end
%             cen = flip(cen,3);
        end
        
        function [image,cen,sub_size] = pinhols(obj,n,centers)
        % function to produce "cen" center of probs according to number of pinhols n(1)x n(2)
        % basically, divide the image into n(1)x n(2) sub-images, returning
        % the centers and sub-size.
        % in addition, will return an image with the probes on it
        % input:
        % n - number of pinhols, row and col ([1x2] int)
        % output:
        % cen - centers of probes ([n(1) x n(2) x 2])
        % sub_size - size of one image of one probe
        % image - the probes on the pinhols grid ([obj.im_size])
            if exist('centers', 'var') 
                [cen,sub_size] = pinhols@RealProbe(obj,n,centers);
            else
                [cen,sub_size] = pinhols@RealProbe(obj,n);
            end
            image = zeros(obj.im_size);
            for ii=1:n(1)
                for jj=1:n(2)
                  obj.center = reshape(cen(ii,jj,:),[1,2]);
                  image = image + obj.probe;
                end
            end
        end
        
    end
    
    methods(Static)
        function image_out = shift_center_fft(image_in,c)
        % Image shift to center throuth fft and linear phase
            s = size(image_in);    
            sr = 0:s(1)-1;
            sc = 0:s(2)-1;
            [sc,sr] = meshgrid(sc,sr);
            image_out = ifft2(ifftshift(fftshift(fft2(image_in)).*exp(-1i*2*pi*((c(1)-s(1)/2)*sr/s(1)+(c(2)-s(2)/2)*sc/s(2)))));

        end
    end
end
