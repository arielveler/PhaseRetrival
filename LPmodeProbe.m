classdef LPmodeProbe < RealProbe
    properties
      order
      fiber_params % - CoreRad, CoreClad, lambda, n1, n2, v
    end
    properties (Dependent)
      probe
    end
    methods
        function obj = LPmodeProbe(order,params,pixel_pitch,im_size,center) % Constructor
            addpath('C:\Users\arielv\Desktop\Master\Fibers'); % for UC and params
            switch nargin
                case 0
                    order = DefaultOrder(); % LP_01 mode
                    params = DefaultParams(); % default - corning SMF-28 core radius
                    pixel_pitch = [3.45,3.45].*1e-6; % default - acA2440-35um - Basler ace
                    im_size = [2048,2048]; % default
                    center  = [1024,1024]; % default
                case 1
                    if class(order) == 'double'
                        order = DefaultOrder(order);
                    end
                    params = DefaultParams(); %default - corning SMF-28
                    pixel_pitch = [3.45,3.45].*1e-6; % default - acA2440-35um - Basler ace
                    im_size = [2048,2048]; % default
                    center  = [1024,1024]; % default
                case 2
                    pixel_pitch = [3.45,3.45].*1e-6; % default - acA2440-35um - Basler ace
                    im_size = [2048,2048]; % default
                    center  = [1024,1024]; % default
                case 3
                    im_size = [2048,2048]; % default
                    center  = [1024,1024]; % default
                case 4
                    center  = [1024,1024]; % default
                otherwise
                    if nargin > 5
                        error('Too many input parameters!');
                    end
            end
            if order.l < 0 || order.m < 1 || any(mod([order.l order.m],1) > 0)
              error('Order of the mode (l,m) must be a nonegative integer (l>=0, m>=1)');
            end
            if ~(order.polarization == 'x' || order.polarization == 'y')
              error('Polarization of the mode may be only "x" or "y" ');
            end
            if params.CoreRad < 0 || params.CladRad < 0
              error('radius must be nonegative');
            end
            if params.lambda < 210e-9 || params.lambda >10e-6
                error('central wavelength must be between 0.21-10 microns');
            end
            if params.n1 < 1 || params.n2 < 1 
                error('Refractive indexes must be larger than 1');
            end
            % special case: order: [0, m>1] & polarization 'y' is 0:
            if order.l == 0 && order.m > 1 && order.polarization == 'y'
                warning('radial symmetric modes cannot have "y" polarization, changing to "x" polarization');
                order.polarization = 'x';
            end
            
            obj@RealProbe(pixel_pitch,params.lambda,im_size,center);
            obj.order = order;
            obj.fiber_params = params;
            
        end
        
        % get methods
        function ind = get.probe(obj)
            
          %%%%%%%%%%%%%%%%%% LP modes %%%%%%%%%%%%%%%
          % The code is taken mainly from the article "Weakly guided fibers" by D. Glodge,
          % APPLIED OPTICS / Vol. 10, No. 10 / October 1971
          % And from "Optical Electronics in Modern Communications 4th edition" by Amnon Yariv
          epsilon = eps*1e11;
          cam_size = obj.im_size.*fliplr(obj.pix);
          l = obj.order.l; m = obj.order.m; rc = obj.fiber_params.CoreRad;
          NA = sqrt(obj.fiber_params.n1.^2+obj.fiber_params.n1.^2);
          v = obj.fiber_params.CoreRad.*(2*pi./obj.fiber_params.lambda).*NA;
          UC = load('UC.mat'); UC = UC.UC;% loading zeros of bessele(l-1,v)
          uc = UC{l+1};

          % V_cutoff = @(m,l)m*pi+(l-3/2)*(pi/2);         
          if m*pi+(l-3/2)*(pi/2) > v
            error(['This mode - LP{l=',num2str(l),', m=',num2str(m),'} is not possible in this fiber'])
          end
          
          w = sqrt(v^2-uc(m)^2);
          k = 1-(w^2+l^2+1)^(-1/2);
          s = sqrt(uc(m)^2-l^2-1);
          if l==0 && m==1
            u = (1+sqrt(2))*v/(1+(4+v^4)^0.25);
          else
            u = uc(m)*exp((asec(v/s)-asec(uc(m)/s))/s);
          end
          
          dx = obj.pix(2); dy = obj.pix(1);
          [x,y] = meshgrid(0:dx:dx*(obj.im_size(2)-1),0:dy:dy*(obj.im_size(1)-1));
          
          r = sqrt((x-cam_size(2)/2).^2+(y-cam_size(1)/2).^2);
          phi = atan((y-cam_size(1)/2)./(x-cam_size(2)/2)); phi(1,1) = 0;
          A = 1; % amplitude
          
          switch obj.order.polarization
              case 'x'
                  % x-polarized mode:
                  ind = (double(r.^2 < rc.^2).*A.*besselj(l,u.*r./rc) + ...
                         double(r.^2 > rc.^2).*A.*(besselj(l,u)./(epsilon+besselk(l,w))).*besselk(l,w.*r./rc)).*cos(l.*phi);
              case 'y'
                  % y-polarized mode:
                  ind = (double(r.^2 < rc.^2).*A.*besselj(l,u.*r./rc) + ...
                         double(r.^2 > rc.^2).*A.*(besselj(l,u)./(epsilon+besselk(l,w))).*besselk(l,w.*r./rc)).*sin(l.*phi);
          end
          if find(isnan(ind))
              ind = inpaint_nans(ind);
          end
          ind = obj.shift_center_fft(ind,obj.center);
          ind = ind./max(ind(:));
%           ind(ind < 0.0001) = 0;   %%%%%%%%%%%%%%%%%%%%%%%%%
%           probe_support = ind>max(ind(:)).*0.5; % x_cut = FWHM
%           ind = ind.*probe_support;
%           ind = ind./max(ind(:));
        end  
        
        % set methods
        function obj = set.order(obj,order)
            switch class(order)
                case 'double'
                    new_order.l = order(1);
                    new_order.m = order(2);
                    new_order.polarization = 'x'; % default polarization
                case 'struct'
                    new_order = order;
                otherwise
                        error('wrong parameter "order" ')
            end
            order = new_order;
            if order.l < 0 || order.m < 1 || any(mod([order.l order.m],1) > 0)
              error('Order of the mode (l,m) must be a nonegative complete');
            end
            if ~(order.polarization == 'x' || order.polarization == 'y')
              error('Polarization of the mode may be only "x" or "y" ');
            end
            % special case: order: [0, m>1] & polarization 'y' is 0:
            if order.l == 0 && order.m > 1 && order.polarization == 'y'
              warning('radial symmetric modes cannot have "y" polarization, changing to "x" polarization');
              order.polarization = 'x';
            end
            obj.order = order;
        end
        
        function obj = set.fiber_params(obj,params)
            if params.CoreRad < 0 || params.CladRad < 0
              error('radius must be nonegative');
            end
            if params.lambda < 210e-9 || params.lambda >10e-6
                error('central wavelength must be between 0.21-10 microns');
            end
            if params.n1 < 1 || params.n2 < 1 
                error('Refractive indexes must be larger than 1');
            end
            
            if nargin == 0
                params = DefaultParams();
            end
            obj.fiber_params = params;
        end

        % end of get and set methods
        
        % find methods
        function cen = find_next_center(obj, overlap)
        % find distance between center of adjacent probs, based on the overlaping
        % area (percentage or part) and prob radius
        % input:
        % area - desired percentage of overlapping area
        % radius - the radius of the probs
        % output:
        % cen - dintance of next center
            if overlap >1          %if in percentage
                overlap = overlap/100;
            end
            
%             xc = linspace(0,1,obj.im_size(1)); yc = linspace(0,1,obj.im_size(2));           
%             fun = @(c) trapz(xc,trapz(yc,obj.probe.*conj(abs(obj.shift_center_fft(obj.probe,[0 c]))))).^2./(trapz(xc,trapz(yc, abs(obj.probe).^2, 1)) * trapz(xc,trapz(yc, abs(abs(obj.shift_center_fft(obj.probe,[0 c]))).^2, 1 ))) - overlap;
%             cen = abs(fzero(fun,0.00001));
            
%             fun = @(c) sum(min(abs(obj.probe),abs(obj.shift_center_fft(obj.probe,obj.center+[0 c]))),'all')/sum(obj.probe,'all')-overlap;
%             cen = abs(fzero(fun,0.00001));

             fun = @(r,c) (2*acos(max(min(c,2*r),0)./(2*r))*r.^2-(max(c,0)/2).*sqrt(max(4*r.^2-c.^2,0)))./(pi*r.^2);
             cen_fun = @(c) fun(obj.fiber_params.CoreRad,c)-overlap;
             cen = fzero(cen_fun,0.1)./obj.pix(2);   %%%%%%%% check if better way %%%%%%%

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

function order = DefaultOrder(ord)
    % default parameters for Multimode fiber - OM2-OM5
    if nargin == 0
        order.l = 0; % Basic mode -  LP01
        order.m = 1; % Basic mode -  LP01
    else
        order.l = ord(1);
        order.m = ord(2);
    end
    order.polarization = 'x'; % 'x' polarization as default
end

function params = DefaultParams()
    % default parameters for Multimode fiber - OM2-OM5
    params.CoreRad = 50e-6; % 50 microns core radius
    params.CladRad = 125e-6; % 50 microns clad radius
    params.lambda = 850e-9; % 850 nm central wavelength
    params.n1 = refrIndex('SILICA', params.lambda); % Silica is core material
    params.n2 = refrIndex('SILICACLAD', params.lambda); 
end
