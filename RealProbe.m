classdef RealProbe < Probe
    properties
        pix
        lam0
    end
    methods
        function obj = RealProbe(pixel_pitch,lambda,s,c) % Constructor
%             addpath('C:\Users\arielv\Desktop\תואר שני\Speckle-fibers'); % for UC and params
            switch nargin
                case 0
                    pixel_pitch = [3.45,3.45].*1e-6; % default - acA2440-35um - Basler ace
                    lambda = 850e-9; % [m] default
                    s = [2048,2048]; % default
                    c  = [1024,1024]; % default
                case 1
                    lambda = 850e-9; % [m] default
                    s = [2048,2048]; % default
                    c  = [1024,1024]; % default
                case 2
                    s = [2048,2048]; % default
                    c  = [1024,1024]; % default
                case 3
                    c  = [1024,1024]; % default
                otherwise
                    if nargin > 4
                        error('Too many input parameters!');
                    end
            end
            if any(pixel_pitch < 0) || any(~isa(pixel_pitch,'double'))
              error('Pixel pitch must be a nonegative double');
            end
            if lambda < 210e-9 || lambda >10e-6
                error('central wavelength must be between 0.21-10 microns');
            end
            
            obj@Probe(s,c);
            obj.pix = pixel_pitch;
            obj.lam0 = lambda;
            
        end
        
        % set methods
        function obj = set.pix(obj,pixel_pitch)
            if any(pixel_pitch < 0) || any(~isa(pixel_pitch,'double'))
              error('Pixel pitch must be a nonegative double');
            end
            obj.pix = pixel_pitch;
        end
        
        function obj = set.lam0(obj,lambda)
            if lambda < 210e-9 || lambda >10e-6
                error('central wavelength must be between 0.21-10 microns');
            end
            obj.lam0 = lambda;
        end
        
        %create pinhole image
        function [cen,sub_size] = pinhols(obj,n,centers)
        % function to produce "cen" center of probs according to number of pinhols n(1)x n(2)
        % basically, divide the image into n(1)x n(2) sub-images, returning
        % the centers and sub-size
        % input:
        % n - number of pinhols, row and col ([1x2] int)
        % output:
        % cen - centers of probes ([n(1) x n(2) x 2])
        % sub_size - size of one image of one probe 
            if exist('centers', 'var') 
                cen = centers;
                cen_x = mean(centers(:,:,2),1); cen_y = mean(centers(:,:,1),2); 
                sub_size = [mean(cen_x(2:end)-cen_x(1:end-1)),mean(cen_y(2:end)-cen_y(1:end-1))];
            else
                sub_size = obj.im_size./n;
                [x,y] = meshgrid(sub_size(2)/2:sub_size(2):(n(2)-0.5)*sub_size(2), sub_size(1)/2:sub_size(1):(n(1)-0.5)*sub_size(1));
                cen = cat(3,y,x);
            end
        end
        
    end
end
