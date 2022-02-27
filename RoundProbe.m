classdef RoundProbe < Probe
    properties
        radius
    end
     properties (Dependent)
      probe
    end
    methods
        function obj = RoundProbe(r,s,c) % Constructor
            switch nargin
                case 0
                    r = 10; % default
                    s = [256,256]; % default
                    c  = [128,128]; % default
                case 1
                    s = [256,256]; % default
                    c  = [128,128]; % default
                case 2
                     c = [128,128]; % default
                otherwise
                    if nargin > 3
                        error('Too many input parameters!');
                    end
            end
            if r<0
                error('Radius must be nonegative');
            end
            obj@Probe(s,c);
            obj.radius = r;
        end
        
        % get methods
        function ind = get.probe(obj)
          [x,y] = meshgrid(1:obj.im_size(2), 1:obj.im_size(1));
          ind = double((x-obj.center(1)).^2+(y-obj.center(2)).^2 < obj.radius.^2);
        end  
        
        % set methods
        function obj = set.radius(obj,r)
            if r < 0
              error('Radius must be nonegative')
            end
            obj.radius = r;
        end

        % end of get and set methods
        
        
        function p = create_next_probe_horz(obj,overlap, err)
            c = [obj.center(1)+obj.find_next_center(overlap+err),obj.center(2)];
            p = Probe(obj.radius,round(c),obj.im_size);
        end
        function p = create_next_probe_vert(obj,overlap, err)
            c = [obj.center(1),obj.center(2)+obj.find_next_center(overlap+err)];
            p = Probe(obj.radius,round(c),obj.im_size);
        end
        
        
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

            fun = @(r,c) (2*acos(max(min(c,2*r),0)./(2*r))*r.^2-(max(c,0)/2).*sqrt(max(4*r.^2-c.^2,0)))./(pi*r.^2);
            cen_fun = @(c) fun(obj.radius,c)-overlap;
            cen = fzero(cen_fun,0.1);

        end
        
        
        function per = find_prob_overlap(obj1,obj2)
        % find overlap index between binary images ind1 & ind2, return
        % ratio of overlapping pixels (between 0 to 1) of the probes
            per = sum(sum(ismember(find(obj1.probe),find(obj2.probe))))/sum(obj1.probe(:));
        end
        
        function center = find_first_center(obj, overlap, n)
        %find the first probe center, according to the desired overlap area and number of probes at each axis 
            scan_area = round([obj.find_next_center(overlap)*(n(1)-1), obj.find_next_center(overlap)*(n(2)-1)]);
            if any(scan_area > obj.im_size)
                error('Scan area is bigger than the image');
            end
            center = round((obj.im_size-fliplr(scan_area))/2);
        end
        
        % make methods
        function [probes,cen] = make_probes(obj,n,p,rf)
        % function to produce "n" probs  of image size "s" with prob radius "r" 
        % and "p" percent overlaps 
        % input:
        % s - size of the image ([1x2] int)
        % r - radius of the probs (double)
        % n - number of probes, row and col ([1x2] int)
        % p - percent of overlap between probes (double)
        % rf - randomize factor in probe overlaping plus minus in percent (double)
        % output:
        % matrix of probes ([s(1) x s(2) x n])
        
            probes = zeros(obj.im_size(1),obj.im_size(2),n(1),n(2));
            cen = zeros(n(1),n(2),2);
            obj.center = obj.find_first_center(p,n);
            cen_dist = obj.find_next_center(p);
            [x,y] = meshgrid(obj.center(1):cen_dist:obj.center(1)+(n(2)-1)*cen_dist,obj.center(2):cen_dist:obj.center(2)+(n(1)-1)*cen_dist);
            cen = cat(3,x,y)+rand(size(cen))*rf/(2*sqrt(2));
            for ii=1:n(1)
                for jj=1:n(2)
                    obj.center = cen(ii,jj,:);
                    probes(:,:,ii,jj) = obj.probe;
                end
            end
            cen = flip(cen,3);
        end
    end
end
