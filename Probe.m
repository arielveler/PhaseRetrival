classdef Probe
    properties
        im_size
        center
    end
    methods
        function obj = Probe(s,c) % Constructor
            switch nargin
                case 0
                    s = [256,256]; % default
                    c = [128,128]; % default
                case 1
                    c = [128,128]; % default             
                otherwise
                    if nargin > 2
                        error('Too many input parameters!');
                    end
            end
            if any(s<0)
                error('Image size indexes must be [1x2] positive integers')
            end
            if any(c<0)
                error('Center indexes must be [1x2] positive integers')
            end
            obj.im_size = s;
            obj.center = c;
        end
        
        % get methods

        % set methods     
        function obj = set.center(obj,c)
            if numel(c) == 2
                if c(1) >0 && c(2) > 0
                    obj.center = c;
                else
                    error('Center must be [1x2] positive double')
                end
            end
        end
        
        function obj = set.im_size(obj,s)
            if numel(s) == 2
                if s(1) >0 && s(2) > 0
                    obj.im_size = s;
                else
                    error('Image size must be [1x2] positive integers')
                end
            end
        end 
        % end of get and set methods

        % find methods
        
        % other methods
        
    end
end
        
