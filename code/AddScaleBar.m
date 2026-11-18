function img = AddScaleBar(varargin)
switch nargin
    case 1
        im=varargin{1};
        pix=1;
        scale = 0.2;
        method = 'BottomRight';
    case 2
        im=varargin{1};
        pix=varargin{2};
        scale = 0.2;
        method = 'BottomRight';
    case 3
        im=varargin{1};
        pix=varargin{2};
        scale = varargin{3};
        method='BottomRight';
    case 4
        im=varargin{1};
        pix=varargin{2};
        scale = varargin{3};
        method=varargin{4};
end

opt = {'m' 'mm' 'um' 'nm' 'pm' 'fm'};
mag = 3;
s = size(im);
a=floor(pix*size(im,2)*scale);
c=0;
while ~a
    c = c+1;
    a = floor(pix*size(im,2)*scale*10^c);
end
l = floor(a/10^c/pix);
w = floor(a/10^c/pix/20);
sc = opt{ceil(c/3)+1};
if ~(mod(c,3) == 0)
    a = a*10^(3-mod(c,3));
end
switch method
    case 'TopRight'
        indx = floor(0.05*s(1)):floor(0.05*s(1))+w; indy = s(2)-floor(0.05*s(2))-l:s(2)-floor(0.05*s(2));
    case 'TopLeft'
        indx = floor(0.05*s(1)):floor(0.05*s(1))+w; indy = floor(0.05*s(2)):floor(0.05*s(2))+l;
    case 'BottomRight'
        indx = s(1)-floor(0.05*s(1))-w:s(1)-floor(0.05*s(1)); indy = s(2)-floor(0.05*s(2))-l:s(2)-floor(0.05*s(2));
    case 'BottomLeft'
        indx = s(1)-floor(0.05*s(1))-w:s(1)-floor(0.05*s(1)); indy = floor(0.05*s(2)):floor(0.05*s(2))+l;
end
NBITS = @(X) nextpow2(double(X)+1);
maxbits = NBITS(max(im(:)));
im(indx,indy) = 2^maxbits;
num_im = zeros(s);
if any(s> 400)
    if strcmp(method,'TopRight') || strcmp(method,'TopLeft')
        position_text = [indy(1) indx(end)-numel(indx)];
    else
        position_text = [indy(1) indx(1)-numel(indx)-1.2*min(floor(numel(indx)*4),200)];
    end
    text_str = {[num2str(a) sc]};
    num_im = insertText(num_im,position_text,text_str,'FontSize',min(floor(numel(indx)*4),200),'TextColor','white','Font','Cambria Bold Italic','BoxColor',...
        'w','BoxOpacity',0);
    im = im + rgb2gray(num_im)*2^maxbits;
    im(im>2^maxbits) = 2^maxbits;
%     im(rgb2gray(num_im)>0) = 2^maxbits/2;
    else
    warning('Image is too small to display letters');
end
% position_line = [indy(1) indx(ceil(end/2)) indy(end) indx(ceil(end/2))]; 
% im = insertShape(im,'Line',position_line,'LineWidth',w,'Color','white','Opacity',1,'SmoothEdges',true);

img = im;
end