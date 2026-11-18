function varargout=phplot(varargin)
%PHPLOT(FIELD)
%Plots the phase of FIELD in a continuous color scale (hue) and represents
%the normalized amplitude as brightness (r+g+b)*amplitude.
%PHPLOT(FIELD,AMP,FLAG)
%If AMP = 0 the amplitude is not plot
%If ICON = 1 the function creates a figure with a dial scale (from 0 to
%2*pi) and radial brightness (from 0 to one)
%A=PHPLOT(...) creates a 3D uint8 array that can be saved as an image with
%IMWRITE(A,'filename','fmt').
%Iacopo Mochi, Lawrence Berkeley National Laboratory 06/6/2010
switch nargin
    case 1
        field=varargin{1};
        Amp=1;
        icon=0;
        method = 'TopRight';
    case 2
        field=varargin{1};
        Amp=varargin{2};
        icon=0;
        method = 'TopRight';
    case 3
        field=varargin{1};
        Amp=varargin{2};
        icon=varargin{3};
        method = 'TopRight';
    case 4
        field=varargin{1};
        Amp=varargin{2};
        icon=varargin{3};
        method = varargin{4};
    case 0
        print('PHPLOT requires at least 1 input argument');
        return
end
phase=angle(field);
amplitude=abs(field); %amplitude(amplitude>20) = max(amplitude(:));
amplitude=amplitude/max(amplitude(:));
if Amp==0
    amplitude=ones(size(amplitude));
end
gain = 1; % gain for the intensity image
A=zeros(size(field,1),size(field,2),3);     %Declare RGB array
A(:,:,1)=gain*(sin(phase)+1).*amplitude;     %Red
A(:,:,2)=gain*(sin(phase+pi/2)+1).*amplitude;%Green
A(:,:,3)=gain*(-sin(phase)+1).*amplitude;    %Blue
s = size(A,[1 2]);
if icon==1
%     figure
    phase=-pi:2*pi/255:pi;
    r=0.5*(sin(phase)+1);     %Red
    g=0.5*(sin(phase+pi/2)+1);%Green
    b=0.5*(-sin(phase)+1);    %Blue
    
    warphase=[r(:)';g(:)';b(:)'];
    
    [x,y]=meshgrid(-1:2/255:1);
    a=(1i*y+x);

    colormap(warphase');
    a((x.^2+y.^2)>1)=nan;
    [~,a] = phplot(a,1,2,method);
    scale = 6;
    a = imresize(a,s/scale,"bicubic");
    sa = min(size(a));

    switch method
        case 'TopRight'
            indx = floor(size(a,1)/2):floor(size(a,1)/2)+size(a,1)-1; indy = floor((1-1/scale)*s(2))-floor(size(a,2)/2)+1:floor((1-1/scale)*s(2))-floor(size(a,2)/2)+size(a,2);
        case 'TopLeft'
            indx = floor(size(a,1)/2):floor(size(a,1)/2)+size(a,1)-1; indy = floor(size(a,2)/2):floor(size(a,1)/2)+size(a,2)-1;
        case 'BottomRight'
            indx = floor((1-1/scale)*s(1))-floor(size(a,1)/2)+1:floor((1-1/scale)*s(1))-floor(size(a,1)/2)+size(a,1); indy = floor((1-1/scale)*s(2))-floor(size(a,2)/2)+1:floor((1-1/scale)*s(2))-floor(size(a,2)/2)+size(a,2);
        case 'BottomLeft'
            indx = floor((1-1/scale)*s(1))-floor(size(a,1)/2)+1:floor((1-1/scale)*s(1))-floor(size(a,1)/2)+size(a,1); indy = floor(size(a,2)/2):floor(size(a,1)/2)+size(a,2)-1;
    end

%     a(floor(end/2),1:floor(end/2),:) = 0;
    A(indx,indy,:) = double(a)/255;
    
end

switch nargout
    case 2
        varargout{2}=uint8(A*255);
    otherwise
end

if icon == 2
    return
end
if icon ==1 
    if any(s> 400)
        position_text = [indy(1)-floor(numel(indy)/3.5) indx(floor(end/2))-floor(numel(indx)/2); indy(end)-floor(numel(indy)/32) indx(floor(end/2))-floor(numel(indx)/2)];
        text_str = {'0',native2unicode([207 128],'UTF-8')};
        A = insertText(A,position_text,text_str,'FontSize',min(floor(numel(indx)/3),200),'Font','Cambria Bold Italic','BoxColor',...
            'w','BoxOpacity',0,'TextColor','white');
        position_line = [indy(1)-floor(numel(indy)/10) indx(ceil(end/2)) indy(end)+floor(numel(indy)/10) indx(ceil(end/2)); indy(ceil(end/2)) indx(1)-floor(numel(indx)/10) indy(ceil(end/2)) indx(end)+floor(numel(indx)/10)];
        A = insertShape(A,'Line',position_line,'LineWidth',ceil(min(numel(indx),numel(indy))/80),'Color','white','Opacity',1,'SmoothEdges',true);
    else
        warning('Image is too small to display letters');
    end
end
if nargout >= 1
    varargout{1} = image(A);
else
    image(A);
end
