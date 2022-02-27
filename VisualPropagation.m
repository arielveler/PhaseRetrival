function [Eout,imnumout] = VisualPropagation(type,params,zin,steps,imnumin)
fields = fieldnames(params);
imnumout = 0;
if imnumin
    if ~exist('.\movie','dir')
        mkdir '.\movie';
    end
end
for i=1:numel(fields)
    eval([fields{i},' = params.',fields{i},';']);
end
%check type of function

if strcmp(type,'Free') % free propagation - steps
    for i=1:steps
        Eout = FreePropagation(Ein,z*i/steps,lambda,k_x,k_y);
        imagesc(abs(Eout));impixelinfo;
        title(['Field at z = ',num2str(1e3*(zin+z*i/steps)),'[mm]']);
        drawnow;
        if imnumin %if save image
            imnumin = imnumin + 1;
            saveas(gcf,['.\movie\im_' ,num2str(imnumin) ,'.png']);
            clf
        end
    end
    imnumout = imnumin;

else % no distant propagation
    switch type
        case 'Lens'
            Eout = LensPropagation(Ein,f,lambda,x,y);

        case 'Multispot'
            Eout = MultispotPropagation(Ein,theta,n,lambda,x,y);

        case 'MicroLens'
            Eout = MicroLensPropagation(Ein,n,f,lambda,x,y);  
    end
    
    imagesc(abs(Eout));impixelinfo;
    title(['Field at z = ',num2str(1e3*zin),'[mm]']);
    drawnow;
    if imnumin %if save image
        imnumin = imnumin + 1;
        saveas(gcf,['.\movie\im_' ,num2str(imnumin) ,'.png']);
        clf
    end
    imnumout = imnumin;
end

end