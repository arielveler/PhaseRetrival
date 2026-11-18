function imout = noise_removal(imin, imdark, method, beta)
% noise removal according to the article - "Background noise removal in
% x-ray ptychography", March 10 2017 / Applied Optics.

    % default
    if ~exist('method','var');  method = 'minimization';  end
    if ~exist('beta','var');    beta = 1;                 end

    % Binary image:
    BW = true(size(imin));
%     margins = round(size(imin)*0.03);
    margins = round(size(imin)*0.03);
    BW(margins(1):end-margins(1),margins(2):end-margins(2)) = false;

    % Minimization method:
    alpha = sum(BW.*imin.*imdark,'all')/sum(BW.*imdark.*imdark,'all');
    imeff = imin-alpha*imdark;  
    imeff(imeff<0) = 0;
    
    % Tresholding methode:
    switch method
        case 'minimization'
            imout = imeff;
        case 'thresholding'
%             imout = imeff.*(imeff>(max(BW(:).*imeff(:))-2*std(imeff(BW),0,'all')*beta));
            imout = imeff.*(imeff>(max(BW(:).*imeff(:))*beta));
    end

end