function sp = SymPointsSelection(im_orig)
%This function take as input the image and return a list of some pairs of
%symmetric points

    dx = [-1 0 1; -1 0 1; -1 0 1];   % Derivative masks
    dy = dx';

    Ix = conv2(im_orig, dx, 'same');      % Image derivatives
    Iy = conv2(im_orig, dy, 'same');

    % set the parameter for Gaussian convolution used in Harris Corner Detector
    SIGMA_gaussian=4;
    g = fspecial('gaussian',max(1,fix(3*SIGMA_gaussian)+1), SIGMA_gaussian);

    Ix2 = conv2(Ix.^2, g, 'same'); % Smoothed squared image derivatives
    Iy2 = conv2(Iy.^2, g, 'same');
    Ixy = conv2(Ix.*Iy, g, 'same');


    % cim = det(M) - k trace(M)^2.
    % cim is large means that the eigenvalues of matrix M is large
    k = 0.04;
    cim = (Ix2.*Iy2 - Ixy.^2) - k * (Ix2 + Iy2);

    T=mean(cim(:));
    CIM=cim;
    CIM(find(cim<T))=0;

    support=true(30);
    % compute maximum over a square neighbor of size 30 x 30
    maxima=ordfilt2(CIM,sum(support(:)),support);
    % determine the locations where the max over the neigh or 30 x 30 corresponds to the cim values
    [loc_i,loc_j]=find((cim==maxima).*(CIM>0));

    %Inside symmetricPoints.mat I have some pairs of symmetric points that have been
    %taken manually. I'll look for the nearest points exstracted with the
    %harris approach.
    load('symmetricPoints.mat'); %load x and y variales
    symP=[x y];
    eps= 200;
    ind=[];

    %I'm parsing all the corner found, and I'll save the ones very close to the points loaded
    for(i=1:size(symP,1))
        for(j=1:size(loc_i,1))
            point=symP(i,:);
            error = (point(2)-loc_i(j)).^2+(point(1)-loc_j(j)).^2;
            if(error<eps)
                ind=[j; ind];
            end
        end
    end

    %I save in the structure below the result points
    new_loc_i=[];
    new_loc_j=[];

    for (i = 1:size(ind,1))
        new_loc_i=[loc_i(ind(i)); new_loc_i];
        new_loc_j=[loc_j(ind(i)); new_loc_j];
    end


    sp = [new_loc_j, new_loc_i];

end
