function im_ok = SetupImage(imRead)

    im_orig=rgb2gray(im2double(imRead));

    % MAKE THE PICTURE SQUARED
    sizeIm = size(im_orig);
    hor_size = sizeIm(2);
    ver_size = sizeIm(1);
    del = (hor_size-ver_size)/2 ;
    im_ok=im_orig(:, del : (hor_size-del)-1);
    
end



