function x = get_features(patch, features, cell_size, cos_window, w, h, channel)
        
    x = double(fhog(single(patch) / 255, cell_size, features.hog_orientations));
    im_patch = mexResize(patch, [w,h], 'auto');
    
%     im_patch = imresize(patch, [w h]);
    if channel > 1
        x(:, :, end) = (single(rgb2gray(im_patch))/255) - 0.5;
    else
        x(:, :, end) = (single(im_patch)/255) - 0.5;
    end
    x = bsxfun(@times, x, cos_window);
end