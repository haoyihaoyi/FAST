function x = get_features_mhy(patch, features, cell_size, cos_window, w2c, w, h, channel)
        
    x = double(fhog(single(patch) / 255, cell_size, features.hog_orientations));
    
%     im_patch = imresize(patch, [w h]);
    im_patch = mexResize(patch, [w,h], 'auto');

    if channel > 1
        x(:, :, end) = (single(rgb2gray(im_patch))/255) - 0.5;
    else
        x(:, :, end) = (single(im_patch)/255) - 0.5;
    end
    x = bsxfun(@times, x, cos_window);
end