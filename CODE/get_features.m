function out_pca = get_features(im, pos, tmp_sz, window_sz, cos_window, w, h, channel, ori_dim)

    xs = floor(pos(2)) + (1:tmp_sz(2)) - floor(tmp_sz(2)/2);
    ys = floor(pos(1)) + (1:tmp_sz(1)) - floor(tmp_sz(1)/2);
%% check for out-of-border value and set them to the values at the border
    xs(xs < 1) = 1;
    ys(ys < 1) = 1;
    xs(xs > size(im,2)) = size(im,2);
    ys(ys > size(im,1)) = size(im,1);
%% crop the image
    im_patch = im(ys, xs, :);
%% resize the image to model size
    %patch = mexResize(im_patch, window_sz, 'linear');
    patch = imresize(im_patch, window_sz,'bilinear');
%% get the features
    temp_pca = double(fhog(single(patch)/255, 4, 9));
%     temp_pca(:,:,32) = cell_grayscale(patch);
    im_patch = imresize(patch, [w h]);
    if channel > 1
        temp_pca(:, :, end) = (single(rgb2gray(im_patch))/255) - 0.5;
    else
        temp_pca(:, :, end) = (single(im_patch)/255) - 0.5;
    end
    temp_pca = bsxfun(@times, temp_pca, cos_window);
    out_pca = reshape(temp_pca, [w*h, ori_dim]);

end