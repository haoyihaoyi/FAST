function patch = get_subwindow_mhy(im, pos, tmp_sz, window_sz)
        
    xs = floor(pos(2)) + (1:tmp_sz(2)) - floor(tmp_sz(2)/2);
    ys = floor(pos(1)) + (1:tmp_sz(1)) - floor(tmp_sz(1)/2);
%% check for out-of-border value and set them to the values at the border
    xs(xs < 1) = 1;
    ys(ys < 1) = 1;
    xs(xs > size(im,2)) = size(im,2);
    ys(ys > size(im,1)) = size(im,1);
    im_patch = im(ys, xs, :);
%% resize the image to model size
    patch = mexResize(im_patch, window_sz, 'linear');
%     patch = imresize(im_patch, window_sz,'bilinear');
end