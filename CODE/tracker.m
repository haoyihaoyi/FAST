function [rect_results, time] = tracker(video_path, img_files, pos, target_sz, ...
	padding, lambda, output_sigma_factor, interp_factor, cell_size, ...
	search_size_3 ,search_size_4, features, show_visualization)

	resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
    if resize_image
		pos = floor(pos / 2);
		target_sz = floor(target_sz / 2);
    end
    
	window_sz = floor(target_sz * (1 + padding));
	output_sigma = sqrt(prod(target_sz)) * output_sigma_factor / cell_size;
	yf = fft2(gaussian_shaped_labels(output_sigma, floor(window_sz / cell_size)));
	cos_window = hann(size(yf,1)) * hann(size(yf,2))';	
	
    if show_visualization  %create video interface
		update_visualization = show_video(img_files, video_path, resize_image);
	end
	
	time = 0;  %to calculate FPS
	positions = zeros(numel(img_files), 2);  %to calculate precision
	rect_results = zeros(numel(img_files), 4);  %to calculate 
    search_size = cat(2,search_size_3,search_size_4);
    [w,h] = size(cos_window);
    response = zeros(w,h,size(search_size,2));

    im = imread([video_path img_files{1}]);
    [~,~,channel] = size(im);
    
    
    temp = load('w2crs');
    w2c = temp.w2crs;
%     w2c = [];
    szid = 1;

	for frame = 1:numel(img_files)
		im = imread([video_path img_files{frame}]);

        tic()
		if resize_image
			im = imresize(im, 0.5);
		end

		if frame > 1
            for i=1:size(search_size, 2)
                tmp_sz = floor((target_sz * (1 + padding))*search_size(i));
                patch = get_subwindow_mhy(im, pos, tmp_sz, window_sz);
                zf = fft2(get_features(patch, features, cell_size, cos_window, w2c, w, h, channel));
                kzf = linear_correlation(zf, model_xf);
                response(:,:,i) = real(ifft2(model_alphaf .* kzf));  %equation for fast detection
            end
            
            APCE = Compute_APCE(response,w,h);
            szid = find(APCE == max(APCE(:)));
            szid = szid(1);
            response_sel = response(:,:,szid);
            [vert_delta, horiz_delta] = find(response_sel == max(max(response_sel))); 
            vert_delta  = vert_delta(1);
            horiz_delta = horiz_delta(1);
			if vert_delta > w / 2  %wrap around to negative half-space of vertical axis
				vert_delta = vert_delta - w;
            end
            
			if horiz_delta > h / 2  %same for horizontal axis
				horiz_delta = horiz_delta - h;
            end

            tmp_sz = floor((target_sz * (1 + padding))*search_size(szid));
            current_size = tmp_sz(2)/window_sz(2);
			pos = pos + current_size*cell_size * [vert_delta - 1, horiz_delta - 1];
		end

        target_sz = target_sz * search_size(szid);
        tmp_sz = floor((target_sz * (1 + padding)));
        patch = get_subwindow_mhy(im, pos, tmp_sz, window_sz);
		xf = fft2(get_features(patch, features, cell_size, cos_window,w2c, w, h, channel));
        
        kf = linear_correlation(xf, xf);
        alphaf = yf ./ (kf + lambda);
        
        if frame == 1
            model_alphaf = alphaf;
            model_xf = xf;
        else
			model_alphaf = (1 - interp_factor) * model_alphaf + interp_factor * alphaf;
			model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
        end


		%save position and timing
		positions(frame,:) = pos;
		time = time + toc();

		box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        rect_results(frame,:)=box;
   		%visualization
%         savepath_name = '/home/haoyima/makevideo/';
%         figname = num2str(frame,'%03i.png');
%         final_path_name = [savepath_name figname];
	    if show_visualization
                stop = update_visualization(frame, box);
%                 saveas(gcf,final_path_name);
			if stop, break, end  %user pressed Esc, stop early
			drawnow
% 			pause(0.05)  %uncomment to run slower
		end
		
	end

	if resize_image
		positions = positions * 2;
        rect_results = rect_results*2;
    end
end

