function [rect_results, time] = tracker_good_mosse(video_path, img_files, pos, target_sz, ...
	padding, lambda, output_sigma_factor, interp_factor, cell_size, ...
	search_size_3 ,search_size_4, features, show_visualization)

	resize_image = (sqrt(prod(target_sz)) >= 80);  %diagonal size >= threshold
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
%     search_size = cat(2,search_size_3,search_size_4);
    search_size = 1;
    [w,h] = size(cos_window);
    response = zeros(w,h,size(search_size,2));

    im = imread([video_path img_files{1}]);
    [~,~,channel] = size(im);
    
    
%     temp = load('w2crs');
%     w2c = temp.w2crs;
    w2c = [];
    ori_dim = 32;
    rd_dim = 18;
    szid = 1;


	for frame = 1:numel(img_files)

		im = imread([video_path img_files{frame}]);
        tic()
        if resize_image
            im = imresize(im, 0.5);
            %             im = mexResize(im, 0.5*[w h], 'auto');
        end

		if frame > 1
            for i=1:size(search_size, 2)
                tmp_sz = floor((target_sz * (1 + padding))*search_size(i));
                patch = get_subwindow_mhy(im, pos, tmp_sz, window_sz);
                z = get_features_mhy(patch, features, cell_size, cos_window, w2c, w, h, channel);
                zf = fft2(feature_projection_mhy(reshape(z, [w*h, ori_dim]), projection_matrix, rd_dim, w, h));
%                 zf = fft2(z);
                response(:,:,i) = real(ifft2(sum(hf_num .* zf, 3) ./ (hf_den + lambda)));
                 %equation for fast detection
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
		x_pca = get_features_mhy(patch, features, cell_size, cos_window,w2c, w, h, channel);
        
        if frame == 1
            h_num_pca = x_pca;
        else
            h_num_pca = (1 - interp_factor) * h_num_pca + interp_factor * x_pca;
        end
        
        data_matrix = reshape(h_num_pca,[w*h, ori_dim]);
        [pca_basis, ~, ~] = svd(data_matrix' * data_matrix);
        projection_matrix = pca_basis(:, 1:rd_dim);
        
        hf_proj = fft2(feature_projection_mhy(data_matrix, projection_matrix, rd_dim, w, h));
        hf_num = bsxfun(@times, yf, conj(hf_proj));
        xlf = fft2(feature_projection_mhy(reshape(x_pca,[w*h, ori_dim]), projection_matrix, rd_dim, w, h));
        new_hf_den = sum(xlf .* conj(xlf), 3);

%         hf_proj = fft2(h_num_pca);
%         hf_num = bsxfun(@times, yf, conj(hf_proj));
%         xlf = fft2(x_pca);
%         new_hf_den = sum(xlf .* conj(xlf), 3);
%        
        
        
    
        if frame == 1
            hf_den = new_hf_den;
        else
            hf_den = (1 - interp_factor) * hf_den + interp_factor * new_hf_den;
        end
        

		%save position and timing
		positions(frame,:) = pos;
		time = time + toc();

		box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        rect_results(frame,:)=box;
%    		%visualization
%         savepath_name = '/home/haoyima/makevideo/';
%         figname = num2str(frame,'%03i.png');
%         final_path_name = [savepath_name [seq_name '/'] figname]; 
%         mkdir(final_path_name);
% 	    if show_visualization
%                 stop = update_visualization(frame, box);
%                 saveas(gcf,final_path_name);
% 			if stop, break, end  %user pressed Esc, stop early
% 			drawnow
% % 			pause(0.05)  %uncomment to run slower
% 		end
		
	end

	if resize_image
		positions = positions * 2;
        rect_results = rect_results*2;
    end
end

