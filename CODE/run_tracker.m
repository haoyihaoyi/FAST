%
%  High-Speed Tracking with Kernelized Correlation Filters
%
%  Joao F. Henriques, 2014
%  http://www.isr.uc.pt/~henriques/
%
%  Main interface for Kernelized/Dual Correlation Filters (KCF/DCF).
%  This function takes care of setting up parameters, loading video
%  information and computing precisions. For the actual tracking code,
%  check out the TRACKER function.

%  revised by: Haoyi Ma, Jan, 2020

function [rect_result, fps] = run_tracker(video, feature_type, show_visualization, show_plots)

	%path to the videos (you'll be able to choose one with the GUI).
	base_path ='/home/haoyi/DATA/TB100/';

	%default settings
	if nargin < 1, video = 'choose'; end
	if nargin < 2, feature_type = 'hog'; end
	if nargin < 3, show_visualization = ~strcmp(video, 'all'); end
	if nargin < 4, show_plots = ~strcmp(video, 'all'); end
	
	padding = 1.5;  %extra area surrounding the target
	lambda = 0.08;  %regularization
	output_sigma_factor = 0.1;  %spatial bandwidth (proportional to target)
	
    switch feature_type
        case 'hog'
            interp_factor = 0.004;
            features.hog = true;
            features.hog_orientations = 9;
            cell_size = 4;
        otherwise
            error('Unknown feature.')
    end

	switch video
	case 'choose'
		%ask the user for the video, then call self with that video name.
		video = choose_video(base_path);
        if ~isempty(video)
            [~, fps] = run_tracker(video, ...
                feature_type, show_visualization, show_plots);
            
            if nargout == 0  %don't output precision as an argument
                clear precision
            end
        end
	otherwise
		%we were given the name of a single video to process
		%get image file names, initial state, and ground truth for evaluation
		[img_files, pos, target_sz, ~, video_path] = load_video_info(base_path, video);
        search_size = [1 0.995 1.005 0.985 0.99 1.01 1.015];
		%call tracker function with all the relevant parameters
		[rects, time] = tracker(video_path, img_files, pos, target_sz, ...
            padding, lambda, output_sigma_factor, interp_factor, cell_size, ...
            search_size, features, show_visualization);
		rect_result = rects;
		fps = numel(img_files) / time;
        disp(['fps: ' num2str(fps)])
	end
end
