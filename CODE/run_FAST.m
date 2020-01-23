function results = run_FAST(seq, res_path, bSaveImage)
%  function results = run_SITUPTEST(seq, res_path, bSaveImage)   
	padding = 1.5;  %extra area surrounding the target
    lambda = 0.08;  %regularization
    interp_factor = 0.004;
	output_sigma_factor = 0.1;  %spatial bandwidth (proportional to target)
    search_size_3 = [1 0.995 1.005];
    search_size_4 = [0.985 0.99 1.01 1.015];
%     search_size_4 = [0.985 0.99 1.01 1.015, 1.01,1.01,1.01,1.01,1.01,1.01,1.01,];
%     search_size = [1 0.985 0.99 0.995 1.005 1.01 1.015];
%     lambda = 0.001;
%     interp_factor = 0.004;
    cell_size = 4;
    features.hog_orientations = 9;
    target_sz = seq.init_rect(1,[4,3]);
    pos = seq.init_rect(1,[2,1]) + floor(target_sz/2);
    img_files = seq.s_frames;
    video_path = [];
%     seq_name = seq.name;
    %call tracker function with all the relevant parameters
    addpath(genpath('./mexResize'));
    [rect_results, t]= tracker_good_mosse(video_path, img_files, pos, target_sz, ...
        padding, lambda, output_sigma_factor, interp_factor, cell_size, ...
        search_size_3, search_size_4, features, 0);

    %return results to benchmark, in a workspace variable
    fps = numel(seq.s_frames) / t;
    disp(['fps: ' num2str(fps)])
    disp(['time: ' num2str(t)])
    results.type = 'rect';
    results.res = rect_results;
    results.fps = fps;
    results.time = t;
end