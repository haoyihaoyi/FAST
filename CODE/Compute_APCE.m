function [APCE] = Compute_APCE(response_map,w,h)
    Fmax = max(max(response_map));
    Fmin = min(min(response_map));
    for i = 1:size(response_map,3)
        response_map(:,:,i) = response_map(:,:,i)-Fmin(i);
    end
    Diff_avg = sum(sum(response_map.^2))/(w*h);
    APCE = ((Fmax-Fmin).^2)./Diff_avg;
end