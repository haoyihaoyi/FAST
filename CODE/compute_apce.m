function [APCE] = Compute_APCE(response_map,w,h)
    Fmax = max(max(response_map));
    Fmin = min(min(response_map));
    Diff_avg = sum(sum((response_map-Fmin).^2))/(w*h);
    APCE = ((Fmax-Fmin).^2)./Diff_avg;
end