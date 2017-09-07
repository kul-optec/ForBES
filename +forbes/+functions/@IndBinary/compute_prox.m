function [p, v] = compute_prox(obj, x, gam)
    diff_low = abs(x - obj.low);
    diff_high = abs(x - obj.high);
    ind_low = diff_low <= diff_high;
    p = ind_low .* obj.low + ~ind_low .* obj.high;
    v = 0.0;
end
