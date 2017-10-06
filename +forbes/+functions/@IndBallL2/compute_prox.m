function [p, v] = compute_prox(obj, x, gam)
    normx = norm(x);
    if normx > obj.R
        p = obj.R/normx * x;
    else
        p = x;
    end
    v = 0.0;
end
