function [p, v] = compute_prox(obj, x, gam)
    normx = norm(x);
    p = obj.R/normx * x;
    v = 0.0;
end
