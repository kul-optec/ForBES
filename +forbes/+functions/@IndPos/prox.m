function [p, v] = prox(obj, x, gam)
    p = max(obj.lo, x);
    v = 0;
end
