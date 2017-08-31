function [p, v] = prox(obj, x, gam)
    over = abs(x) > sqrt(2*gam*obj.lam);
    p = x.*over;
    v = obj.lam*nnz(p);
end