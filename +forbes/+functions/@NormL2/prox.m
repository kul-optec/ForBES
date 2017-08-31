function [p, v] = prox(obj, x, gam)
    normx = norm(x, 'fro');
    lamgam = obj.lam*gam;
    if normx <= lamgam
        p = zeros(size(x));
        v = 0;
    else
        scal = (1-lamgam/normx);
        p = (1-lamgam/normx)*x;
        v = obj.lam*scal*normx;
    end
end