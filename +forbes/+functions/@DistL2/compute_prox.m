function [y, v] = compute_prox(obj, x, gam)
    p = obj.ind.compute_prox(x, 1.0);
    d = norm(x-p, 'fro');
    gamlam = (gam*obj.lam);
    if gamlam < d
        gamlamd = gamlam/d;
        y = (1-gamlamd)*x + gamlamd*p;
        v = obj.lam*(d-gamlam);
    else
        y = p;
        v = 0.0;
    end
end
