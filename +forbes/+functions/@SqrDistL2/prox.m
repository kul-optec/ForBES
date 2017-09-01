function [p, v] = prox(obj, x, gam)
    proj_x = obj.ind.prox(x, 1.0);
    lamgam = obj.ind.lam * gam;
    p = (x+lamgam*proj_x)/(1+lamgam);
    v = 0.5/(lamgam+1.0)*norm(proj_x-x, 'fro')^2-0.5/lamgam*norm(p-x, 'fro')^2;
end
