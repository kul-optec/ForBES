function [p, v] = prox(obj, x, gam)
    proj_x = obj.ind.prox(x, 1.0);
    p = (x+gam*proj_x)/(1+gam);
    v = 0.5/(gam+1.0)*norm(proj_x-x)^2-0.5/gam*norm(p-x)^2;
end
