function [p, v] = prox(obj, x, gam)
    uz = max(0.0, abs(x) - obj.lam*gam);
    p = sign(x).*uz;
    v = obj.lam*sum(uz(:));
end
