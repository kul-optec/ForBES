function [p, v] = compute_prox(obj, x, gam)
    a2 = obj.a^2;
    [p1, v] = obj.f.compute_prox(obj.a*x + obj.b, a2*gam);
    p = (p1 - obj.b)/obj.a;
end
