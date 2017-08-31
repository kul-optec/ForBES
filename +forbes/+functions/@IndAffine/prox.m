function [p, v] = prox(obj, x, gam)
    if isempty(obj.L_prox)
        obj.makeprox();
    end
    res = obj.A*x - obj.b;
    p = x - obj.A'*(obj.L_prox'\(obj.L_prox\res));
    v = 0.0;
end
